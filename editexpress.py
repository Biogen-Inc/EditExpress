#!/usr/bin/python
import sys
import re
import os
import argparse
import fnmatch
import subprocess
import string
if int(sys.version[0])==2:
    import ConfigParser
    from ConfigParser import ConfigParser as Parser
elif int(sys.version[0])==3:
    import configparser
    from configparser import ConfigParser as Parser
import shutil
import distutils.spawn
import errno
import itertools 
import multiprocessing as mp
import logging
import signal
import time
import six 
from six.moves import map,zip,range

def merge_mp_logs(path, workflow):
    mp_prefix = 'mp.'
    general_log = os.path.join(path, 'general.' + workflow + '_run.log')
    for files in os.listdir(path):
        with open(general_log, 'a') as f, open(os.path.join(path,files), 'r') as g:
            if fnmatch.fnmatch(files, mp_prefix+workflow+'*'):
                for line in g:
                    f.write(line)
    return

def merge_sp_logs(path, workflow):
    sp_prefix = 'sp.'
    general_log = os.path.join(path, 'general.' + workflow + '_run.log')
    list_of_logs = []
    for files in os.listdir(path):
        if fnmatch.fnmatch(files, sp_prefix+workflow+'*'):
            list_of_logs.append(files)
    if len(list_of_logs)>0:
        with open(general_log, 'w') as f:
            for i, files in enumerate(list_of_logs):
                with  open(os.path.join(path, files), 'r') as g:
                    for line in g:
                        if i==0:
                            f.write(line)
                        elif i>0 and not re.sub('\s+', ' ', line).split(' ')[4]=='Main':
                            f.write(line)                            
    return

def merge_summary_logs(path, workflow):
    summary_log = os.path.join(path, workflow+'_run.summary.log')
    general_log = os.path.join(path, 'general.' + workflow + '_run.log')
    list_of_log = []
    with open(general_log, 'r') as f:
        for line in f:
            if not re.sub('\s+', ' ', line).split(' ')[4]=='Summary':
                list_of_log.append(line)
    with open(general_log, 'w') as f, open(summary_log, 'r') as g:
        for i in list_of_log:
            f.write(i)
        for line in g:
            f.write(line)
    return


def natural_sortkey(string): 
    """short function that returns tuples of ints and strings to serve as key for 
    sorting numbers the 'human' way"""
    return tuple(int(num) if num else alpha for num, alpha in \
re.compile(r'(\d+)|(\D+)').findall(string))

def tableMerger(filespath, extension, output):
    """uses xlsxwriter module to concatenate multiple files as a multisheet excel 
    spreadsheet. function itself is actually fairly generic"""
    import xlsxwriter #lazy load  
    def num(s): 
        """attempts to coerce a string to int/float so that the cell value doesn't
        have any issues"""
        try:
            return int(s)
        except ValueError:
            try:
                return float(s)
            except ValueError:
                return(s)
    file_list = []
    for files in os.listdir(filespath):
        if fnmatch.fnmatch(files, '*'+extension):
            file_list.append(os.path.join(filespath, files))
    
    file_list2=sorted(file_list, key=natural_sortkey)
    sep = '\t'
    workbook  = xlsxwriter.Workbook(output)
    cell_format=workbook.add_format({'bold':True})
    sheet={}
    for files in file_list2:
        name = re.sub(extension, '', os.path.basename(files))
        sheet[name]=workbook.add_worksheet(name)
        out = open(files, 'r')
        i=0
        col_widths=[]
        for line in out:
            line = line.strip('\n')
            line = line.split(sep)
            if i == 0:
                for j in range(len(line)):
                    col_widths.append(0)
            for j in range(len(line)):
                temp=num(line[j])
                sheet[name].write(i, j, temp)
                if len(str(temp))+1>col_widths[j] and len(str(temp))+1<=30: #max width 30. keeps columns from becoming unwieldy (long sequences/mutation lists)
                    col_widths[j]=len(str(temp))+1
            i+=1

        for j in range(len(col_widths)):
            sheet[name].set_column(j,j, col_widths[j])
        
        sheet[name].set_row(0, None, cell_format)
        out.close()
    workbook.close()

def parseArgs():

    parser = argparse.ArgumentParser()
    parser.add_argument('--pipeline_config', required=True)
    parser.add_argument('--skip_program_check', action='store_true')
    parser.add_argument('--fastq1', default='')
    parser.add_argument('--fastq2', default='')
    args = vars(parser.parse_args())
    return args

def get_reference(ref_file): 
    """reads in the reference.txt file (e.g. target, primer(s), 
    amplicon, sgRNA). aims to be flexible (should work for 1 or 2 primers, multiple 
    sgRNAs). still requires specific column names and primer1 listed before primer2"""
    with open(ref_file) as f:
        clean1 = [re.sub('\s+', '\t', line.rstrip()) for line in f]
        clean2 = [re.sub('\t+', '\t', line) for line in clean1]
        clean3 = [line.split('\t') for line in clean2]
   
    if any([i for i in clean3 if len(i)!=max([len(i) for i in clean3])]): #one or more lines is of variable length
        print('Incorrectly formatted reference sequence - variable row lengths')
        sys.exit()
    y=len(clean3)
    x=len(clean3[0])

    target_indices=[i for i, elem in enumerate([str(j).lower() for j in clean3[0]]) if 'target' in elem]
    primer_indices=[i for i, elem in enumerate([str(j).lower() for j in clean3[0]]) if 'primer' in elem] 
    amplicon_indices=[i for i, elem in enumerate([str(j).lower() for j in clean3[0]]) if 'amplicon' in elem]
    sgRNA_indices=[i for i, elem in enumerate([str(j).lower() for j in clean3[0]]) if 'sgrna' in elem]

    primer1=[]
    primer2=[]
    amplicon=[]
    target=[]
    sgRNA=[]
    for i in range(y-1):
        if len(primer_indices)==1: #only forward primer
            primer1_parse=clean3[i+1][primer_indices[0]]
            primer2_parse=''
        elif len(primer_indices)==2: #forward and reverse primer
            primer1_parse=clean3[i+1][primer_indices[0]]
            primer2_parse=clean3[i+1][primer_indices[1]]
        elif len(primer_indices)!=1 and len(primer_indices)!=2: #no primers
            print("Incorrectly formatted reference file - primer error")
            sys.exit()
        if len(amplicon_indices)==1: #one amplicon sequence per line
            amplicon_parse=clean3[i+1][amplicon_indices[0]]
        elif len(amplicon_indices)!=1:
            print("Incorrectly formatted reference file - amplicon error")
            sys.exit()
        if len(target_indices)==1: #one target/genename per line
            target_parse=clean3[i+1][target_indices[0]]
        elif len(target_indices)!=1:
            print("Incorrectly formatted reference file - target error")
            sys.exit() 
        if len(sgRNA_indices)>=1: #one or more guides
            sgRNA_parse=[]
            for index in sgRNA_indices:
                sgRNA_parse.append(clean3[i+1][index])      
        elif len(sgRNA_indices)==0: #no guides
            sgRNA_parse=''

        primer1.append(primer1_parse)
        primer2.append(primer2_parse)
        amplicon.append(amplicon_parse)
        target.append(target_parse)
        sgRNA.append(sgRNA_parse)

    return(target, primer1, primer2, amplicon, sgRNA, y-1)

def reverse_complement(seq): #simple rc function. could have used biopython
    seq = seq.upper()
    alphabet={'A':'T', 'C':'G', 'T':'A', 'G':'C', 'N':'N'}
    seq_out=''
    for i in range(len(seq)):    
        seq_out+=alphabet[seq[i]]
    seq_out=seq_out[::-1]
    return seq_out

def fastq_setup(fastq_dir, R1_suffix, R2_suffix, paired_end, R1_in, R2_in, runmode, log):
    file_R1=[]
    file_R2=[]
    if runmode=='project':
        if paired_end:
            for file in os.listdir(fastq_dir):
                if fnmatch.fnmatch(file, '*'+R1_suffix):
                    file_R1.append(os.path.join(fastq_dir, file))
                elif fnmatch.fnmatch(file, '*'+R2_suffix):
                    file_R2.append(os.path.join(fastq_dir, file))
                else:
                    pass
            file_R1.sort()
            file_R2.sort()
            if len(file_R1)==0:
                log.critical("No fastqs with suffix " + R1_suffix + " found in " + fastq_dir + "! Exiting")
                sys.exit()
            if len(file_R2)==0:
                log.critical("No fastqs with suffix " + R2_suffix + " found in " + fastq_dir + "! Exiting")
                sys.exit()
            if len(file_R1) != len(file_R2):
                log.critical("Inequal number of R1/R2 fastq files! Exiting")
                sys.exit()

        elif not paired_end:
            for file in os.listdir(fastq_dir):
                if fnmatch.fnmatch(file, '*'+R1_suffix):
                    file_R1.append(os.path.join(fastq_dir, file))
                    file_R2.append("Empty")
                else:
                    pass
            if len(file_R1)==0:
                log.critical("No fastqs with suffix " + R1_suffix + " found in " + fastq_dir + "! Exiting")
                sys.exit()
    elif runmode=='single_sample':
        if R1_suffix not in R1_in:
            R1_in=R1_in+R1_suffix
        if os.path.isfile(R1_in):
            file_R1.append(R1_in)
        elif os.path.exists(os.path.join(fastq_dir, R1_in)):
            file_R1.append(os.path.join(fastq_dir, R1_in))
        else:
            log.critical(R1_in + ' not found as specified or in ' + fastq_dir)
        if paired_end:
            if R2_suffix not in R2_in:
                R2_in=R2_in+R2_suffix
            if os.path.isfile(R2_in):
                file_R2.append(R2_in)
            elif os.path.exists(os.path.join(fastq_dir, R2_in)):
                file_R2.append(os.path.join(fastq_dir, R2_in))
            else:
                log.critical(R2_in + ' not found as specified or in ' + fastq_dir)
        elif not paired_end:
            file_R2.append('Empty')
        if len(file_R1)==0 or len(file_R2)==0:
            log.critical('Missing fastqs! Exiting')
            sys.exit()
        
    return(file_R1, file_R2)

def create_amplicon_fasta(ref_file, paired_end, target, primer1, primer2, amplicon, path, log):
    """runs some checks and generates a fasta from the reference txt file, adds
    primers if necessary, determines if as given is actually reverse complement """
    rc='no' 
    amplicon_fasta = target+'.fasta' #fasta derived from reference file
    amplicon_fasta = os.path.join(path, amplicon_fasta)
    header=target
    write_fasta=False
    try:
        with open(amplicon_fasta): pass
    except IOError:
        write_fasta=True
    if primer1.upper() in amplicon.upper(): #is primer1 in the amplicon sequence? if so, the amplicon is primered. if not, then the primer sequences are added
        if paired_end:
            if len(primer2)>0 and primer2.upper() in amplicon.upper(): #primer2 isn't blank and is in the amplicon sequence
                pass
            elif len(primer2)>0 and reverse_complement(primer2) in amplicon.upper(): #primer2 isn't blank and the reverse complement is in the amplicon sequence
                rc='yes'
            else:
                logging.critical("Amplicon is primered but neither primer2 nor its reverse complement were found. Please verify sequence.")
                sys.exit()
        if write_fasta:
            logging.info('Creating ' + amplicon_fasta)
            f_fasta = open(amplicon_fasta, 'w')
            f_fasta.write('>'+header+'\n')
            f_fasta.write(amplicon+'\n')
            f_fasta.close()
    else:
        if write_fasta:
            logging.warning("Amplicon not primered. Primer(s) will be added to the amplicon sequence; if there are downstream issues, try taking the reverse complement of primer2")
            f_fasta = open(amplicon_fasta, 'w')
            f_fasta.write('>'+header+'\n')
            f_fasta.write(primer1+amplicon+primer2+'\n')
            f_fasta.close()
    return(amplicon_fasta, rc)

def default_config():
    """creates a dictionary of default values for most parameters"""
    dictionary = {}
    dictionary['i/o']={'output_directory':'./'}
    dictionary['qc_parameters']={'percent_assembled':50, 'percent_primered':50, 'percent_mapped_amplicon':50, 'percent_mapped_genome':50, 'percent_max_rl':50}
    dictionary['runtime']={'runmode':'project'}
    dictionary['modules']={'workflow':'standard','run_qc':True, 'alignment':True, 'mutation_calling':True}
    dictionary['processing']={'paired_end':True, 'read_merger':'flash', 'r1_suffix':'_L001_R1_001.fastq', 'r2_suffix':'_L001_R2_001.fastq', 'merge_mut_tables':False, 'clean_intermediates':False}
    dictionary['mutation_calling']={'restrict_alignment':False, 'bwa_min_mapq':15, 'summary_threshold':3.0, 'topseqs_min_hits':10, 'variant_filter':0.25, 'parse_substitutions':True}
    dictionary['offtarget']={'bwa_min_mapq':15, 'percent_cutoff':0, 'n_support':1}
    dictionary['multiprocess']={'n_processes':1}
    dictionary['programs']={'featurecounts':'featureCounts', 'samtools':'samtools', 'needle':'needle', 'bwa':'bwa', 'pandaseq':'pandaseq', 'seqtk':'seqtk', 'fastqc':'fastqc', 'bedtools':'bedtools', 'flash':'flash'}
    dictionary['needle']={'gapopen':10, 'gapextend':0.5}
    dictionary['bwa_mem']={'min_seed_length':19, 'band_width':100, 'z_dropoff':100, 'reseeding':1.5, 'mem_discard_threshold':10000, 'matching_score':1, 'mismatch_penalty':4, 'gap_open_penalty':6, 'gap_extension_penalty':6, 'clipping_penalty':5, 'unpaired_read_pair_penalty':9}
    dictionary['pandaseq']={'phred64':False, 'algorithm':'pear', 'no_uncalled_bases':False, 'min_overlap':5, 'max_assembly_length':499, 'min_assembly_length':120, 'threshold_score':0.6}
    dictionary['flash']={'min_overlap':5, 'max_overlap':65,'max_mismatch_density':0.25, 'allow_outies':False, 'phred_offset':33}
    return(dictionary)
    
def validate_config(parser_dict, default_dict):
    """runs some checks to see if parameters as given are of the correct type"""
    for section in default_dict:
        if section not in parser_dict:
            parser_dict[section]=default_dict[section]
        else:
            for param in default_dict[section]:
                if param not in parser_dict[section] or parser_dict[section][param]=="":
                    parser_dict[section][param]=default_dict[section][param]
                    print(('Missing value for ' + str(section) + ':' + str(param) + '. Using default value of ' + str(default_dict[section][param])))
                elif type(parser_dict[section][param])!=type(default_dict[section][param]):
                    error_msg = 'Value for ' + str(section) + ':'+str(param) + ' is of incorrect type, using the default value. Entered: ' + str(parser_dict[section][param]) + ' Default:' + str(default_dict[section][param])
                    if type(default_dict[section][param]) is float:
                        try:
                            parser_dict[section][param] = float(parser_dict[section][param])
                        except ValueError:
                            print(error_msg)
                            parser_dict[section][param]=default_dict[section][param]
                    elif type(default_dict[section][param]) is int:
                        try:
                            parser_dict[section][param] = int(parser_dict[section][param])
                        except ValueError:
                            print(error_msg)
                            parser_dict[section][param]=default_dict[section][param]
                    elif type(default_dict[section][param]) is bool:
                        if parser_dict[section][param].lower() == 'yes':
                            parser_dict[section][param]=True
                        elif parser_dict[section][param].lower() == 'no':
                            parser_dict[section][param]=False
                        else:
                            print(error_msg)
                            parser_dict[section][param]=default_dict[section][param]
    if 'fastq_directory' not in parser_dict['i/o'] and parser_dict['runtime']['runmode'].lower()=='project':
        print('No fastq directory given! Project level runs require this. Exiting')
        sys.exit()

    return(parser_dict)

def validate_programs(parser_dict):
    """checks if programs are valid and executable by the user, then compares that
    against what the user has decided to run, exiting if there are conflicts"""
    def is_tool(path):
        """is the path to a program a valid executable?"""
        return distutils.spawn.find_executable(path) is not None
    featureCounts= is_tool(parser_dict['programs']['featurecounts']) and os.access(parser_dict['programs']['featurecounts'], os.X_OK)
    samtools = is_tool(parser_dict['programs']['samtools']) and os.access(parser_dict['programs']['samtools'], os.X_OK)
    needle = is_tool(parser_dict['programs']['needle']) and os.access(parser_dict['programs']['needle'], os.X_OK)
    bwa = is_tool(parser_dict['programs']['bwa']) and os.access(parser_dict['programs']['bwa'], os.X_OK)
    seqtk = is_tool(parser_dict['programs']['seqtk']) and os.access(parser_dict['programs']['seqtk'], os.X_OK)
    pandaseq = is_tool(parser_dict['programs']['pandaseq']) and os.access(parser_dict['programs']['pandaseq'], os.X_OK)
    bedtools = is_tool(parser_dict['programs']['bedtools']) and os.access(parser_dict['programs']['bedtools'], os.X_OK)
    flash = is_tool(parser_dict['programs']['flash']) and os.access(parser_dict['programs']['flash'], os.X_OK)
    flag=False
    if not needle and parser_dict['modules']['alignment'] and parser_dict['modules']['workflow'].lower()=='standard':
        print('Invalid path to needle aligner. Required for the alignment module.')
        flag=True
    if not needle and parser_dict['modules']['run_qc']:
        print('Invalid path to needle aligner. Required for qc module.')
        flag=True
    if not seqtk and parser_dict['modules']['alignment'] and parser_dict['modules']['workflow'].lower()=='standard':
        print('Invalid path to the Seqtk. Required for alignment module')
        flag=True
    if not pandaseq and parser_dict['modules']['alignment'] and parser_dict['modules']['workflow'].lower()=='standard' and parser_dict['processing']['read_merger'].lower()=='pandaseq':
        print('Invalid path to the selected read merger, pandaseq. Required for alignment module')
        flag=True
    if not flash and parser_dict['modules']['alignment'] and parser_dict['modules']['workflow'].lower()=='standard' and parser_dict['processing']['read_merger'].lower()=='flash':
        print('Invalid path to the selected read merger, flash. Required for alignment module')
        flag=True
    if not bwa and parser_dict['modules']['alignment'] and parser_dict['modules']['workflow'].lower()=='gapped' or not bwa and parser_dict['modules']['alignment'] and parser_dict['modules']['workflow'].lower()=='offtarget':
        print('Invalid path to bwa aligner. Required for the gapped and offtarget alignment modules.')
        flag=True
    if not bwa and parser_dict['modules']['run_qc']:
        print('Invalid path to bwa aligner. Required for qc module.')
        flag=True
    if not samtools and parser_dict['modules']['alignment'] and parser_dict['modules']['workflow'].lower()=='gapped':
        print('Invalid path to samtools. Required for gapped alignment module.')
        flag=True
    if not samtools and parser_dict['modules']['run_qc']:
        print('Invalid path to samtools. Required for qc module.')
        flag=True         
    if not bedtools and parser_dict['modules']['workflow'].lower()=='offtarget' and parser_dict['modules']['mutation_calling']:
        print('Invalid path to bedtools. Bedtools  is required for mutation calling for the offtarget runtime.')
        flag=True
    if flag:
        print('Invalid executables, exiting. If you believe this is in error, try flag "--skip_program_checks"')
        sys.exit()
 
if __name__ == '__main__':
    
    def sample_level_processing(file_R1R2): 
        """ALL sample level steps contained in the same function (alignment, 
        mutation calling, qc, etc) for easy use of the map function. The object
        file_R1R2 is just a list containing all relevant sample level information. 
        At minimum it is just the paths to the fastqs. Inherits parser_dict"""
        file_R1, file_R2 = file_R1R2[0][0] #may not work with offtarget still 
        sample_name = re.sub(R1_suffix, '', os.path.basename(file_R1))
        sampleFilter=LogFilter()
        sampleFilter.sample_name=sample_name
        if file_R1R2[len(file_R1R2)-1]: #a value of True or False indicating multiprocessing
            mpLogger = logging.getLogger()
            mpLogger.addFilter(sampleFilter)
            logFormatter = logging.Formatter(fmt='%(asctime)s %(levelname)-8s %(app_name)8s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
            mpLogger = logging.getLogger()
            for handler in mpLogger.handlers:
                mpLogger.removeHandler(handler)
            mpLogger.setLevel('DEBUG')
            mp_log=os.path.join(log_folder,'mp.'+workflow+'_run.'+sample_name+'.log')
            fileHandler = logging.FileHandler(mp_log, mode='w')
            fileHandler.setFormatter(logFormatter)
            mpLogger.addHandler(fileHandler)
            log=mpLogger
        else:
            log=rootLogger 
            log.addFilter(sampleFilter)
        log.readout('Beginning sample level steps') 
        if workflow=='offtarget':
            offtarget_mut_dir = os.path.join(output_dir, 'offtarget', 'mut')
            mapq = parser_dict['mutation_calling']['bwa_min_mapq']
            bedtools = parser_dict['programs']['bedtools']
            percent_cutoff=parser_dict['offtarget']['percent_cutoff']
            n_support=parser_dict['offtarget']['n_support']
            try:
               os.makedirs(offtarget_mut_dir)
            except OSError as e:
                if e.errno != errno.EEXIST:
                    raise
            bam = os.path.join(genome_bam_dir, sample_name + '.genome.bam')
            sam = os.path.join(genome_bam_dir, sample_name + '.genome.sam')
            try:
                if not os.path.exists(bam) and not os.path.exists(sam):
                    from align_sample_bwa import distance, pre_process, get_params, bwa_align
                    fastqs=pre_process(file_R1, file_R2,paired_end)
                    bwa_mem_params = get_params(parser_dict)
                    log.readout('Aligning')
                    bwa_align(genome_bam_dir, sample_name, bwa_mem_params, genome_fasta, fastqs, bwa, samtools, workflow,log)
                    log.readout('Finished aligning')

                elif os.path.exists(sam):
                    with open(bam, 'wb') as f:
                        cmd = samtools + ' view -hSb ' + sam
                        subprocess.check_call(cmd.split(), shell=False, stdout=f)
                else:
                    pass    #pointless but I want to be explicit.
            
                from offtarget_indel_caller import getIndels, indelCaller
                log.readout('Beginning offtarget mutation analysis')
                indelCaller(bam,mapq,percent_cutoff, n_support, offtarget_mut_dir, samtools, bedtools)
                log.readout('Completed offtarget mutation analysis')
            except OSError:
                log.critical('Unable to process this sample')
                log.critical('Exiting')
                sys.exit()
        
        elif run_qc or workflow!='offtarget':
            mutation_calling = parser_dict['modules']['mutation_calling'] #want to be able to turn it off per sample
            ref_num, target,mut_dir, bam_dir, qc_dir, processed_fasta_dir, amplicon_fasta, primer1, primer2, sgRNA=file_R1R2[0][1]
            if alignment and workflow=='gapped' or run_qc and workflow=='gapped':
                from align_sample_bwa import distance, pre_process, get_params, bwa_align 
                fastqs=pre_process(file_R1, file_R2,paired_end)  
                bwa_mem_params = get_params(parser_dict)
                log.readout('Aligning')
                bwa_align(bam_dir, sample_name, bwa_mem_params, amplicon_fasta, fastqs, bwa, samtools, workflow, log)
                log.readout('Finished aligning')
            elif alignment and workflow=='standard' or run_qc and workflow=='standard':
                from align_sample_needle import distance, primered_reads, process_fastq, needle_align, pre_process_fastqs, get_param
                assembler_param, needle_param, assembler_path, needle, seqtk, read_merger= get_param(parser_dict)
                log.readout('Preprocessing')
                fasta_primered = pre_process_fastqs(file_R1, file_R2, output_dir, target, sample_name, seqtk, read_merger, assembler_path, assembler_param, primer1, primer2, rc, log)
                log.readout('Aligning')
                success=needle_align(bam_dir, sample_name, amplicon_fasta, fasta_primered, needle_param, needle, samtools,log)
                if not success:
                    mutation_calling = False
                    log.readout('Unable to complete alignment')
                else:
                    log.readout('Finished aligning')
            if run_qc:
                from sample_qc import main as main_qc, get_ref_and_runtime, count_read_lengths, distance, primer_rate_R1_R2
                gtf, bwa_mem_params = get_ref_and_runtime(parser_dict)
                log.readout('Running qc steps on ' + sample_name)
                main_qc(qc_dir, parser_dict, file_R1, file_R2, paired_end, ref_file, genome_fasta, genome_bam_dir, gtf, sample_name, bwa_mem_params, log)
                if paired_end and workflow=='standard':
                    primer_rate_R1_R2(parser_dict, amplicon_fasta, bam_dir, processed_fasta_dir, sample_name, file_R1, file_R2, primer1, primer2, rc, log)
            if mutation_calling:
                if workflow=='gapped':
                    try:
                        from parse_bwamem_alignment import get_mut_parser_params, single_sample_mut_parser, setup, parse_cigar, filter_reads, get_frame
                        mapq, parse_snps, noise_filter = get_mut_parser_params(parser_dict)
                        ref_dict = setup(amplicon_fasta)
                        log.readout('Calling mutations on')
                        bam = os.path.join(bam_dir, re.sub(R1_suffix, '.bam', os.path.basename(file_R1)))
                        single_sample_mut_parser(mut_dir, ref_dict, bam, mapq, parse_snps, noise_filter)
                        log.readout('Finished calling mutations')
                    except IOError:
                        logging.critical('bamfile missing!: ' + bam) 
                        return
                elif workflow=='standard':
                    try:
                        from parse_needle_alignment import setup_mut_caller,get_mut_parser_params, single_sample_mut_parser, parse_cigar, filter_reads, get_seq, get_frame, setup_restrict, parse_restrict_alignment 
                        q_cutoff, parse_snps, noise_filter = get_mut_parser_params(parser_dict)
                        try:
                            coords=parser_dict['alignment_plot']['seq_range']
                        except KeyError:
                            coords=''
                        ref_dict=setup_mut_caller(amplicon_fasta)
                        restrict_alignment, seq_range=setup_restrict(parser_dict, coords, ref_num, log)
                        bam = os.path.join(bam_dir, re.sub(R1_suffix, '.bam', os.path.basename(file_R1)))
                        log.readout('Calling mutations')
                        single_sample_mut_parser(mut_dir, q_cutoff, ref_dict, bam, parse_snps, noise_filter,restrict_alignment, seq_range)
                        log.readout('Finished calling mutations')
                    except IOError:
                        logging.critical('bamfile missing!: ' + bam) #this will happen if mutation calling is run separate from alignment, but the bams are missing
                        return
                    log.info('Building alignment plot...')
                    from build_plot_svg import visualize, parse_sample, parse_mut, setup_alignment_plot, pre_alignment_plot
                    zero_site, seq_range, ref_lines, default = setup_alignment_plot(parser_dict, amplicon_fasta, coords, ref_num, sgRNA, log)
                    pre_alignment_plot(mut_dir, os.path.join(mut_dir, sample_name)+mut_suffix, 'sample', parser_dict, seq_range, default, zero_site, ref_lines, target, log)

                from pie_charts import parse_sample_level, pie_chart
                logging.info('Baking some pie charts')
                parse_sample_level(os.path.join(mut_dir, sample_name)+mut_suffix)
            log.readout('Completed sample level steps')
            log.removeFilter(sampleFilter)
                       
    def summary_steps(per_ref_details):
        """All summary level processing steps are contained in this function.
        The object per_ref_details is just a list containing relevant directories,
        fasta, primer information, etc, for ease of applying the map function"""
        ref_num, target,mut_dir, bam_dir, qc_dir, processed_fasta_dir, amplicon_fasta, primer1, primer2, sgRNA = per_ref_details

        log = rootLogger
        summary = LogFilter()
        summary.sample_name='Summary'
        log.addFilter(summary)
        log.readout('Summarizing %s...', target)
        mutation_calling = parser_dict['modules']['mutation_calling']
        try:
            if not any(['.mut.xls' in i for i in os.listdir(mut_dir)]):
                mutation_calling = False
                log.error('Mutation calling was selected by mutation tables are missing!')
        except OSError:
            log.error('Mutation calling was selected but the mutation directory is missing!')
            mutation_calling = False

        if mutation_calling:
            from summary import setup_files_topSeqs, summary_topSeqs
            log.readout('Generating topSeqs file')
            setup_files_topSeqs(mut_dir, parser_dict, log)
            if merge_mut_tables:
                log.readout('Trying to merge mutation tables...')
                try:
                    tableMerger(mut_dir, '.mut.xls', mut_dir+'/merged_mut_calling.xlsx')
                    log.readout('Successfully merged mutation tables!')
                except ImportError:
                    log.warning('Missing xlsxwriter package! Please install if you want to generate a merged xlsx file')
            if workflow=='standard':
                try:
                    coords=parser_dict['alignment_plot']['seq_range']
                except KeyError:
                    coords=''
                from build_plot_svg import pre_alignment_plot, setup_alignment_plot, parse_mut, parse_topseqs, visualize
                zero_site, seq_range, ref_lines, default = setup_alignment_plot(parser_dict, amplicon_fasta, coords, ref_num, sgRNA, log)
                topSeqs = os.path.join(mut_dir, 'topSeqs_coverage' + str(parser_dict['mutation_calling']['summary_threshold'])+ '_percent.xls')
                pre_alignment_plot(output_dir, topSeqs, 'topSeqs', parser_dict, seq_range, default, zero_site, ref_lines, target, log)

        if run_qc:
            from summary import main as summary_main
            log.readout('Creating summary pdf')
            summary_main(parser_dict, workflow, target, qc_dir, genome_bam_dir, paired_end, fastq_dir, R1_suffix, R2_suffix, samtools, log)
        if clean_intermediates:
            from summary import cleanup_intermediate_folders
            log.readout('Removing intermediate files and folders')
            cleanup_intermediate_folders(output_dir,parser_dict,target)
        log.readout('Completed summary of %s', target)
        log.removeFilter(summary)
        return

    args = parseArgs()

    pipeline_config = args['pipeline_config'] 
    parser=Parser()
    if int(sys.version[0])==3:
        try:
            parser.read_file(open(pipeline_config)) #reading in the config file
        except IOError as e:
            print('Missing config file!')
            sys.exit()
        try:
            parser.read_file(open(parser.get('i/o', 'static_config')))
        except configparser.NoOptionError:
            pass
    elif int(sys.version[0])==2:
        try:
            parser.readfp(open(pipeline_config))
        except IOError as e:
            print('Missing config file!')
            sys.exit()
        try:
            parser.readfp(open(parser.get('i/o', 'static_config')))
        except ConfigParser.NoOptionError:
            pass

    parser_dict = {s:dict(parser.items(s)) for s in parser.sections()}
    default_dict = default_config()
    parser_dict = validate_config(parser_dict, default_dict)
    if not args['skip_program_check']:
        validate_programs(parser_dict)
    #########################
    ##paths and parameters###
    #########################

    ###various parameter setup to save me from referencing parser_dict all the time 
    workflow = parser_dict['modules']['workflow'].lower()
    runmode = parser_dict['runtime']['runmode'].lower()
    run_qc = parser_dict['modules']['run_qc'] 
    alignment = parser_dict['modules']['alignment']
    paired_end = parser_dict['processing']['paired_end']
    clean_intermediates = parser_dict['processing']['clean_intermediates']
    output_dir = parser_dict['i/o']['output_directory']
    samtools = parser_dict['programs']['samtools']
    fastq_dir = parser_dict['i/o']['fastq_directory'] #currently loops through fastqs found in specified folder and matches suffix
    R1_suffix = parser_dict['processing']['r1_suffix']
    R2_suffix = parser_dict['processing']['r2_suffix']
    merge_mut_tables = parser_dict['processing']['merge_mut_tables']
    n_processes=parser_dict['multiprocess']['n_processes']

    if workflow not in ['standard', 'gapped', 'offtarget']:
        print(parser_dict['modules']['workflow'] + ' is not a valid workflow. Please select from standard, gapped and offtarget.')
        sys.exit()
    if runmode not in ['project', 'single_sample', 'summary']:
        print(parser_dict['runtime']['runmode'] + ' is not a recognized runmode. Please select from project, single_sample and summary.')
        sys.exit()
    if runmode =='single_sample' and args['fastq1']=='' and args['fastq2']=='':
        print('Runmode set to single sample but no R1/R2 fastqs given')
        sys.exit()

    log_folder = os.path.join(output_dir, 'logs')
    try:
        os.makedirs(log_folder)
    except OSError as e:
        if e.errno!=errno.EEXIST:
            raise

    if workflow == 'gapped' and not paired_end:
        print('Please use the standard workflow for single end reads.')
        sys.exit()

    if runmode == 'project':
        log_file = os.path.join(log_folder, 'general.' + workflow + '_run.log')
    elif runmode == 'single_sample':
        log_file = os.path.join(log_folder, 'sp.'+workflow + '_run.' + re.sub(R1_suffix, '', os.path.basename(args['fastq1'])))
    elif runmode == 'summary':
        log_file = os.path.join(log_folder, workflow + '_run.summary.log')

    print('Beginning logging...')

    class LogFilter(logging.Filter):
        def __init__(self):
            self.sample_name=''
        def filter(self, record):
            record.app_name = self.sample_name
            return True
    main=LogFilter()
    main.sample_name='Main'
    fmt = '%(asctime)s %(levelname)-8s %(app_name)s %(message)s' 

    logging.addLevelName(35, 'READOUT')
    def readout(self, message, *args, **kws):
        if self.isEnabledFor(35):
            self._log(35, message, args, **kws) 
    logging.Logger.readout = readout

    rootLogger = logging.getLogger()
    rootLogger.setLevel('DEBUG')
    rootLogger.addFilter(main)
    logFormatter = logging.Formatter(fmt=fmt, datefmt='%Y/%m/%d %I:%M:%S %p')
    
    fileHandler = logging.FileHandler(log_file, mode='w')
    fileHandler.setFormatter(logFormatter)
    rootLogger.addHandler(fileHandler)
    
    consoleHandler = logging.StreamHandler(sys.stdout)
    consoleHandler.setFormatter(logFormatter)
    consoleHandler.setLevel(logging.WARNING)
    rootLogger.addHandler(consoleHandler)

    mut_suffix = '.mut.xls'

    if 'fastq_directory' not in parser_dict['i/o']:
        if runmode=='single_sample':
            parser_dict['i/o']['fastq_directory']=''

    if runmode !='summary':
        file_R1, file_R2 = fastq_setup(fastq_dir, R1_suffix, R2_suffix, paired_end, args['fastq1'], args['fastq2'], runmode, rootLogger)
        file_R1R2=list(zip(file_R1, file_R2)) #'base' level of this object will be just with files R1 and R2 zipped. will permutate with target/primers/fasta as necessary
 
    try: 
        os.makedirs(output_dir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
    if run_qc or workflow=='offtarget':
        genome_bam_dir = os.path.join(output_dir, 'genome_bams')
        if runmode!='summary':
            logging.info('Genome bams will be located in %s', genome_bam_dir)
    else:
        genome_bam_dir = ''
    if run_qc and not runmode=='summary' or workflow=='offtarget' and not runmode=='summary':
        bwa = parser_dict['programs']['bwa']
        genome_fasta = parser_dict['reference_files']['genome_fasta']
        if genome_bam_dir!='':
            try:
                os.makedirs(genome_bam_dir)
            except OSError as e:
                if e.errno != errno.EEXIST:
                    raise

        index_files = ['.pac', '.amb', '.ann', '.bwt', '.sa']
        index=False
        for ext in index_files:
            if not os.path.isfile(genome_fasta+ext):
                index=True
        if index:
            link_fasta = os.path.join(output_dir, os.path.basename(genome_fasta))
            try:
                os.link(genome_fasta, link_fasta)
            except OSError as e:
                if e.errno==errno.EXDEV:
                    try:
                        s='Genome fasta is on an external device, copying to ' + link_fasta + ' instead of forming a symbolic link. This will take longer.'
                        logging.warning('%s',s)
                        shutil.copy(genome_fasta, link_fasta)
                    except shutil.Error:
                        raise
                if e.errno != errno.EEXIST:
                    logging.critical(genome_fasta + ' does not exist!')
                    raise
            index_files = ['.pac', '.amb', '.ann', '.bwt', '.sa']
            index=False
            for ext in index_files:
                if not os.path.isfile(link_fasta + ext):
                    index=True
            if index and runmode=='single_sample':
                s='Unable to avoid running bwa index multiple times if you parallelize across samples. To avoid conflicts and excessive runtimes, please copy/link the index files to ' + out_dir +  ' or run the following command before rerunning this script: ' + bwa + ' index ' +  link_fasta
                logging.error('%s', s)
                sys.exit()
            if index:
                with open(os.devnull, 'w') as f:
                    cmd = bwa + ' index ' + link_fasta
                    subprocess.check_call(cmd.split(), stdout=f, stderr=f)
            genome_fasta = link_fasta
    
    if workflow=='offtarget':
        file_R1R2 = [[i] for i in file_R1R2]

    if not workflow=='offtarget':
        ref_file = parser_dict['reference_files']['amplicon_reference']
        try:
            target, primer1, primer2, amplicon, sgRNA, n_ref = get_reference(ref_file)
            #if not paired_end:
            #    primer2=['' for i in primer2] #probably can't guarantee that that primer2 should always be empty, i.e. if they were using pre-merged fastqs
        except IOError:
            logging.error(ref_file + ' not found!')
            sys.exit()
            
        per_ref_details = []
        for j in range(n_ref): #for each amplicon/line in the reference file
            try: #should subvert issues with running wrapper in parallel and trying to create the same directory
                os.makedirs(os.path.join(output_dir, target[j]))
            except OSError as e:
                if e.errno != errno.EEXIST:
                    raise
        
            mut_dir=os.path.join(output_dir, target[j], 'mut')
            bam_dir = os.path.join(output_dir, target[j], 'bam')
            qc_dir = os.path.join(output_dir, target[j], 'qc')
            processed_fasta_dir = os.path.join(output_dir, target[j], 'assembled_reads')
            if workflow=='standard' and paired_end:
                try:
                    os.makedirs(processed_fasta_dir)
                except OSError as e:
                    if e.errno != errno.EEXIST:
                        raise
            amplicon_fasta, rc = create_amplicon_fasta(ref_file, paired_end, target[j], primer1[j], primer2[j], amplicon[j], os.path.join(output_dir, target[j]), runmode)
            ##index amplicon fasta
            if not os.path.exists(amplicon_fasta+'.fai'):
                rootLogger.readout('Indexing amplicon_fasta %s with samtools faidx', amplicon_fasta)
                cmd = samtools + ' faidx ' + amplicon_fasta
                subprocess.check_call(cmd.split())
        
            if workflow=='gapped':
                rootLogger.readout('Indexing amplicon_fasta %s with bwa index', amplicon_fasta)
                bwa = parser_dict['programs']['bwa']
                index_files = ['.pac', '.amb', '.ann', '.bwt', '.sa']
                index=False
                for ext in index_files:
                    if not os.path.isfile(amplicon_fasta + ext):
                        index=True
                if index:
                    rootLogger.readout('Missing necessary index files for bwa, indexing amplicon_fasta %s with bwa index', amplicon_fasta)
                    with open(os.devnull, 'w') as f:
                        cmd = bwa + ' index ' + amplicon_fasta 
                        subprocess.check_call(cmd.split(), stdout=f, stderr=f)
                else:
                    rootLogger.readout('Necessary index files for bwa found, no need to index')
            
            for path in [mut_dir, bam_dir, qc_dir]:
                try:
                    os.makedirs(path)
                except OSError as e: 
                    if e.errno != errno.EEXIST:
                        raise
            per_ref_details.append((j, target[j], mut_dir, bam_dir, qc_dir, processed_fasta_dir, amplicon_fasta,primer1[j], primer2[j], sgRNA[j]))
        if not runmode=='summary':    
            file_R1R2=list(itertools.product(file_R1R2, per_ref_details))
            file_R1R2=sorted(file_R1R2, key=lambda x: x[1][0])
    if runmode=='project':
        if mp.cpu_count()==1 and n_processes>1:
            logging.warning('Only one core detected, multiprocessing not available') 
            n_processes=1   
        if n_processes==1:
            file_R1R2=list(itertools.product(file_R1R2, [False]))
            rootLogger.readout('Looping over sample level steps...')
            list(map(sample_level_processing, file_R1R2))
            if not workflow=='offtarget':
                list(map(summary_steps, per_ref_details))
        elif n_processes>1:
            file_R1R2=list(itertools.product(file_R1R2, [True]))
            original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
            timeout = 60000
            pool = mp.Pool(min(n_processes, len(file_R1R2), mp.cpu_count()))
            signal.signal(signal.SIGINT, original_sigint_handler)
            try:
                rootLogger.readout('Multiprocessing sample level steps...')
                pool.map_async(sample_level_processing, file_R1R2).get(timeout)
            except KeyboardInterrupt:
                pool.terminate()
                sys.exit()
            else:
                pool.close()
            pool.join()
            try:
                merge_mp_logs(log_folder,workflow)
            except IOError:
                rootLogger.warning('Unable to merge multiprocessed sample log files!')
            rootLogger.removeHandler(fileHandler)
            fileHandler.flush()
            fileHandler.close()
            fileHandler = logging.FileHandler(log_file, mode='a')
            fileHandler.setFormatter(logFormatter)
            rootLogger.addHandler(fileHandler)
            rootLogger.readout('Successfully completed multiprocessed steps!')

            if not workflow=='offtarget':
                list(map(summary_steps, per_ref_details))
    if runmode=='single_sample':
        file_R1R2=list(itertools.product(file_R1R2,[False]))
        list(map(sample_level_processing,file_R1R2))
    if runmode=='summary':
        try:
            merge_sp_logs(log_folder,workflow)
        except IOError:
            rootLogger.warning('Unable to merge single_sample log files!')
        list(map(summary_steps, per_ref_details))
        try:
            merge_summary_logs(log_folder,workflow)
        except IOError:
            rootLogger.warning('Unable to merge general and summary logs!')
