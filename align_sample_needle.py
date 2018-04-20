import sys
import re
import os
import subprocess
import fnmatch
import string
import numpy
from editexpress import natural_sortkey, get_reference, reverse_complement
import errno
 
def distance(seq1, seq2): 
    """simple distance metric"""
    if len(seq1) != len(seq2):
        return 2
    else:
        n_mismatch = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]: n_mismatch += 1
    return n_mismatch

def primered_reads(fasta, fasta_primered, primer1, primer2, rc, paired_end): 
    """filter out reads not containing primer sequence(s) or with more than 1 mismatch"""

    fa_out = open(fasta_primered, 'w')
    fa_in = open(fasta, 'r')

    ID = ''
    n_primered = 0
    n_total = 0
    if rc =='yes':
        primer2 = reverse_complement(primer2)
    primer1 = primer1.upper()
    primer2 = primer2.upper()
    for line in fa_in:
        if not re.match('^>', line):
            line = line.strip('\n')
            n_mis1 = distance(primer1, line[:len(primer1)])
            if len(primer2)>0 and paired_end:
                n_mis2 = distance(primer2, line[-len(primer2):])
            else:
                n_mis2 = 0
            if n_mis1 <= 1 and n_mis2 <= 1:
                n_primered += 1
                fa_out.write(ID)
                fa_out.write(line + "\n")
        else:
            n_total += 1
            ID = line

    fa_in.close()
    fa_out.close()

    return 


def process_fastq(assembled_dir, file_R1, file_R2, assembler_param, sample_name, seqtk, read_assembler, assembler_path, primer1, primer2, rc, log): 
    """performs read merging with w/e read merger
    if possible and then converts fastq to fasta"""
    merged_fastq_prefix = os.path.join(assembled_dir, sample_name + '.merged')
    merged_fastq = merged_fastq_prefix + '.assembled.fastq'
    fasta = os.path.join(assembled_dir, sample_name + '.fasta')    
    fasta_primered = re.sub('.fasta', '.primered.fasta', fasta)

    if read_assembler == 'pandaseq':
        try:
            cmd = assembler_path + ' -f ' + file_R1 + ' -r ' + file_R2 + ' -F -w ' + merged_fastq  + ' -g ' + os.devnull + assembler_param
            subprocess.check_call(cmd.split())
            cmd = seqtk + ' seq -a ' + merged_fastq
            with open(fasta, 'w') as f:
                try:
                    subprocess.check_call(cmd.split(), stdout=f)
                except subprocess.CalledProcessError as e:
                    log.error('Unable to convert ' + merged_fastq + ' to fasta')
                    log.error(e)
                    log.warning('Proceding without this sample')
        except subprocess.CalledProcessError as e:
            with open(fasta, 'w') as f:
                log.error('Unable to merge ' + file_R1 + ' and ' + file_R2 + ' using PANDAseq')
                log.error(e)
                log.warning('Proceding without this sample')
    elif read_assembler == 'flash':
        try:
            cmd = ' '.join([assembler_path, assembler_param, '-q', '-d', assembled_dir, '-o', sample_name+'.merged.assembled', '-t 1', file_R1, file_R2])
            with open(os.devnull) as f:
                subprocess.check_call(cmd.split(), stdout=f, stderr=f)
            os.rename(os.path.join(merged_fastq_prefix+'.assembled.extendedFrags.fastq'), merged_fastq)

            cmd = seqtk + ' seq -a ' + merged_fastq
            with open(fasta, 'w') as f:
                try:
                    subprocess.check_call(cmd.split(), stdout=f)
                except subprocess.CalledProcessError as e:
                    log.error('Unable to convert ' + merged_fastq + ' to fasta')
                    log.error(e)
                    log.warning('Proceding without this sample')
        except subprocess.CalledProcessError as e:
            with open(fasta, 'w') as f:
                log.error('Unable to merge ' + file_R1 + ' and ' + file_R2 + ' using FLASH')
                log.error(e)
                log.warning('Proceding without this sample')

    primered_reads(fasta, fasta_primered, primer1, primer2, rc, True)
    
    return

def needle_align(bam_dir, sample, amplicon_fasta, fasta_primered, needle_param, needle, samtools, log):
    """Align by Needleman-Wunsch alignment algorithm"""
    sam_out = os.path.join(bam_dir, sample + '.sam')
    bam_out = os.path.join(bam_dir, sample + '.bam')
    try:
        with open(os.devnull) as f:
            log.info('Aligning ' + sample + ' with Needle...')
            cmd = needle + ' -asequence ' + amplicon_fasta +  ' -bsequence ' + fasta_primered + ' ' + needle_param + ' -aformat sam --outfile ' + sam_out
            subprocess.check_call(cmd.split(), stdout=f, stderr=f)
            log.info('Successfully aligned ' + sample)
        cmd = samtools + ' view -t ' + amplicon_fasta + '.fai -Sb -o ' + bam_out + ' ' + sam_out  #needle sam files are lacking a header file and need the ref fasta to be converted to bam
        subprocess.call(cmd.split())
        log.info('Converted' + sam_out + ' to bam format...')
        os.remove(sam_out)
        return True
    except subprocess.CalledProcessError:
        with open(bam_out, 'w') as f:
            log.error('Error aligning fasta ' + fasta_primered + ', see details in qc summary report. Mutation calling will be set to false.')
        return False

def get_param(parser_dict):
    """just sets up the parameters for needle and whatever 
    assembler was chosen"""
    seqtk = parser_dict['programs']['seqtk']
    if parser_dict['processing']['read_merger']=='pandaseq':
        read_merger = 'pandaseq'
        assembler = parser_dict['programs']['pandaseq']
        if parser_dict['pandaseq']['phred64']:
            phred = ' -6 '
        else:
            phred = ''
        if parser_dict['pandaseq']['no_uncalled_bases']:
            N = ' -N '
        else:
            N = ''
        assembler_param = phred + N + ' -A ' + parser_dict['pandaseq']['algorithm'] + ' -o ' + str(parser_dict['pandaseq']['min_overlap']) + ' -L ' + str(parser_dict['pandaseq']['max_assembly_length']) + ' -l ' + str(parser_dict['pandaseq']['min_assembly_length']) + ' -t ' + str(parser_dict['pandaseq']['threshold_score'])
    elif parser_dict['processing']['read_merger']=='flash':
        read_merger = 'flash'
        assembler = parser_dict['programs']['flash']
        if parser_dict['flash']['allow_outies']:
            outies = '-O '
        else:
            outies = ''
        assembler_param = outies + ' -m ' + str(parser_dict['flash']['min_overlap']) + ' -M ' + str(parser_dict['flash']['max_overlap']) + ' -x ' + str(parser_dict['flash']['max_mismatch_density']) 
    needle = parser_dict['programs']['needle']
    needle_param = '-gapopen ' + str(parser_dict['needle']['gapopen']) + ' -gapextend ' + str(parser_dict['needle']['gapextend'])

    return(assembler_param, needle_param, assembler, needle, seqtk, read_merger)

def pre_process_fastqs(file_R1, file_R2, output_dir, target, sample, seqtk, read_merger, assembler_path, assembler_param, primer1, primer2, rc, log):
    """some checks on what to do for single/paired end"""
    #if paired end, merge R1 and R2 then filter by primers. else just filter by primer
    log.info('Running preprocessing steps...')
    if file_R2!="Empty": #could also use parameter "paired_end"
        # Merge R1 and R2 reads
        assembled_dir = os.path.join(output_dir, target, 'assembled_reads/')
        merged_fastq_prefix = os.path.join(assembled_dir, sample + '.merged')
        merged_fastq = merged_fastq_prefix + '.assembled.fastq'
        try:
            os.makedirs(assembled_dir)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
        log.info('Completed preprocessing steps')
        process_fastq(assembled_dir, file_R1, file_R2, assembler_param, sample, seqtk, read_merger, assembler_path, primer1, primer2, rc, log)
        # Align to reference amplicon
        fasta_primered = os.path.join(assembled_dir, sample+ '.primered.fasta')
    else:
        primer_dir = os.path.join(output_dir, target, 'primered_reads/')
        fasta = os.path.join(primer_dir, sample + '.fasta')
        try:
            os.makedirs(primer_dir)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

        cmd = seqtk + ' seq -a ' + file_R1
        with open(fasta, 'w') as f:
            subprocess.check_call(cmd.split(), stdout=f)
        fasta_primered = os.path.join(primer_dir, sample + '.primered.fasta')
        log.info('Completed preprocessing steps')
        primered_reads(fasta, fasta_primered, primer1, primer2, rc, False)
    return(fasta_primered)

