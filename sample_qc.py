#!/usr/bin/python
import sys
import re
import os
import subprocess
import string
from align_sample_needle import needle_align, primered_reads, reverse_complement
import errno
from contextlib import contextmanager
import io
import shutil
import gzip
 
def distance(seq1, seq2): 
    """simple distance determination"""

    if len(seq1) != len(seq2):
        #print('Seq lengths do not match')
        #sys.exit()
        return 2
    else:
        n_mismatch = 0
        for i in range(len(seq1)):
            if seq1[i] != seq2[i]: n_mismatch += 1
           
    return n_mismatch

def count_read_lengths(infile, outfile): 
    """count length of all reads in fastq"""
    rl_out=open(outfile, 'w')
    n=0
    for line in open(infile,'r'):
        n+=1
        if n%2==0 and n%4!=0:
            rl_out.write(str(len(line)-1)+ '\n')
    rl_out.close()

    return

def main(out_dir, parser_dict, fastq1, fastq2, paired_end, amplicon_fasta, genome_fasta,genome_bam_dir, gtf, sample_name, bwa_mem_params, log):
    """counts R1/R2 fastqs, runs bwa mem against genome, featurecounts, and 
    generates a simpler count table"""
    read_lengths = os.path.join(out_dir, 'read_lengths/')

    try:
        os.makedirs(read_lengths)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
    delete=[]
    if fastq1.endswith('.gz'):
        fastq1_old=fastq1
        fastq1=os.path.join(out_dir, sample_name+'_tempR1.fq')
        delete.append(fastq1)
        with open(fastq1, 'wb') as f_out, gzip.open(fastq1_old, 'rb') as f_in:
            shutil.copyfileobj(f_in, f_out)
    if paired_end :
        if fastq2.endswith('.gz'):
            fastq2_old=fastq2
            fastq2=os.path.join(out_dir, sample_name+'_tempR2.fq')
            delete.append(fastq2)
            with open(fastq2, 'wb') as f_out, gzip.open(fastq2_old, 'rb') as f_in:
                shutil.copyfileobj(f_in, f_out)

    count_read_lengths(fastq1, os.path.join(read_lengths, sample_name+'_R1.lengths'))

    if paired_end:
        count_read_lengths(fastq2, os.path.join(read_lengths, sample_name+"_R2.lengths"))

    count_dir = os.path.join(out_dir, 'counts/')

    try:
        os.makedirs(count_dir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    bwa = parser_dict['programs']['bwa']
    featureCounts = parser_dict['programs']['featurecounts']

    sam_out = os.path.join(genome_bam_dir, sample_name + '.genome.sam')
    bam_out = os.path.join(genome_bam_dir, sample_name + '.genome.bam')
    if os.path.isfile(bam_out) and not os.path.isfile(sam_out):
        cmd = samtools + ' view -H ' + bam_out + ' > ' + sam_out
        subprocess.check_call(cmd.split())
    if not os.path.isfile(sam_out) or os.path.getsize(sam_out)==0:
        if paired_end: #always runs if paired end
            fastqs = fastq1 + ' ' + fastq2
        elif not paired_end:
            fastqs = fastq1
        with open(sam_out, 'w') as f: 
            log.info('Aligning ' + sample_name + ' to ' + genome_fasta + ' using bwa mem')
            with open(os.devnull, 'w') as g:
                cmd = bwa + ' mem ' + bwa_mem_params  + ' -v 0 -a -M ' + genome_fasta + ' ' + fastqs
                subprocess.check_call(cmd.split(), stdout=f, stderr=g)
    #featurecounts
    counts_out = os.path.join(out_dir, 'counts/', sample_name + '.genome.count')
    with open(os.devnull, 'w') as f:
        cmd = featureCounts + ' -a '+ gtf + ' -g gene_name -p -C -o ' + counts_out + ' ' + sam_out
        subprocess.check_call(cmd.split(), stdout=f, stderr=f)
    #grab only geneid and count columns
    with open(counts_out, 'r') as f:
        lines_out=[]
        next(f)
        next(f)
        for lines in f:
            line=lines.rstrip().split('\t')
            lines_out.append(line[0] + '\t' + line[6]+'\n')
    with open(counts_out, 'w') as f:
        f.write('Geneid' + '\t' + sam_out + '\n')
        for line in lines_out:
            f.write(line)

    for item in delete:
        os.remove(item)

    return

def primer_rate_R1_R2(parser_dict, amplicon_fasta, bam_dir, merged_fastq_dir, sample, fastq_R1, fastq_R2, primer1, primer2, rc, log):
    """returns primer rates for unmerged reads (i.e. how many R1 reads have 
    primer1, how many R2 reads have primer 2"""    
    def fasta_length_filter(fasta, filtered_fasta):
        """finds longest sequence in fasta file"""   
        longest_seq=len(max(open(fasta, 'r'), key=len).rstrip()) #the header line had better never be longer...
        
        f_filt = open(filtered_fasta, 'w')
        with open(fasta, 'r') as f:
            for line in f:
                line=line.rstrip()
                if line.startswith('>'):
                    readname = line
                else:
                    seq = line
                    if len(seq)>=0.9*longest_seq:
                        f_filt.write(readname + '\n' + seq + '\n')
        f_filt.close()
        return

    seqtk = parser_dict['programs']['seqtk']
   
    fasta_R1 = os.path.join(merged_fastq_dir, sample+'_R1.fasta')
    fasta_R2 = os.path.join(merged_fastq_dir, sample+'_R2.fasta')
    fasta_primered_R1=os.path.join(merged_fastq_dir, sample+'_R1.primered.fasta')
    fasta_primered_R2=os.path.join(merged_fastq_dir, sample+'_R2.primered.fasta')
    fasta_filtered_R1=os.path.join(merged_fastq_dir, sample+'_R1.filtered.primered.fasta')
    fasta_filtered_R2=os.path.join(merged_fastq_dir, sample+'_R2.filtered.primered.fasta')

    try:
        cmd = seqtk + ' seq -a ' + fastq_R1
        with open(fasta_R1, 'w') as f:
            subprocess.check_call(cmd.split(), stdout=f)
        primered_reads(fasta_R1, fasta_primered_R1, primer1, '', 'no', False)
    except subprocess.CalledProcessError as e:
        with open(fasta_primered_R1, 'w'):
            log.error('Unable to convert ' + fastq_R1 + ' to fasta')
            log.error(e)
            log.warning('Proceding with empty primered fasta')
        
    try:
        cmd = seqtk + ' seq -a ' + fastq_R2
        with open(fasta_R2, 'w') as f:
            subprocess.check_call(cmd.split(), stdout=f)
        if rc=='no': #this seems backwards but the rc parameter is relative to R1. R2 reads will start with the rc of primer2
            primered_reads(fasta_R2, fasta_primered_R2, reverse_complement(primer2),'', 'no', False)
        elif rc=='yes':
            primered_reads(fasta_R2, fasta_primered_R2, primer2,'', 'no', False) 
    except subprocess.CalledProcessError as e:
        with open(fasta_primered_R2, 'w'):
            log.error('Unable to convert ' + fastq_R2 + ' to fasta')
            log.error(e)
            log.warning('Proceding with empty primered fasta')

    try:
        fasta_length_filter(fasta_primered_R1, fasta_filtered_R1)
    except ValueError:
        with open(fasta_filtered_R1, 'w') as f:
            log.warning(fasta_primered_R1 + ' is empty! So shall ' + fasta_filtered_R1 + ' be!')
    try:
        fasta_length_filter(fasta_primered_R2, fasta_filtered_R2)
    except ValueError:
        with open(fasta_filtered_R2, 'w') as f:
            log.warning(fasta_primered_R2 + ' is empty! So shall ' + fasta_filtered_R2 + ' be!')
    return

def get_ref_and_runtime(parser_dict):
    """sets up some parameters for convenience"""

    gtf = parser_dict['reference_files']['genome_gtf']
    path = parser_dict['i/o']['output_directory']

    bwa_mem_params = '-k ' + str(parser_dict['bwa_mem']['min_seed_length']) + ' -w ' + str(parser_dict['bwa_mem']['band_width']) + ' -d ' + str(parser_dict['bwa_mem']['z_dropoff']) + ' -r ' + str(parser_dict['bwa_mem']['reseeding']) + ' -c ' + str(parser_dict['bwa_mem']['mem_discard_threshold']) + ' -A ' + str(parser_dict['bwa_mem']['matching_score']) + ' -B ' + str(parser_dict['bwa_mem']['mismatch_penalty']) + ' -O' + str(parser_dict['bwa_mem']['gap_open_penalty']) + ' -E ' + str(parser_dict['bwa_mem']['gap_extension_penalty']) + ' -L ' + str(parser_dict['bwa_mem']['clipping_penalty']) + ' -U ' + str(parser_dict['bwa_mem']['unpaired_read_pair_penalty'])

    return(gtf, bwa_mem_params)

