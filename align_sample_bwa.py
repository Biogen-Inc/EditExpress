#!/usr/bin/python
import sys
import re
import os
import subprocess

def distance(seq1, seq2): 
    """simple distance metric"""
    if len(seq1) != len(seq2):
        sys.exit()
    else:
        n_mismatch = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]: n_mismatch += 1
    return n_mismatch

def pre_process(fastq1, fastq2, paired_end): 
    """used to be more complex. rn just concatenates fastqs""" 
    if not paired_end:
        fastqs=fastq1
    else:
        fastqs=fastq1+ ' ' + fastq2
    return(fastqs)

def get_params(parser_dict): 
    """just fetches all bwa_mem alignment parameters from the parser_dict object"""
    bwa_mem_param = '-k ' + str(parser_dict['bwa_mem']['min_seed_length']) + ' -w ' + str(parser_dict['bwa_mem']['band_width']) + ' -d ' + str(parser_dict['bwa_mem']['z_dropoff']) + ' -r ' + str(parser_dict['bwa_mem']['reseeding']) + ' -c ' + str(parser_dict['bwa_mem']['mem_discard_threshold']) + ' -A ' + str(parser_dict['bwa_mem']['matching_score']) + ' -B ' + str(parser_dict['bwa_mem']['mismatch_penalty']) + ' -O' + str(parser_dict['bwa_mem']['gap_open_penalty']) + ' -E ' + str(parser_dict['bwa_mem']['gap_extension_penalty']) + ' -L ' + str(parser_dict['bwa_mem']['clipping_penalty']) + ' -U ' + str(parser_dict['bwa_mem']['unpaired_read_pair_penalty'])
    return(bwa_mem_param)

def bwa_align(path, sample, bwa_mem_param, ref_fasta, fastqs, bwa, samtools, workflow, log): 
    """runs bwa mem on either a single or paired end fastqs. also makes a sam file"""
    if workflow=='gapped':
        bam_out = os.path.join(path, sample + '.bam')
        sam_out = os.path.join(path, sample + '.sam')
    elif workflow=='offtarget':
        bam_out = os.path.join(path, sample + '.genome.bam')
        sam_out = os.path.join(path, sample + '.genome.sam')
    with open(sam_out, 'w') as f:
        log.info('Aligning ' + sample + ' using bwa mem...')
        with open(os.devnull, 'w') as g:
            cmd = bwa + ' mem ' + bwa_mem_param + ' -v 0 -a -M ' + ref_fasta + ' ' + fastqs 
            subprocess.check_call(cmd.split(), stdout=f, stderr=g)
        log.info('Successfully aligned ' + sample)
    cmd = samtools + ' view -o ' + bam_out + ' -Sb ' +sam_out
    with open(os.devnull, 'w') as f:
        subprocess.check_call(cmd.split(), stdout=f, stderr=f)
    if workflow=='gapped':
        os.remove(sam_out)
    return

