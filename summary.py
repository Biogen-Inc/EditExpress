#!/usr/bin/python
import sys
import re
import os
import subprocess
import numpy
import matplotlib
matplotlib.use('Agg')
import fnmatch 
import matplotlib.pyplot as pyplot
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
import reportlab
from reportlab.lib import colors
from reportlab.lib.enums import TA_JUSTIFY, TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import letter
from reportlab.platypus import *
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.pdfbase import pdfmetrics
from editexpress import natural_sortkey, get_reference, reverse_complement
import pandas
import math
import shutil
import string
import gzip
import six

def main(parser_dict, workflow, target, directory, genome_bam_dir, paired_end, fastq_directory, R1_suffix, R2_suffix, samtools, log):
    """contains all QC summary functions, namely the summary report generated 
    using reportlab"""    
    def line_count(inputfile, zipped=False):
        """counts lines"""
        lines=0
        if not zipped:
            for line in open(inputfile):
                lines+=1
        elif zipped:
            for line in gzip.open(inputfile):
                lines+=1
        return lines


    percent_primered=parser_dict['qc_parameters']['percent_primered']
    percent_mapped_amplicon = parser_dict['qc_parameters']['percent_mapped_amplicon']
    percent_mapped_genome=parser_dict['qc_parameters']['percent_mapped_genome']
    percent_max_rl=parser_dict['qc_parameters']['percent_max_rl']
    percent_assembled=parser_dict['qc_parameters']['percent_assembled']
    paired_end = parser_dict['processing']['paired_end']
    
    count_suffix = '.genome.count'

    pyplot.close('all') #should ensure no weird instances where other plots were not closed properly

    counts=os.path.join(directory,'counts')
    read_lengths=os.path.join(directory, 'read_lengths')
    bam_dir = re.sub('/qc', '/bam', directory)

    assembled_reads = re.sub('/qc', '/assembled_reads', directory)
    primered_dir = re.sub('/qc', '/primered_reads', directory)

    sam = os.path.join(directory, 'sam')
    samples = []

    log.info('Combining counts tables...')
    counts_df = None
    try:
        for files in os.listdir(counts):
            if fnmatch.fnmatch(files, '*'+count_suffix):
                samples.append(os.path.join(counts, files))
                temp=pandas.read_csv(os.path.join(counts, files), header=0, sep='\t', index_col=None)
                temp.columns=('Geneid', re.sub(count_suffix, '', os.path.basename(files)))
                if counts_df is None:
                    counts_df = temp
                else:
                    counts_df = pandas.merge(counts_df, temp, on="Geneid")
    except OSError:
        log.critical('Error! Missing counts directory. Please regenerate to perform the QC summary')
        sys.exit()
    rowsum=pandas.DataFrame(counts_df.sum(1))
    counts_df2=counts_df.merge(rowsum, left_index=True, right_index=True)
    counts_df2.sort_values(by=0, ascending=False)
    del(counts_df2[0])
    counts_df2.to_csv(os.path.join(directory, 'merged.counts'), sep='\t', index=False)
    
    
    samples=sorted(samples, key=natural_sortkey)

    #####
    #setting up reportlab - general format is to create a "Story" or list and append paragraphs, tables, figures, spacings, etc
    #####
    log.info('Beginning Reportlab document...')
    doc = SimpleDocTemplate(os.path.join(directory, "summary.pdf"),pagesize=letter,rightMargin=72,leftMargin=72,topMargin=72,bottomMargin=18)
    Story=[]
    styles=getSampleStyleSheet()
    ##Page 1 - always written. Header, list of samples, parameter descriptions. If there are tons of samples then 'page1' can be multiple pages long
    p = Paragraph("Input Parameters", styles['Title'])
    Story.append(p)
    Story.append(Spacer(1,12))
    samples2=[]
    for i in range(len(samples)):
        samples2.append(re.sub(count_suffix, '', os.path.basename(samples[i])))
    samples2=(', ').join(samples2) #comma separated list to output on the main page (just a list of all samples). used to be <br/> separated but got unwieldy when there were many samples
    
    p = Paragraph(str('Sample list:'), styles['Heading4'])
    Story.append(p)
    p = Paragraph(samples2, styles['Normal'])
    Story.append(p)
    p=Paragraph(("<br/>All sample level files are located in " + os.path.abspath(directory) + " including sam files, counts, assembled fastas and read lengths"), styles['Heading4'])
    Story.append(p)
    p=Paragraph(("<br/>Parameters<br/>"), styles['Heading4'])
    Story.append(p)
    p=Paragraph((("Greater than ") + str(percent_mapped_amplicon)+ "% of processed reads mapped to the amplicon.<br/>At least " + str(percent_max_rl) + "% of R1/R2 reads were within 90% of the maximum read length<br/>"), styles['Bullet'])
    Story.append(p)
    if os.path.exists(assembled_reads) and workflow=='standard' and paired_end:
        p=Paragraph(("After merging R1/R2, " + str(percent_primered) + "% of resulting reads were correctly primered<br/>"), styles['Normal'])
        Story.append(p)
    elif os.path.exists(primered_dir):
        p=Paragraph((str(percent_primered) + "% of reads were correctly primered<br/>"), styles['Normal'])
        Story.append(p)

    p=Paragraph(("After mapping to the reference genome at least " + str(percent_mapped_genome) + "% of reads were mapped<br/>"), styles['Normal'])
    Story.append(p)
    p=Paragraph("<br/>Samples which failed these parameters, if any, appear on subsequent pages", styles['Heading3'])
    Story.append(p)

    Story.append(PageBreak())
    ########
    
    if R1_suffix.endswith('.gz'):
        zipped=True
    else:
        zipped=False
    #additional pages: only written if a parameter flag is set to true
    for sample in samples:
        sample_name = re.sub(count_suffix, '', os.path.basename(sample))
        log.info('Trying ' + sample_name + '...')
        flag=False
        genome_flag=False
        amplicon_flag = False
        primer_flag=False
        assembly_flag=False
        rl_flag=False
        primer1_flag=False
        primer2_flag=False
        #use samtools to count mapped reads and total reads (-F4 flag is filter out unmapped)
        bam = os.path.join(bam_dir, sample_name + '.bam')
        if not os.path.getsize(bam)==0:
 
            if workflow =='gapped':
                cmd = samtools + ' view -c -F260 -f64 ' + bam #filter out secondary, unmapped. filter to first in pair
            else:
                cmd = samtools + ' view -c -F260 ' + bam #-f64 flag doesn't work with output from needle since it's unpaired
            amplicon_mapping = int(subprocess.check_output(cmd.split()))
        else:
            amplicon_mapping = 0

        r1_reads = int(line_count(os.path.join(fastq_directory,sample_name+R1_suffix),zipped))//4
        total_reads=r1_reads
        if float(100.0*amplicon_mapping/total_reads) < percent_mapped_amplicon:
            amplicon_flag = True
            flag = True
            log.info(sample_name + ' has less than ' + str(percent_mapped_amplicon) + '% reads mapped to amplicon')
        if os.path.exists(assembled_reads) and workflow=='standard' and paired_end: #assembly rate and primer rate for merged reads

            #various read counts
            r2_reads = int(line_count(os.path.join(fastq_directory,sample_name+R2_suffix),zipped))//4

            assembled=os.path.join(assembled_reads, sample_name+'.merged.assembled.fastq')
            n_assembled = line_count(assembled)/4
            assembly_rate = float(n_assembled)/r1_reads * 100.0
            try: 
                primered = int(line_count(os.path.join(assembled_reads, sample_name + '.primered.fasta')))//2
                primer_out = 100.0*primered/n_assembled
            except ZeroDivisionError:
                primer_out = 0.0

            if primer_out < percent_primered:
                flag = True
                primer_flag=True
                log.info(sample_name + ' has less than ' + str(percent_primered) + '% reads containing primer sequences')
            if assembly_rate < percent_assembled:
                flag = True
                assembly_flag=True
                log.info(sample_name + ' has less than ' + str(percent_assembled) + '% reads assembled')

            r1_filtered = os.path.join(assembled_reads, sample_name+'_R1.filtered.primered.fasta')
            r2_filtered = os.path.join(assembled_reads, sample_name+'_R2.filtered.primered.fasta')
            
            n_r1_primered = int(line_count(r1_filtered))//2
            n_r2_primered = int(line_count(r2_filtered))//2

            r1_primer_rate = float(n_r1_primered)/r1_reads*100.0
            r2_primer_rate = float(n_r2_primered)/r2_reads*100.0
            if r1_primer_rate < percent_primered:
                flag=True
                primer1_flag=True
                log.info(sample_name + ' has less than ' + str(percent_primered) + "% R1 reads containing 5' primer sequence")
            if r2_primer_rate < percent_primered:
                flag=True
                primer2_flag=True
                log.info(sample_name + ' has less than ' + str(percent_primered) + "% R2 reads containing 3' primer sequence")


        if os.path.exists(primered_dir): #primer rate for unmerged reads
            primered = int(line_count(os.path.join(primered_dir, sample_name + '.primered.fasta')))//2
            primer_out = 100.0*primered/total_reads
            if primer_out < percent_primered:
                flag=True
                primer_flag=True
                log.info(sample_name + ' has less than ' + str(percent_primered) + '% reads containing primer sequence')


        sam = os.path.join(genome_bam_dir, sample_name + '.genome.sam') 
        if paired_end:
            try:
                cmd = samtools + ' view -Sc -F256 -f64 ' +  sam #filter out secondary alignments, filter to first in pair. should give a count of alignments with unique read names. probably not necessary as the count should be equal to the number of R1 or R2 reads, but just in case it isn't
                total=int(subprocess.check_output(cmd.split()))        
                cmd = samtools + ' view -Sc -F260 -f64 ' +  sam  #same as above but only mapped reads
                matches = int(subprocess.check_output(cmd.split()))
                percent_match = 100.0*matches/total
                if percent_match<percent_mapped_genome:
                    flag=True
                    genome_flag=True
                    log.info(sample_name + ' has less than ' + str(percent_mapped_genome) + '% reads mapping to the given genome')
            except ValueError:
                total='NA'
                matches = 'NA'
                flag=True
                genome_flag=True
                log.info(sample_name + ' produced a SAM file with errors when aligning to given genome')
        else: #first in pair flag is not set for single end reads...
            try:
                cmd = samtools + ' view -Sc -F256 ' +  sam #filter out secondary alignments, filter to first in pair. should give a count of alignments with unique read names. probably not necessary as the count should be equal to the number of R1 or R2 reads, but just in case it isn't
                total=int(subprocess.check_output(cmd.split()))
                cmd = samtools + ' view -Sc -F260 ' +  sam  #same as above but only mapped reads
                matches = int(subprocess.check_output(cmd.split()))
                percent_match = 100.0*matches/total
                if percent_match<percent_mapped_genome:
                    flag=True
                    genome_flag=True
                    log.info(sample_name + ' has less than ' + str(percent_mapped_genome) + '% reads mapping to the given genome')
            except ValueError:
                total='NA'
                matches = 'NA'
                flag=True
                genome_flag=True
                log.info(sample_name + ' produced a SAM file with errors when aligning to given genome')

        #read in read lengths file(s)
        R1 = numpy.genfromtxt(os.path.join(read_lengths, sample_name + '_R1.lengths'), dtype='int')
        R1_percentile = 100.0*((0.01*90*max(R1)<R1) & (R1 <= max(R1))).sum() / len(R1) #checks to see if most reads are within 90% of maximum read length. did not use max in case there were a few instances of erroneously long reads, and did not want to ignore reads that are short only a few bp

        if paired_end:
            R2 = numpy.genfromtxt(os.path.join(read_lengths, sample_name + '_R2.lengths'), dtype='int')
            R2_percentile = 100.0*((0.01*90*max(R2)<R2) & (R2 <= max(R2))).sum() / len(R2)
        if R1_percentile < percent_max_rl or paired_end and R2_percentile < percent_max_rl:
            flag=True
            rl_flag=True
            log.info(sample_name + ' has less than  ' + str(percent_max_rl) + '% R1 or R2 reads close to the maximum read length')

        #generate read length histogram(s)
        bins=(max(R1)-min(R1))//2
        if bins==0:
            bins=10
        pyplot.hist(R1, bins, facecolor='green')
        pyplot.title(sample_name + '_R1 length distribution')
        pyplot.xlabel('Read length', fontsize=16)
        pyplot.ylabel('Read count', fontsize=16)
        pyplot.gcf()
        jpg_error=False
        try:
            pyplot.savefig(os.path.join(read_lengths, sample_name + '_R1.jpg'))
        except ValueError:
            jpg_error=True
            pyplot.savefig(os.path.join(read_lengths, sample_name + '_R1.png'))

        pyplot.close()
        if paired_end:
            bins=(max(R2)-min(R2))//2
            if bins==0:
                bins=10
            pyplot.hist(R2, bins, facecolor='red')
            pyplot.title(sample_name + '_R2 length distribution')
            pyplot.xlabel('Read length', fontsize=16)
            pyplot.ylabel('Read count', fontsize=16)
            pyplot.gcf()
            try:
                pyplot.savefig(os.path.join(read_lengths, sample_name + '_R2.jpg'))
            except ValueError:
                pyplot.savefig(os.path.join(read_lengths, sample_name + '_R2.png'))
            pyplot.close()
        if jpg_error:
            log.warning('Unable to save jpg, Pillow is likely not installed. Saving as png instead. Pngs cannot be added to the summary pdf')

        if flag: #did anything fail a parameter?
            
            R1_image=os.path.join(read_lengths, sample_name + '_R1.jpg')
            R1_len=str(len(R1))
            R1_percentile = str(round(R1_percentile, 2))
            if paired_end:
                R2_image=os.path.join(read_lengths, sample_name + '_R2.jpg')
                R2_len = str(len(R2))
                R2_percentile = str(round(R2_percentile, 2))
            if not paired_end:
                R2_len = 'N/A'
                R2_percentile = 'N/A'

            p = Paragraph(sample_name, styles['Title'])
            #list of lists are converted to tables for reportlab
            data = [["Read length summary", "N. reads R1", "N. reads R2", "% R1 full read length", "% R2 full read length"],["",R1_len, R2_len, R1_percentile, R2_percentile]]
            assembly_adj = 0
            #add additional rows if they exist/depending on workflow
            if os.path.exists(assembled_reads) and workflow=='standard' and paired_end: 
                data.append(["R1/R2 primer summary", "N. R1 primered", "R1 primer rate", "N. R2 primered", "R2 primer rate"])
                data.append(["", str(n_r1_primered), str(round(r1_primer_rate,2)), str(n_r2_primered), str(round(r2_primer_rate,2))])
                data.append(["Read assembly summary", "Assembled reads", "Assembly rate", "Primered reads", "Assembled primer rate"])
                data.append(["", str(n_assembled), str(round(assembly_rate,2)), str(primered), str(round(primer_out, 2))])
                data.append(["Amplicon mapping statistics", "Target", "N. mapped", "N. total", "Mapping rate"])
                data.append(["", "Amplicon", str(amplicon_mapping), str(total_reads), str(round(float(100.0*amplicon_mapping/total_reads),2))])
                assembly_adj = 4
            elif os.path.exists(primered_dir):
                data.append(["Primer rate summary", 'Assembled reads', 'Assembly rate',"Primered reads", "Primer rate"])
                data.append(["", "N/A", "N/A", str(primered), str(round(primer_out, 2))])
                data.append(["Amplicon mapping statistics", "Target", "N. mapped", "N. total", "Mapping rate"])
                data.append(["", "Amplicon", str(amplicon_mapping), str(total_reads), str(round(float(100.0*amplicon_mapping/total_reads),2))])
                assembly_adj = 2
            else:
                data.append(["Amplicon mapping statistics", "Target", "N. mapped", "N. Total", "Mapping rate"])
                data.append(["", "Amplicon", str(amplicon_mapping), str(total_reads), str(round(float(100.0*amplicon_mapping/total_reads),2))])
            data.append(["Genome mapping statistics", "Target", "N. mapped", "N. total", "Mapping rate"])
            data.append(["", "Whole genome", str(matches), str(total), str(round(percent_match,2))])
            t=Table(data, hAlign="LEFT")
            #add grid lines, default color to black
            t.setStyle(TableStyle([('INNERGRID', (0,0), (-1,-1), 0.25, colors.black),('BOX', (0,0), (-1,-1), 0.25, colors.black),('FONTNAME', (0, 0), (-1, 0), 'Times-Bold'),('FONTNAME', (0, 2+assembly_adj), (-1, 2+assembly_adj), 'Times-Bold'),]))
            #set color for optional fields to black
            t.setStyle(TableStyle([('FONTNAME', (0,4+assembly_adj), (-1,4+assembly_adj), 'Times-Bold')]))

            if os.path.exists(assembled_reads) and workflow=='standard' and paired_end:
                t.setStyle(TableStyle([('FONTNAME', (0,2), (-1,2), 'Times-Bold')]))
                t.setStyle(TableStyle([('FONTNAME', (0,4), (-1,4), 'Times-Bold')]))

            if os.path.exists(primered_dir):
                t.setStyle(TableStyle([('FONTNAME', (0,2), (-1,2), 'Times-Bold')]))

            #change color to red if flags failed
            if rl_flag:
                t.setStyle(TableStyle([('TEXTCOLOR', (1,1), (-1,1), colors.red)]))
            if primer1_flag:
                t.setStyle(TableStyle([('TEXTCOLOR', (1,3), (2,3), colors.red)]))
            if primer2_flag:
                t.setStyle(TableStyle([('TEXTCOLOR', (3,3), (-1,3), colors.red)]))
            if amplicon_flag:
                t.setStyle(TableStyle([('TEXTCOLOR', (1,3+assembly_adj), (-1,3+assembly_adj), colors.red)]))
            if primer_flag:
                t.setStyle(TableStyle([('TEXTCOLOR', (3,1+assembly_adj), (-1,1+assembly_adj), colors.red)]))
            if assembly_flag:
                t.setStyle(TableStyle([('TEXTCOLOR', (1,5), (2,5), colors.red)]))
            if genome_flag:
                t.setStyle(TableStyle([('TEXTCOLOR', (1,5+assembly_adj), (-1,5+assembly_adj), colors.red)]))

            Story.append(p)
            Story.append(t)
            Story.append(Spacer(1,12))
            if os.path.exists(R1_image):
                im1 = Image(R1_image, 4.5*inch, 3*inch)
                Story.append(im1)
            if paired_end and os.path.exists(R2_image):
                im2 = Image(R2_image, 4.5*inch, 3*inch)
                Story.append(im2)
            Story.append(PageBreak())
        else:
            log.info(sample_name + ' did not fail any parameters!')
    doc.build(Story)


def setup_files_topSeqs(mut_dir, parser_dict, log):
    """grabs a few parameters, writes header for topSeqs file based on 
    gapped/normal runtime, loops through .mut.xls files and runs summary_topSeqs"""
    workflow = parser_dict['modules']['workflow']
    min_cov_percent = parser_dict['mutation_calling']['summary_threshold']
    n_hits = parser_dict['mutation_calling']['topseqs_min_hits']
    topSeqs = open(os.path.join(mut_dir, 'topSeqs_coverage'+str(min_cov_percent)+'_percent.xls'), 'w')

    suffix = '.mut.xls'
    
    if workflow=='gapped':
        topSeqs.write('\t'.join(['Sample','Ref','N_mapped','Exp_WT','Exp_WT_%','Mut','N_mut','Mut_%','Frameshift', '...'])+'\n')
    else:
        topSeqs.write('\t'.join(['Sample','Ref','N_mapped','WT','WT_%','Mut','N_mut','Mut_%','Frameshift','Alignment','...'])+'\n')
    
    log.info('Creating topSeqs mutation summary...')
    samples = []
    try:
        for mut_file in os.listdir(mut_dir):
            if fnmatch.fnmatch(mut_file, '*'+suffix):
                samples.append(os.path.join(mut_dir, mut_file))
    except OSError:
        log.critical('Error! ' + mut_dir + ' not found! Please regenerate in order to summarize mutation tables')
        sys.exit()
    samples=sorted(samples, key=natural_sortkey)
    if len(samples)==0:
        log.critical('Error! No mutation tables found in ' + mut_dir + '. Please regenerate in order to summarize mutation tables')
        sys.exit()
    for sample in samples:
        sample_name=re.sub(suffix, '', os.path.basename(sample))
        log.info('Looking for topSeqs in ' + sample_name)
        summary_topSeqs(sample, sample_name, topSeqs, min_cov_percent, n_hits, workflow)
    log.info('Finished topSeqs mutation summary')
    return

def summary_topSeqs(sample, sample_name, topSeqs, min_cov_percent, n_hits, workflow):
    """Finds all sequences above user specified coverage and hit thresholds and
    outputs them in the topSeqs file. If a sample does not have any qualifying
    sequences, nothing is output in the topSeqs file"""
    summary_line=[]
    summary_line.append(sample_name)
    summary_line.append('ref')
    summary_line.append('N_total')
    summary_line.append(0)
    summary_line.append(0.0)
    muts={}
    seqs={}
    with open(sample, 'r') as out_file:
        next(out_file)
        for line in out_file:
            summary=line.rstrip().split('\t')
            total_mapped=summary[4]
            summary_line[2]=total_mapped
            seqs[summary[1]]=summary[2]
            if not workflow=='gapped':
                seqs[summary[1]]+='\t'+summary[6]
                WT = 'WT'
            else:
                WT = 'exp_WT'
            muts[summary[0]+"\t"+summary[1]]=int(summary[3])    
        ordered_out = sorted(list(six.iteritems(muts)), key=lambda x: x[1], reverse=True)
        for i in range(len(ordered_out)):
            out = ordered_out[i][0].split('\t')
            if i==0:
                summary_line[1]=out[0]
            if out[1]==WT:
                summary_line[3]=ordered_out[i][1]
                summary_line[4]=100.0*ordered_out[i][1]/int(total_mapped)
            elif (100.0*ordered_out[i][1]/int(total_mapped))>min_cov_percent and ordered_out[i][1]>=n_hits:
                summary_line.append(out[1])
                summary_line.append(str(ordered_out[i][1]))
                temp = 100.0*ordered_out[i][1]/int(total_mapped)
                summary_line.append("%.2f"%(100.0*ordered_out[i][1]/int(total_mapped)))
                summary_line.append(seqs[out[1]])
        summary_line[3]=str(summary_line[3])
        summary_line[4]="%.2f" % summary_line[4]
        if len(summary_line)>5:
            topSeqs.write('\t'.join(summary_line) + '\n')
    return 

def cleanup_intermediate_folders(path, parser_dict, target):
    """deletes intermediate files/folders if specified"""
    bam_dir = os.path.join(path,target, 'bam')
    qc_dir = os.path.join(path,target, 'qc')
    mut_dir = os.path.join(path,target, 'mut')
    amplicon_fasta = target+'.fasta'
    amplicon_fasta = os.path.join(path, target,amplicon_fasta)
    os.remove(amplicon_fasta)
    os.remove(amplicon_fasta+'.fai')
    shutil.rmtree(bam_dir)
    if os.path.isdir(os.path.join(path,target, 'assembled_reads')):
        shutil.rmtree(os.path.join(path,target, 'assembled_reads'))
    if os.path.isdir(os.path.join(path,'genome_bams')):
        shutil.rmtree(os.path.join(path,'genome_bams'))
    if os.path.isdir(os.path.join(qc_dir, 'read_lengths')):
        for list_file in os.listdir(os.path.join(qc_dir, 'read_lengths')):
            if fnmatch.fnmatch(list_file, '*.lengths'):
                os.remove(os.path.join(qc_dir, 'read_lengths' + '/' + list_file))
    if os.path.isdir(os.path.join(qc_dir, 'counts')):
        shutil.rmtree(os.path.join(qc_dir, 'counts'))
    if os.path.isfile(qc_dir + '/primered_reads'):
        os.remove(qc_dir + '/primered_reads')
    return

