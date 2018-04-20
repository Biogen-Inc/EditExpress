import sys
import re
import os
import fnmatch
import subprocess
from Bio import SeqIO
from copy import deepcopy
import pysam
from editexpress import natural_sortkey, get_reference, reverse_complement
import errno
import six

def getIndels(bam, q_cutoff):

    pos_I = {}
    pos_D = {}
    qname_ID = {}
    ID_qname = {}
    # Record inserted seq if it is an insertion
    insert_seq = {}
    # For checking whether the indel is in the begining or end of reads
    read_pos_I = {}
    read_pos_D = {}
    n_reads_total = 0;
    n_reads_used = 0
    n_reads_ASfailed = 0
    reads = pysam.AlignmentFile(bam, 'rb')
    for read in reads:
        n_reads_total += 1
        if read.tid > -1:
            chr = read.reference_name
            start = read.reference_start
            mapq = read.mapping_quality
            cigar = read.cigarstring
            seq = read.query_sequence
            AS = read.get_tag('AS')
            qname = read.query_name
            if mapq >= q_cutoff and re.search('[ID]', cigar) != None:
                digits = [int(i) for i in re.findall('\d+',cigar)]
                char = re.findall('\D+', cigar)
                cnt_D = 0
                digits_t = deepcopy(digits)
                digits_tt = deepcopy(digits)
                # The alignment start position does not include the soft clip (S) and hard clip (H). 
                #'I' should be not be counted for genomic position
                for i in range(len(char)):
                    if char[i] in ['S', 'I', 'H']:
                        digits_t[i] = 0
                # The sequence should not include the deleted sequence.
                for i in range(len(char)):
                    if char[i] in ['D', 'H']:
                        digits_tt[i] = 0
                for i in range(len(char)):
                    # Insertion or deletion
                    if char[i] == 'I' or char[i] == 'D':
                        if i == 0:
                            start_t = start
                        else:
                            start_t = start + sum(digits_t[:i])
                           
                        end_t = start_t + digits[i] - 1

                        key = chr + '-' + str(start_t) + '-' + str(end_t)
                        key_size = end_t - start_t + 1
        
                        read_start = sum(digits_tt[:i]) - cnt_D
                        read_end = sum(digits_tt[:i+1]) - cnt_D

                        if char[i] == 'I':
                            if key in pos_I:
                                pos_I[key] += 1
                                insert_seq[key].append(seq[read_start:read_end])
                                read_pos_I[key].append(read_start)
                                ID_qname[key].append(qname)
                            else:
                                pos_I[key] = 1
                                insert_seq[key] = [seq[read_start:read_end]]
                                read_pos_I[key] = [read_start]
                                ID_qname[key] = [qname]
                            # Keep track of read name
                            if key in list(six.itervalues(qname_ID)):
                                qname_ID[qname].extend(key)
                            else:
                                qname_ID[qname] = [key]
                        else:
                            cnt_D = cnt_D + digits[i]
                            if key in pos_D:
                                pos_D[key] += 1
                                read_pos_D[key].append(read_start)
                                ID_qname[key].append(qname)
                            else:
                                pos_D[key] = 1
                                read_pos_D[key] = [read_start]
                                ID_qname[key] = [qname]
                            # Keep track of read name
                            if key in list(six.itervalues(qname_ID)):
                                qname_ID[qname].extend(key)
                            else:
                                qname_ID[qname] = [key]

    reads.close()
    # Get the most common insert base if they are different
    for k,v in insert_seq.items():
        insert_seq[k] = max(set(v), key=v.count)

    # Get the mean read position of the indel
    for k,v in read_pos_I.items():
        read_pos_I[k] = sum(v)/len(v)

    for k,v in read_pos_D.items():
        read_pos_D[k] = sum(v)/len(v)

    return(pos_I, pos_D, read_pos_I, read_pos_D, insert_seq, qname_ID, ID_qname)


def indelCaller(bam, q_cutoff, percent_cutoff, n_support, out_dir, samtools, bedtools):
    # Create tmp folder for intermediate files generated
    tmp_dir = os.path.join(out_dir, 'tmp')
    try:
        os.makedirs(tmp_dir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    [I_dict, D_dict, read_pos_I, read_pos_D, insert_seq, qname_ID, ID_qname] = getIndels(bam, q_cutoff)

    # Write indel locations into bed file, include both start and end
    # Note that the first base is numbered 0 in bed but 1 in sam, chromend is not included in bed
    sample_name = os.path.basename(bam).split('.')[0]
    indel_bed = os.path.join(tmp_dir, sample_name + '.indels.bed')
    fout = open(indel_bed, 'w')
    for k in I_dict:
        splits = k.split('-')
        splits[1] = str(splits[1])
        splits[2] = str(splits[2])
        l = "\t".join(splits) + '\t' + str(I_dict[k]) + '\t' + 'Insert'
        six.print_(l, file=fout)

    for k in D_dict:
        splits = k.split('-')
        splits[1] = str(splits[1])
        splits[2] = str(splits[2])
        l = "\t".join(splits) + '\t' + str(D_dict[k]) + '\t' + 'Del'
        six.print_(l, file=fout)
    fout.close()

    # Calculate coverage for indel region
    cov = os.path.join(tmp_dir, sample_name+'.cov.xls')
    temp_bam = re.sub('\.bam', '.temp.bam', bam)
    cmd = samtools + ' view -hbu -q ' + str(q_cutoff) + ' -o ' + temp_bam + ' ' +  bam
    subprocess.check_call(cmd.split())

    with open(cov, 'w') as f:
        cmd = bedtools + ' coverage -counts -a ' + indel_bed + ' -b ' + temp_bam
        subprocess.check_call(cmd.split(), stdout=f)

    os.remove(temp_bam)
    # Call indels
    f_called = open(os.path.join(tmp_dir, sample_name + '.' + 'indels_called.xls'), 'w')
    f_notcalled = open(os.path.join(tmp_dir, sample_name + '.' + 'indels_notcalled.xls'), 'w')
    
    fcov = open(cov, 'r+')
    for line in fcov:
        splits = line.rstrip().split('\t')
        chr = splits.pop(0)
        indel_type = splits.pop(3)
        splits = [int(x) for x in splits]
        indel_size = splits[1] - splits[0] + 1
        indel_key = chr + '-' + str(splits[0]) + '-' + str(splits[1])
       
        # Inserted sequence
        insert_seq_str = 'None'
        if indel_type == 'Insert':
            insert_seq_str = insert_seq[indel_key]
 
        # Indel position
        indel_start = splits[0] - 1 # requries a base before the indel
        indel_end = splits[1]
            
        # Number of reads spanning the indel and total coverage for this target
        count_ref = splits[3]-splits[2]
        count_indel = splits[2]
        count_total = splits[3]

        error_catch_file = os.path.join(tmp_dir, sample_name + '.' + 'error.xls')
        error_f = open(error_catch_file, 'w') 
 
        if count_indel >= n_support:
            if count_indel <= count_total:
                percent = float(count_indel)*100/count_total
                if count_ref < 0:
                    six.print_("error: negative reference coverage! " + str(count_ref), file=error_f)

                if percent >= percent_cutoff:
                    qnames = ID_qname[indel_key]
                    ids = []
                    for q in qnames:
                        ids.extend(qname_ID[q])
                    ids = set(ids)
                    six.print_(chr + '\t' + str(indel_start) + '\t' + str(indel_size) + '\t' + indel_key + '\t' + str(count_indel) +'\t' + str(count_total) + '\t' + str(round(percent,4)) + '\t' + indel_type + '\t' + insert_seq_str + '\t' + ','.join(ids), file=f_called)
            else:
                splits = [str(x) for x in splits]
                six.print_(chr + '\t' + "\t".join(splits) + '\t' + indel_type + '\t' + insert_seq_str, file=f_notcalled)
    
    f_called.close()
    f_notcalled.close()

    # Sort output .xls file based on %_Coverage
    result_file = os.path.join(out_dir, sample_name + '.indel.xls')
    f_result=open(result_file, 'w')
    f_result.write('\t'.join(['Chr', 'Pos', 'Indel_Size', 'ID', 'N_Alt', 'N_Total', 'Alt_cov_%', 'Type', 'Inserted_Seq', 'Indels_Connected'])+'\n')
    f_result.close()
    with open(result_file, 'a+') as f:
        cmd = 'sort -nr -k 7 ' + os.path.join(tmp_dir, sample_name + '.' + 'indels_called.xls')
        subprocess.check_call(cmd.split(), stdout=f)

