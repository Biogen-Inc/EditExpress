import sys
import re
import os
import fnmatch
import subprocess
from copy import deepcopy
import pysam
import difflib
import string
#from configparser import SafeConfigParser
import numpy
from editexpress import natural_sortkey, get_reference, reverse_complement
import six

def get_seq(mut_string, ref_lines, seq_in): 
    """regenerate an alignment string post filtering using the reference 
    and old seq"""
    def accum_sum(list_in):
        """simple running sum"""
        total = 0
        for i in list_in:
            total += i
            yield total
    I_splice=[]
    temp=mut_string.split(';')
    parse_digits=[]
    parse_chars=[]
    for muts in temp:
        parse_digits.append(int(muts.split(":")[0]))
        parse_chars.append(muts.split(":")[1])
    seq_out = ref_lines
    for j in range(len(parse_digits)):
        if parse_chars[j] in ['A', 'C', 'G', 'T', 'N']:
            seq_out = seq_out[0:parse_digits[j]-1] + parse_chars[j] + seq_out[parse_digits[j]:len(seq_out)]
        elif 'D' in parse_chars[j]:
            seq_out = seq_out[0:parse_digits[j]-1] + int(parse_chars[j].split('D')[0])*'-'+seq_out[parse_digits[j]+int(parse_chars[j].split('D')[0])-1:len(seq_out)]
    if any([r for r in temp if any(z in r for z in ['I'])]):
        parse_I=[]
        I=''
        pattern = r'\[\w+\]'
        noI=re.split(pattern, seq_in)
        I=re.findall(pattern, seq_in)
        noI_num = [len(i) for i in noI]
        I_pos = list(accum_sum(noI_num))[:-1]
        I_splice = [(n,m)for n,m in zip(I_pos, I)]
        I_splice = [(h,i) for (h,i) in I_splice if h in parse_digits]
        adj = 0
        for splice in I_splice:
            seq_out = seq_out[0:adj+splice[0]] + splice[1] + seq_out[adj+splice[0]:len(seq_out)]
            adj+=len(splice[1])
    return(seq_out)

def get_frame(mut_string): 
    """calculate frameshift (0,1, or 2) based on mut string"""
    ins_temp=[r for r in mut_string.split(';') if any(z in r for z in ['I'])]
    ins_sum=0
    if len(ins_temp)>0:
        for k in ins_temp:
            ins_sum+=[int(j) for j in re.findall('\d+', k.split(':')[1])][0]
    del_temp=[r for r in mut_string.split(';') if any(z in r for z in ['D'])]
    del_sum=0
    if len(del_temp)>0:
        for k in del_temp:
            del_sum+=[int(j) for j in re.findall('\d+', k.split(':')[1])][0]
    frame_flag=ins_sum-del_sum
    frame_flag = abs(frame_flag % 3)
    return(frame_flag)

def setup_restrict(parser_dict, coords, ref_num, log):
    """sets up various parameters and performs checks to see if alignment
    restriction (narrowing the alignment output to N bases around the cleavage
    site, as determined by guide rnas or directly)"""
    if not parser_dict['mutation_calling']['restrict_alignment']:
        return False, 0
               
    seq_range=coords.split(',')[ref_num]
    seq_range = seq_range.split('-')
    if not len(seq_range)==2:
        return False, 0
    try:
        seq_range=[int(i) for i in seq_range]
    except ValueError:
        return False, 0

    try:
        if any([number<0 for number in seq_range]):
            return False, 0
        else:
            seq_range = (min(seq_range), max(seq_range))
    except KeyError:
        return False, 0

    return(True, seq_range)

def parse_restrict_alignment(seq_range, ordered_list):
    """generates aformentioned restricted alignments and new mutation strings. 
    deletions that overlap the cutoff are retained"""
    mut_dict={}
    upstream_restrict=seq_range[0]
    downstream_restrict=seq_range[1]
    for i in range(len(ordered_list)): #pass 1: id muts and count occurences
        parse = ordered_list[i][0].split('\t')
        ref_col = parse[0]
        mut_col = parse[1]
        frame_col = parse[2]
        seq_col = parse[3]
        count = ordered_list[i][1]
        mut_list = mut_col.split(';')
        if 'I' in mut_col:
            parse_I=[]
            I=''
            pattern = r'\[\w+\]'
            noI=re.split(pattern, seq_col)
            I=re.findall(pattern, seq_col)
            noI_num = [len(i) for i in noI]
            ins_seq_parse = 0
            for j in range(len(noI)):
                for k in range(len(noI[j])):
                    if ins_seq_parse==upstream_restrict:
                        upstream_ins_parse = (j,k)
                    elif ins_seq_parse==downstream_restrict-1:
                        downstream_ins_parse = (j,k)
                    ins_seq_parse+=1
            noI[upstream_ins_parse[0]]=noI[upstream_ins_parse[0]][upstream_ins_parse[1]:len(noI[upstream_ins_parse[0]])]
            noI[downstream_ins_parse[0]]=noI[downstream_ins_parse[0]][0:downstream_ins_parse[1]+1]
            for j in range(upstream_ins_parse[0], downstream_ins_parse[0]):
                noI[j] = noI[j]+I[j] 
            seq_col = ''.join(noI[upstream_ins_parse[0]:downstream_ins_parse[0]+1])
        else:
            seq_col = seq_col[upstream_restrict:downstream_restrict]
        
        if mut_col =='WT':
            parse_mut='WT'
        else:
            temp=mut_col.split(';')
            parse_digits=[]
            parse_chars=[]
            for j in range(len(temp)):
                parse_digits.append(int(temp[j].split(":")[0]))
                parse_chars.append(temp[j].split(":")[1])
            to_delete=[]
            for j in range(len(parse_digits)): 
                if parse_digits[j]<upstream_restrict:
                    if parse_chars[j] in ['A', 'C', 'G', 'T', 'N'] or 'I' in parse_chars[j]:
                        to_delete.append(j)
                    elif 'D' in parse_chars[j]:
                        if not parse_digits[j]+int(parse_chars[j].split('D')[0])>=upstream_restrict:
                            to_delete.append(j)
                if parse_digits[j]>downstream_restrict:
                    to_delete.append(j)
            for j in sorted(to_delete, reverse=True):
                del(parse_chars[j])
                del(parse_digits[j])
            if len(parse_digits)>0:
                parse_mut=';'.join([str(m)+':'+n for m,n in zip(parse_digits,parse_chars)])
            else:
                parse_mut='WT'

        frame_flag = get_frame(parse_mut)

        if frame_flag != 0:
            parse_frame = 'yes'
        else:
            parse_frame = 'no'

        if '\t'.join([ref_col, parse_mut, parse_frame, seq_col]) not in mut_dict:
            mut_dict['\t'.join([ref_col, parse_mut, parse_frame, seq_col])]=count
        else:
            mut_dict['\t'.join([ref_col, parse_mut, parse_frame, seq_col])]+=count

    mut_dict_sorted = sorted(list(six.iteritems(mut_dict)), key=lambda x: x[1], reverse=True)
    return(mut_dict_sorted) 

def filter_reads(ordered_list, read_cnt, parse_snps, threshold, ref_seq): 
    """ A 'noise' filter. parses through output and finds total occurences of 
    mutations at a particular site. if they are below a user defined threshold 
    (default of 0.25% error rate), then they will be filtered out of the read. 
    any matching hits will then be summed"""

    mut_dict = {}
    filter_muts = []
    filter_dict = {}
    filter_out = []
    for i in range(len(ordered_list)): #pass 1: id muts and count occurences
        parse = ordered_list[i][0].split('\t')
        ref_col = parse[0]
        mut_col = parse[1]
        frame_col = parse[2]
        seq_col = parse[3]
        count = ordered_list[i][1]
        mut_list = mut_col.split(';')
        for j in range(len(mut_list)):
            if mut_list[j] in mut_dict:
                mut_dict[mut_list[j]]+=count
            else:
                mut_dict[mut_list[j]]=count   

    mut_dict_sorted = sorted(list(six.iteritems(mut_dict)), key=lambda x: x[1], reverse=True)
    for key in mut_dict:
        if (mut_dict[key]/float(read_cnt) < threshold) or (not parse_snps and any(k for k in ['A', 'C', 'G', 'T', 'N'] if k in key)):
            filter_muts.append(key)
    mut_dict={}
    seq_dict={} 
    for i in range(len(ordered_list)): #pass 2: filter low threshold muts
        parse = ordered_list[i][0].split('\t')
        ref_col = parse[0]
        mut_col = parse[1]
        frame_col = parse[2]
        seq_col = parse[3]
        count = ordered_list[i][1]
        mut_list = mut_col.split(';')
        if any([r for r in mut_list if r not in filter_muts]):
            mut_col=";".join([r for r in mut_list if r not in filter_muts])
        else:
            mut_col='WT'
        
        frame_flag=get_frame(mut_col)

        if frame_flag != 0:
            frame_col = 'yes'
        else:
            frame_col = 'no'
  
        if '\t'.join([ref_col,mut_col]) not in seq_dict:
            if mut_col=='WT':
                seq_out = ref_seq
            else:
                seq_out = get_seq(mut_col, ref_seq, seq_col)
            seq_dict['\t'.join([ref_col,mut_col])]= '\t'.join([frame_col, seq_out])
        if '\t'.join([ref_col, mut_col]) not in mut_dict:
            mut_dict['\t'.join([ref_col, mut_col])]=count
        else:
            mut_dict['\t'.join([ref_col, mut_col])]+=count
    mut_dict_sorted = sorted(list(six.iteritems(mut_dict)), key=lambda x: x[1], reverse=True)
    for i in range(len(mut_dict_sorted)):
        mut_dict_sorted[i]=(mut_dict_sorted[i][0]+'\t'+seq_dict[mut_dict_sorted[i][0]], mut_dict_sorted[i][1])
    return(mut_dict_sorted)

def parse_cigar(cigar, start, seq,ref_name, ref_dict): 
    """parses cigar string. M's are compared to reference sequence to identify 
    snps and indels are marked by position and length. optionally the alignment 
    can be restricted to certain width based on the cleavage site. snps can also 
    be ignored. also outputs the sequence and frame information; note if reads 
    are filtered, the frame flag is not updated"""
    Cigchar = re.findall('\D+', cigar)
    digits = [int(i) for i in re.findall('\d+', cigar)]
    refPos = start -1
    cigarPos = 0
    mut = []
    mutPos = []
    ref_seq=(ref_dict[ref_name]).upper()
    
    frame_flag=0
    seq_tracker=''
    for i in range(len(Cigchar)):
        if Cigchar[i] == 'S':
            cigarPos += digits[i]
        elif Cigchar[i] == 'M':
            parseRef = refPos + digits[i]
            parseCig = cigarPos + digits[i]
            for j in range(len(ref_seq[refPos:parseRef])):
                seq_tracker+=(seq[cigarPos:parseCig][j])
                if ref_seq[refPos:parseRef][j]!=seq[cigarPos:parseCig][j]:
                    mutPos.append(str(j+refPos+1))
                    mut.append(seq[cigarPos:parseCig][j])
            cigarPos+=digits[i]
            refPos+=digits[i]
        elif Cigchar[i]=='D':
            mutPos.append(str(refPos+1))
            mut.append(str(digits[i])+'D')
            frame_flag-=(digits[i])
            seq_tracker+='-'*(digits[i])  
            refPos+=digits[i]
        elif Cigchar[i]=='I':
            mutPos.append(str(refPos))
            mut.append(str(digits[i])+'I')
            seq_tracker+=('[' +seq[cigarPos:cigarPos+digits[i]] + ']')
            frame_flag+=digits[i]
            cigarPos+=digits[i]

    MutOut = ''
    frame_flag = abs(frame_flag % 3)
    if frame_flag != 0:
        frame = 'yes'
    else:
        frame = 'no' 
    if len(mut) > 0:
        Mut = ["%s:%s" % t for t in zip(mutPos, mut)]
        MutOut = ";".join(Mut)
    elif len(mut) == 0:
        MutOut = 'WT'
   
    return(MutOut, frame, seq_tracker)                     
                           
def setup_mut_caller(ref_fasta):
    """Converts fasta to dictionary of names/sequences. Supports multisequence 
    fasta but current pipeline implementation should only ever see single seq 
    fastas"""
    ref = open(ref_fasta, 'r')
    ref_lines = ref.readlines()
    for i in range(len(ref_lines)):
        ref_lines[i]=ref_lines[i].rstrip().strip('>')
        #print(ref_lines)
    ref_dict={}
    for i in range(len(ref_lines)//2):
        ref_dict[ref_lines[2*i]]=ref_lines[2*i+1]
    ref.close()
    return(ref_dict)
 
def single_sample_mut_parser(out_dir, q_cutoff, ref_dict, bam, parse_snps, noise_filter,restrict_alignment, seq_range):
    """parses bam file"""
    sample_name = re.sub('.bam', '', os.path.basename(bam))
    out_file = os.path.join(out_dir, sample_name + '.mut.xls')
    #parse_sam(bam, out_file, ref_dict, q_cutoff, noise_filter, parse_snps, restrict_alignment, upstream_restrict, downstream_restrict, N_bases)
    read_info = {}
    N_reads = 0
    bamfile=pysam.AlignmentFile(bam, 'rb')
    filter_pass_reads = 0
    for line in bamfile:
        N_reads += 1
        if line.cigarstring!=None and line.has_tag('AS')  and line.get_tag('AS') >= q_cutoff and line.seq!=None:
            filter_pass_reads +=1
            seq=line.seq
            PAstart = line.reference_start+1
            PAcigar = line.cigarstring
            PAend = line.reference_end
            refname = line.reference_name
            PAmuts, frame, seq_tracker = parse_cigar(PAcigar, PAstart, seq, refname, ref_dict)
            PArange = str(PAstart) + ':' + str(PAend)
            concat = '\t'.join([refname,PAmuts,frame,seq_tracker])
            if concat in read_info:
                read_info[concat]+=1
            else:
                read_info[concat]=1

    bamfile.close()
    ordered_out = sorted(list(six.iteritems(read_info)), key=lambda x: x[1], reverse=True)
    f_out=open(out_file, 'w')
    f_out.write('\t'.join(['Ref', 'Mut', 'Frameshift', 'N_mut', 'N_mapped', 'Mut_%', 'Alignment'])+'\n')

    if noise_filter>0.0 or not parse_snps:
        ordered_out = filter_reads(ordered_out, filter_pass_reads, parse_snps, noise_filter, ref_dict[refname])
    if restrict_alignment:
        ordered_out=parse_restrict_alignment(seq_range, ordered_out)
    for i in range(len(ordered_out)):
        out = ordered_out[i][0].split('\t')
        f_out.write('\t'.join([out[0], out[1], out[2], str(ordered_out[i][1]), str(filter_pass_reads), "%.2f" % (100.0*ordered_out[i][1]/filter_pass_reads), out[3]])+'\n')
    f_out.close()
    return


def get_mut_parser_params(parser_dict):
    """sets some values from parser_dict"""
    q_cutoff = 500 #static because needle's alignment score is esoteric
    parse_snps=parser_dict['mutation_calling']['parse_substitutions']
    noise_filter = parser_dict['mutation_calling']['variant_filter']
    noise_filter = noise_filter/100.0
    return(q_cutoff, parse_snps, noise_filter)


#if __name__ == '__main__': #not currently functional
#
#    parser = SafeConfigParser()
#    
#    ##parameter setup
#    pipeline_config = sys.argv[1]
#    ref_fasta = sys.argv[2]
#    amplicon_num = int(sys.argv[3])
#    bam = sys.argv[4]
#    parser.readfp(open(pipeline_config))
#
#    target, primer1, primer2, amplicon, sgRNA, n_amplicons=get_reference(parser.get('reference_files', 'amplicon_reference'))
#    target = target[amplicon_num]
#    sgRNA = sgRNA[amplicon_num]
#    out_dir = os.path.join(parser.get('i/o', 'output_directory'), target, 'mut')
#    
#    q_cutoff, parse_snps, noise_filter = get_mut_parser_params(parser)
#    restrict_alignment, upstream_restrict, downstream_restrict, N_bases=setup_restrict(parser, sgRNA,ref_fasta)
#    ref_dict=setup_mut_caller(ref_fasta)
#    single_sample_mut_parser(out_dir, q_cutoff, ref_dict, bam, parse_snps, noise_filter, restrict_alignment, upstream_restrict, downstream_restrict, N_bases)
  
