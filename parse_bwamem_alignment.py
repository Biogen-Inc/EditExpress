import sys
import re
import os
import fnmatch
import subprocess
from copy import deepcopy
import pysam
import difflib
#from configparser import SafeConfigParser
from editexpress import natural_sortkey, get_reference, reverse_complement
import six
def get_frame(mut_string):
    """returns the frame (0, 1 or 2 shift) based on positions and length of indels"""
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


def filter_reads(ordered_list, read_cnt, threshold, parse_snps):
    """A 'noise' filter. parses through output and finds total occurences of 
    mutations at a particular site. if they are below a definable threshold, 
    then they will be "ignored"; in the sample level files, they will still 
    appear on separate lines with their sequences intact, but their reported 
    mutations will be altered. these will then be summed in the topSeqs output"""

    mut_dict = {}
    filter_muts = []
    filter_dict = {}
    filter_out = []
    for i in range(len(ordered_list)): # pass 1: grab mutations and determine in how many reads they appear
        parse = ordered_list[i][0].split('\t')
        ref_col = parse[0]
        mut_col = parse[1]
        count = ordered_list[i][1]
        mut_list = mut_col.split(';')
        for j in range(len(mut_list)):
            if mut_list[j] in mut_dict:
                mut_dict[mut_list[j]]+=count
            else:
                mut_dict[mut_list[j]]=count
    for key in mut_dict:
        if mut_dict[key]/float(read_cnt) < threshold and key!='exp_WT' or not parse_snps and any(i for i in ['A', 'C', 'G', 'T', 'N'] if i in key):
            filter_muts.append(key)
    
    mut_dict={}
    for i in range(len(ordered_list)): #pass 2: filter mutations failing threshold
        parse = ordered_list[i][0].split('\t')
        ref_col = parse[0]
        mut_col = parse[1]
        count = ordered_list[i][1]
        mut_list = mut_col.split(';')
        if any([k for k in mut_list if k not in filter_muts]):
            mut_col=";".join([k for k in mut_list if k not in filter_muts])
        else:
            mut_col='exp_WT'

        frame_flag = get_frame(mut_col)
        if (ref_col+'\t'+mut_col + '\t' + str(frame_flag)) in filter_dict:
            filter_dict[ref_col+'\t'+mut_col + '\t' + str(frame_flag)]+=count
        else:
            filter_dict[ref_col+'\t'+mut_col + '\t' + str(frame_flag)]=count
            
    filter_out = sorted(list(six.iteritems(filter_dict)), key=lambda x: x[1], reverse=True)
    return(filter_out)

def parse_cigar(cigar, start, seq,ref_name, ref_dict): 
    """parses cigar string. M's are compared to reference sequence to identify 
    snps and indels are marked by position and length. optionally the alignment 
    can be restricted to certain width based on the cleavage site. snps can 
    also be ignored. for big gap, the sequence is not output because it is 
    potentially quite long"""
    Cigchar = re.findall('\D+', cigar)
    digits = [int(i) for i in re.findall('\d+', cigar)]
    refPos = start -1
    cigarPos = 0
    mut = []
    mutPos = []
    ref_seq=ref_dict[ref_name]

    for i in range(len(Cigchar)):
        if Cigchar[i] == 'S':
            cigarPos += digits[i]
        elif Cigchar[i] == 'M':
            parseRef = refPos + digits[i]
            parseCig = cigarPos + digits[i]
            for j in range(len(ref_seq[refPos:parseRef])):        
                if ref_seq[refPos:parseRef][j]!=seq[cigarPos:parseCig][j]:
                    mutPos.append(str(j+refPos+1))
                    mut.append(seq[cigarPos:parseCig][j])
            cigarPos+=digits[i]
            refPos+=digits[i]
        elif Cigchar[i]=='D':
            mutPos.append(str(refPos+1))
            mut.append(str(digits[i])+'D')
            refPos+=digits[i]
        elif Cigchar[i]=='I':
            mutPos.append(str(refPos))
            mut.append(str(digits[i])+'I')
            cigarPos+=digits[i]
    MutOut = ''
    if len(mut) > 0:
        Mut = ["%s:%s" % t for t in zip(mutPos, mut)]
        MutOut = ";".join(Mut)
    return(MutOut)                       

#def parse_sam(bam, out_file, ref_dict, mapq, noise_filter, parse_snps):
#    """Loops through bam file parses cigar string against sequence to obtain full
#    mutation string. Supplementary alignment tags ('SA' fields) are interepreted
#    to be gapped sequences and combined with the primary alignment; the gap is
#    the distance between PA and SA"""
#    read_info = {}
#    N_reads = 0
#    bamfile=pysam.AlignmentFile(bam, 'rb')
#    for line in bamfile:
#        if line.cigarstring!=None and "H" not in line.cigarstring and line.mapping_quality >= mapq: #basically equivalent to filtering out flags 4 and 256/2048, unmapped and supplementary/secondary (samtools view -F 2308 -q mapq should be identical)
#            qname = line.query_name
#            seq=line.seq
#            PAstart = line.reference_start+1
#            PAcigar = line.cigarstring
#            PAend = line.reference_end
#            refname = line.reference_name
#            PAmuts = parse_cigar(PAcigar, PAstart, seq, refname, ref_dict)
#
#            if qname not in read_info:
#                N_reads += 1
#                read_info[qname] = 'parsing'
#            if line.has_tag("SA"):
#                SA = line.get_tag("SA").split(",")
#                SA_q = SA[4] 
#                if SA_q >= mapq: 
#                    SA = line.get_tag("SA").split(",")
#                    SAstart = int(SA[1])
#                    SAcigar = SA[3]
#                    SAdigits = [int(i) for i in re.findall('\d+',SAcigar)]
#                    SAchar = re.findall('\D+', SAcigar)
#                    SAend = SAstart
#                    for i in range(len(SAchar)):
#                        if SAchar[i] == 'M' or SAchar[i]=='D':
#                            SAend+=int(SAdigits[i])
#                    SAmuts = parse_cigar(SAcigar, SAstart, seq,refname, ref_dict)
#                    PAmuts = parse_cigar(PAcigar, PAstart, seq,refname, ref_dict)
#                    SArange = str(SAstart)+':'+str(SAend)
#                    PArange = str(PAstart) + ':' + str(PAend)
#                    mut_pos = []
#                    mut_char = []
#                    if SAstart>PAend:
#                        Gap = str(SAstart-PAend-1)
#                        if line.is_read1:
#                            Gap+='R1'
#                        else:
#                            Gap+='R2'
#                        all_muts = str(PAend+1)+":"+Gap+"D"
#                        if len(PAmuts)>0:
#                            all_muts = PAmuts + ";" + all_muts 
#                        if len(SAmuts)>0: 
#                            all_muts = all_muts + ";" + SAmuts
#                        mut_pos = [int(i.strip(':')) for i in re.findall('\d+:', all_muts)]
#                        mut_char = re.findall('\D+', all_muts)
#                    else:
#                        Gap = str(PAstart-SAend)
#                        if line.is_read1:
#                            Gap+='R1'
#                        else:
#                            Gap+='R2'
#                        all_muts = str(SAend)+":"+Gap+"D"
#                        if len(SAmuts)>0:
#                            all_muts = SAmuts + ";" + all_muts
#                        if len(PAmuts)>0:
#                            all_muts = all_muts + ";" + PAmuts
#                        mut_pos = [int(i.strip(':')) for i in re.findall('\d+:', all_muts)]
#                        mut_char = re.findall('\D+', all_muts)
#
#                    if read_info[qname]=='parsing' or read_info[qname]==refname + '\t' + 'exp_WT'+'\t'+str(0):
#                        read_info[qname]=refname + '\t' + all_muts + '\t' + str(get_frame(all_muts))
#                    else:
#                        temp_muts = read_info[qname].split('\t')[1]
#                        temp_muts = temp_muts+';'+all_muts
#                        temp_muts = list(set(temp_muts.split(';')))
#                        R_indices = [i for i, elem in enumerate([x for x in temp_muts]) if 'R' in elem] 
#                        if len(R_indices)==2:
#                            for i in range(len(R_indices)):
#                                if 'R1' in temp_muts[R_indices[i]]:
#                                    temp_muts[R_indices[i]]=re.sub('R1', '', temp_muts[R_indices[i]])
#                                elif 'R2' in temp_muts[R_indices[i]]:
#                                    R2=R_indices[i]
#                            temp_muts.pop(R2)
#                        elif len(R_indices)==1:
#                            temp_muts[R_indices[0]]=re.sub('R\d', '', temp_muts[R_indices[0]])
#                        temp_muts = ";".join(sorted(temp_muts, key=natural_sortkey))
#                        frame = get_frame(temp_muts)
#                        read_info[qname]=refname+'\t'+temp_muts + '\t' + str(frame)
#            else:
#                if PAmuts =='':
#                    PAmuts = 'exp_WT'
#                if read_info[qname]=='parsing' or read_info[qname]==refname+'\t'+'exp_WT'+'\t'+str(0):
#                    read_info[qname]=refname + '\t' + PAmuts + '\t' + str(get_frame(PAmuts))
#                elif PAmuts!='exp_WT':
#                    temp_muts = read_info[qname].split('\t')[1]
#                    temp_frame = int(read_info[qname].split('\t')[2])
#                    temp_muts = temp_muts+';'+PAmuts
#                    temp_muts = list(set(temp_muts.split(';')))
#                    R_indices = [i for i, elem in enumerate([x for x in temp_muts]) if 'R' in elem]
#                    if len(R_indices)==2:
#                        for i in range(len(R_indices)):
#                            if 'R1' in temp_muts[R_indices[i]]:
#                                temp_muts[R_indices[i]]=re.sub('R\d', '', temp_muts[R_indices[i]])
#                            elif 'R2' in temp_muts[R_indices[i]]:
#                                R2=R_indices[i]
#                        temp_muts.pop(R2)
#                    elif len(R_indices)==1:
#                       temp_muts[R_indices[0]]=re.sub('R\d', '', temp_muts[R_indices[0]])
#                    temp_muts = ";".join(sorted(temp_muts, key=natural_sortkey))
#                    frame = get_frame(temp_muts)
#                    read_info[qname]=refname+'\t'+temp_muts + '\t' + str(frame%3)
#
#    bamfile.close()
#    read_info2 = {} #removing read names
#    for key in read_info:
#        if read_info[key]!='parsing':
#            temp = read_info[key].split('\t')
#            if 'R' in temp[1]:
#                temp[1] = re.sub('R\d', '', temp[1])
#            key2 = temp[0] + '\t' + temp[1] + '\t' + temp[2]
#            if key2 in read_info2:
#                read_info2[key2]+=1
#            else:
#                read_info2[key2]=1
#    
#    ordered_out = sorted(read_info2.items(), key=lambda x: x[1], reverse=True)
#    
#    f_out = open(out_file, 'w')
#    f_out.write("Ref"+'\t' + 'Mut' + '\t' + 'Frameshift' + '\t' + 'N_mut' + '\t' + 'N_mapped' + '\t' + 'Mut_%'  + "\n")
#    if noise_filter>0 or not parse_snps:
#        ordered_out = filter_reads(ordered_out, N_reads, noise_filter, parse_snps)
#    for i in range(len(ordered_out)):
#        out = ordered_out[i][0].split('\t')
#        if int(out[2])!=0:
#            frame_out = 'yes'
#        else:
#            frame_out = 'no'
#        f_out.write(out[0]+'\t'+ out[1] + '\t' + frame_out + '\t' + str(ordered_out[i][1]) + '\t' + str(N_reads) + '\t' + "%.2f" % (100.0*ordered_out[i][1]/N_reads) + '\n')
#    f_out.close()

def setup(ref_fasta):
    """Converts fasta to dictionary of names/sequences. Supports multisequence fasta
    but current pipeline implementation should only ever see single seq fastas"""
    ref = open(ref_fasta, 'r')
    ref_lines = ref.readlines()
    for i in range(len(ref_lines)):
        ref_lines[i]=ref_lines[i].rstrip().strip('>')
    ref_dict={}
    for i in range(len(ref_lines)//2):
        ref_dict[ref_lines[2*i]]=ref_lines[2*i+1]
    ref.close()
    return(ref_dict)

def get_mut_parser_params(parser_dict):
    """grabs some parameters for"""
    mapq = parser_dict['mutation_calling']['bwa_min_mapq']
    parse_snps=parser_dict['mutation_calling']['parse_substitutions']
    noise_filter=parser_dict['mutation_calling']['variant_filter']
    noise_filter = noise_filter/100.0
    return(mapq, parse_snps, noise_filter)

def single_sample_mut_parser(out_dir, ref_dict, bam, mapq, parse_snps, noise_filter):
    """sets up sample name and outfile, runs the parsing function"""
    sample_name = re.sub('.bam', '', os.path.basename(bam))
    out_file = os.path.join(out_dir, sample_name + '.mut.xls')
    read_info = {}
    N_reads = 0
    bamfile=pysam.AlignmentFile(bam, 'rb')
    for line in bamfile:
        if line.cigarstring!=None and "H" not in line.cigarstring and line.mapping_quality >= mapq: #basically equivalent to filtering out flags 4 and 256/2048, unmapped and supplementary/secondary (samtools view -F 2308 -q mapq should be identical)
            qname = line.query_name
            seq=line.seq
            PAstart = line.reference_start+1
            PAcigar = line.cigarstring
            PAend = line.reference_end
            refname = line.reference_name
            PAmuts = parse_cigar(PAcigar, PAstart, seq, refname, ref_dict)

            if qname not in read_info:
                N_reads += 1
                read_info[qname] = 'parsing'
            if line.has_tag("SA"):
                SA = line.get_tag("SA").split(",")
                SA_q = int(SA[4])
                if SA_q >= mapq: 
                    SA = line.get_tag("SA").split(",")
                    SAstart = int(SA[1])
                    SAcigar = SA[3]
                    SAdigits = [int(i) for i in re.findall('\d+',SAcigar)]
                    SAchar = re.findall('\D+', SAcigar)
                    SAend = SAstart
                    for i in range(len(SAchar)):
                        if SAchar[i] == 'M' or SAchar[i]=='D':
                            SAend+=int(SAdigits[i])
                    SAmuts = parse_cigar(SAcigar, SAstart, seq,refname, ref_dict)
                    PAmuts = parse_cigar(PAcigar, PAstart, seq,refname, ref_dict)
                    SArange = str(SAstart)+':'+str(SAend)
                    PArange = str(PAstart) + ':' + str(PAend)
                    mut_pos = []
                    mut_char = []
                    if SAstart>PAend:
                        Gap = str(SAstart-PAend-1)
                        if line.is_read1:
                            Gap+='R1'
                        else:
                            Gap+='R2'
                        all_muts = str(PAend+1)+":"+Gap+"D"
                        if len(PAmuts)>0:
                            all_muts = PAmuts + ";" + all_muts 
                        if len(SAmuts)>0: 
                            all_muts = all_muts + ";" + SAmuts
                        mut_pos = [int(i.strip(':')) for i in re.findall('\d+:', all_muts)]
                        mut_char = re.findall('\D+', all_muts)
                    else:
                        Gap = str(PAstart-SAend)
                        if line.is_read1:
                            Gap+='R1'
                        else:
                            Gap+='R2'
                        all_muts = str(SAend)+":"+Gap+"D"
                        if len(SAmuts)>0:
                            all_muts = SAmuts + ";" + all_muts
                        if len(PAmuts)>0:
                            all_muts = all_muts + ";" + PAmuts
                        mut_pos = [int(i.strip(':')) for i in re.findall('\d+:', all_muts)]
                        mut_char = re.findall('\D+', all_muts)

                    if read_info[qname]=='parsing' or read_info[qname]==refname + '\t' + 'exp_WT'+'\t'+str(0):
                        read_info[qname]=refname + '\t' + all_muts + '\t' + str(get_frame(all_muts))
                    else:
                        temp_muts = read_info[qname].split('\t')[1]
                        temp_muts = temp_muts+';'+all_muts
                        temp_muts = list(set(temp_muts.split(';')))
                        R_indices = [i for i, elem in enumerate([x for x in temp_muts]) if 'R' in elem] 
                        if len(R_indices)==2:
                            for i in range(len(R_indices)):
                                if 'R1' in temp_muts[R_indices[i]]:
                                    temp_muts[R_indices[i]]=re.sub('R1', '', temp_muts[R_indices[i]])
                                elif 'R2' in temp_muts[R_indices[i]]:
                                    R2=R_indices[i]
                            temp_muts.pop(R2)
                        elif len(R_indices)==1:
                            temp_muts[R_indices[0]]=re.sub('R\d', '', temp_muts[R_indices[0]])
                        temp_muts = ";".join(sorted(temp_muts, key=natural_sortkey))
                        frame = get_frame(temp_muts)
                        read_info[qname]=refname+'\t'+temp_muts + '\t' + str(frame)
            else:
                if PAmuts =='':
                    PAmuts = 'exp_WT'
                if read_info[qname]=='parsing' or read_info[qname]==refname+'\t'+'exp_WT'+'\t'+str(0):
                    read_info[qname]=refname + '\t' + PAmuts + '\t' + str(get_frame(PAmuts))
                elif PAmuts!='exp_WT':
                    temp_muts = read_info[qname].split('\t')[1]
                    temp_frame = int(read_info[qname].split('\t')[2])
                    temp_muts = temp_muts+';'+PAmuts
                    temp_muts = list(set(temp_muts.split(';')))
                    R_indices = [i for i, elem in enumerate([x for x in temp_muts]) if 'R' in elem]
                    if len(R_indices)==2:
                        for i in range(len(R_indices)):
                            if 'R1' in temp_muts[R_indices[i]]:
                                temp_muts[R_indices[i]]=re.sub('R\d', '', temp_muts[R_indices[i]])
                            elif 'R2' in temp_muts[R_indices[i]]:
                                R2=R_indices[i]
                        temp_muts.pop(R2)
                    elif len(R_indices)==1:
                       temp_muts[R_indices[0]]=re.sub('R\d', '', temp_muts[R_indices[0]])
                    temp_muts = ";".join(sorted(temp_muts, key=natural_sortkey))
                    frame = get_frame(temp_muts)
                    read_info[qname]=refname+'\t'+temp_muts + '\t' + str(frame%3)

    bamfile.close()
    read_info2 = {} #removing read names
    for key in read_info:
        if read_info[key]!='parsing':
            temp = read_info[key].split('\t')
            if 'R' in temp[1]:
                temp[1] = re.sub('R\d', '', temp[1])
            key2 = temp[0] + '\t' + temp[1] + '\t' + temp[2]
            if key2 in read_info2:
                read_info2[key2]+=1
            else:
                read_info2[key2]=1
    
    ordered_out = sorted(list(six.iteritems(read_info2)), key=lambda x: x[1], reverse=True)
    
    f_out = open(out_file, 'w')
    f_out.write("Ref"+'\t' + 'Mut' + '\t' + 'Frameshift' + '\t' + 'N_mut' + '\t' + 'N_mapped' + '\t' + 'Mut_%'  + "\n")
    if noise_filter>0 or not parse_snps:
        ordered_out = filter_reads(ordered_out, N_reads, noise_filter, parse_snps)
    for i in range(len(ordered_out)):
        out = ordered_out[i][0].split('\t')
        if int(out[2])!=0:
            frame_out = 'yes'
        else:
            frame_out = 'no'
        f_out.write(out[0]+'\t'+ out[1] + '\t' + frame_out + '\t' + str(ordered_out[i][1]) + '\t' + str(N_reads) + '\t' + "%.2f" % (100.0*ordered_out[i][1]/N_reads) + '\n')
    f_out.close()

#if __name__ == '__main__':
#    
#    parser = SafeConfigParser()
#    runtime_conf = sys.argv[1]
#    ref_fasta = sys.argv[2]
#    target = sys.argv[3]
#    bam = sys.argv[4]
#    parser.readfp(open(runtime_conf))
#    print(runtime_conf, ref_fasta, target, bam)
#    #runtime parameters
#    out_dir = os.path.join(parser.get('i/o', 'output_directory'), target, 'mut')
#    mapq, parse_snps, noise_filter = get_mut_parser_params(parser)
#    ref_dict = setup(ref_fasta)
#    single_sample_mut_parser(out_dir, ref_dict, bam, mapq, parse_snps, noise_filter)

