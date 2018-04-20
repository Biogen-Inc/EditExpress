import sys
import re
import os
import subprocess
import string
import numpy as np
import math
import svgwrite
from editexpress import natural_sortkey, get_reference, reverse_complement
import six
from six.moves import range
colors = {'G': '#1CAE32', 'A': '#FACC21', 'T': '#00D8FF', 'C': '#E50202', '-': '#B9B9B9', 'N':'#000000'}

def visualize(cut_site, ref_seq, muts, title, out_file, additional_offset, upstream, downstream):
    """visualize alignment plot using the svgwrite package. credit to the GuideSeq
    package, this is built up from one of their functions"""
    ref_seq = ref_seq.upper()
    box_size = 15 
    x_max = 0
    y_max = 0
    dwg = svgwrite.Drawing(out_file, profile='full')
    cov_range_offset = 7*additional_offset #provides adjustment for the coverage range/percent coverage string 
    x_offset = 20
    y_offset = 50
    dwg.add(dwg.text(title, insert=(x_offset+cov_range_offset, 30), style="font-size:20px; font-family:Courier New"))
    for j, seq in enumerate(muts): #adds coverage/coverage range
        y = y_offset+10 + j*1.5 * box_size+25
        dwg.add(dwg.text(seq['cov'], insert=(x_offset-5, y), fill='black', style="font-size:10px; font-family:Courier New"))
    x_offset+=cov_range_offset
    tick_locations = [1, len(ref_seq)]
    tick_locations += list(range(len(ref_seq) + 1))[::5][1:]
    for x in tick_locations: #adds ticks. a bit weird if the upper bound isn't even
        dwg.add(dwg.text(str(x+upstream), insert=(x_offset + (x - 1) * box_size, y_offset - 2), style="font-size:10px; font-family:Courier New"))
    for i, c in enumerate(ref_seq): #add ref_seq row
        y = y_offset
        x = x_offset+ i * box_size
        dwg.add(dwg.rect((x, y), (box_size, box_size), fill=colors[c]))
        dwg.add(dwg.text(c, insert=(x + 3, y + box_size - 3), fill='black', style="font-size:15px; font-family:Courier New"))
    dwg.add(dwg.text('Reads', insert=(x_offset + box_size * len(ref_seq) + 16, y_offset + box_size - 3), style="font-size:15px; font-family:Courier New"))
    y_offset += 10  #extra spacing
    ins_out = []
    for j, seq in enumerate(muts): #loops through and draws dots for ref seq matches, letters for mismatches, - for deletions and inserts a plus above the sequence if an insertion begins there
        y = y_offset + (1.5*j * box_size)
        seq['seq'] = (seq['seq']).upper()
        for i, c in enumerate(seq['seq']):
            x = x_offset + i * box_size
            if c == ref_seq[i] or ref_seq[i] == 'N':
                dwg.add(dwg.text(u"\u2022", insert=(x + 4.5, 2 * box_size + y - 4), fill='black', style="font-size:10px; font-family:Courier New"))
                if y_max<2 * box_size + y - 4:
                    y_max = 2 * box_size + y - 4
            elif c=='N' and ref_seq[i] !='N':
                dwg.add(dwg.rect((x, box_size + y), (box_size, box_size), fill=colors[c]))
                dwg.add(dwg.text(c, insert=(x + 3, 2 * box_size + y - 3), fill='black', style="font-size:15px; font-family:Courier New"))
                if y_max<2 * box_size + y - 3:
                    y_max = 2 * box_size + y - 3
            else:
                dwg.add(dwg.rect((x, box_size + y), (box_size, box_size), fill=colors[c]))
                dwg.add(dwg.text(c, insert=(x + 3, 2 * box_size + y - 3), fill='black', style="font-size:15px; font-family:Courier New"))
                if y_max<2 * box_size + y - 3:
                    y_max = 2 * box_size + y - 3
            if len(seq['ins'])>0:
                temp=[ins_seq for ins_seq in seq['ins'] if int(ins_seq[0])==i+upstream]
                ins_tmp={}
                if len(temp)>0:
                    dwg.add(dwg.text(u"\u2795", insert=(x-11, 2*box_size+y-15), fill='black', style="font-size:9px; font-family:Courier New"))
                    if len(ins_tmp)==0:
                        ins_tmp['Mut'+str(j+1)]=[m+n+'+'+o for m, n, o in temp]
                    else:
                        ins_tmp['Mut'+str(j+1)]=ins_tmp['Mut'+str(j+1)] + '; ' + [m+n+'+'+o for m, n, o in temp]
                    ins_out.append(['Mut'+str(j+1),ins_tmp['Mut'+str(j+1)]])
        reads_text_pos = cov_range_offset+box_size * (len(ref_seq) + 1) + 20
        reads_text = dwg.text(str(seq['muts']), insert=(reads_text_pos, y+27), fill='black', style="font-size:15px; font-family:Courier New")
        dwg.add(reads_text)
        if x_max<reads_text_pos:
            x_max=reads_text_pos
    for n in cut_site:
        pos = 1.0*(n-upstream)
        cut_site_line = dwg.rect((pos*box_size+x_offset-1, y_offset-10), (2, y-20), fill='black') #Line indicating the specified cleavage site (directly in the middle of the plot). Actually a thin rectangle
        dwg.add(cut_site_line)
    if len(ins_out)>0:
        dwg.add(dwg.text("Insertion Details", insert = (x_offset, y+50), style="font-size:18px; font-family:Courier New"))
        for j in range(len(ins_out)):
            ins_text = dwg.text(ins_out[j][0] + ':\t'+ '; '.join(ins_out[j][1]), insert = (x_offset, y+65+j*15), style="font-size:15px; font-family:Courier New")  
            dwg.add(ins_text)
        if y_max<y+65+(len(ins_out)-1)*15:
            y_max = y+65+(len(ins_out)-1)*15
    if y_max == 0:
        y_max = 30 
    dwg.stretch()
    dwg.viewbox(0,0,x_max*1.1, y_max*1.1)
    dwg.fit(horiz='center', vert='middle', scale='meet')
    dwg.save()
    return
def parse_topseqs(infile):
    """finds all unique occurences of sequences in a topSeqs file, grabs their mutation
    string and returns a range of their coverage if multiple samples had that
    mutation. A heatmap would be better but I'll save that for if I swap this
    plotting implementation for a matplotlib one"""
    with open(infile, 'r') as topSeqs:
        next(topSeqs)  
        track = {}
        mut_range = {}
        for line in topSeqs:
            parse = line.rstrip().split('\t')
            n_seqs = (len(parse)-5)//5
            for j in range(n_seqs):
                mut=parse[(j+1)*5]
                cov_percent=round(float(parse[(j+1)*5+2]),1)
                seq = parse[(j+1)*5+4].rstrip()
                out = mut+'\t'+seq
                if out in track:
                    track[out]+=1
                    if cov_percent<float(mut_range[out][0]):
                        mut_range[out][0] = str(cov_percent)
                    elif cov_percent>float(mut_range[out][1]):
                        mut_range[out][1] = str(cov_percent)
                else:
                    track[out]=1
                    mut_range[out] = [str(cov_percent), str(cov_percent)]
    return(track, mut_range)

def parse_sample(infile, summary_threshold):
    """for sample level files, grabs coverage and sequence for variants above
    summary threshold"""
    with open(infile, 'r') as mut_file:
        next(mut_file)
        track = {}
        mut_range = {}
        for line in mut_file:
            parse=line.split('\t')
            mut = parse[1]
            cov = round(float(parse[5]),1)
            seq = parse[6]
            out = mut+'\t'+seq
            if cov >= summary_threshold:
                track[out] = cov
    return(track)

def parse_mut(parse_mut, parse_seq):
    """separates out insertions from mutated sequences to make plotting simpler.
    insertions are noted by symbols above their position in the plot and noted 
    below"""
    parse_digits=[]
    parse_chars=[]
    parse_I=[]
    I=''
    if 'I' in parse_mut:
        pattern = r'\[\w+\]'
        noI=re.split(pattern, parse_seq)
        I=re.findall(pattern, parse_seq)
        noI_num = [len(i) for i in noI]
        for i in range(len(I)):
            I[i]=re.sub('\[', '', I[i])
            I[i]=re.sub('\]', '', I[i])
            parse_I.append((str(noI_num[i]), ':', str(I[i])))
        parse_seq=''.join(noI)

    temp=parse_mut.split(';')
    for i in temp:
        parse_digits.append(i.split(":")[0])
        parse_chars.append(i.split(":")[1])

    for i in sorted(list(range(len(parse_digits))), reverse=True): #substitutions may clutter, so they are plotted but not noted on the xaxis. a mut sequence with only subs will be annotated as "substitutions"
        if parse_chars[i] in ['A', 'C', 'G', 'T', 'N']:
            del(parse_digits[i])
            del(parse_chars[i])
    parse_mut=";".join([str(m)+":"+n for m,n in zip(parse_digits,parse_chars)])
    return(parse_seq, parse_mut,parse_I)

def setup_alignment_plot(parser_dict, ref_fasta, coords, ref_num, sgRNA, log):
    """sets up various parameters and performs checks to see if alignment
    restriction (narrowing the alignment output to N bases around the cleavage
    site, as determined by guide rnas or directly)"""
    #changing plot area from N_bases to a simple range, so that multiple ranges are easier to support...I hope
    #this enables empty zero_site! beware
    
    default=False
    seq_range=coords.split(',')[ref_num]
    seq_range = seq_range.split('-')
    if not len(seq_range)==2:
        default=True
    try:
        seq_range=[int(i) for i in seq_range]
    except ValueError:
        default=True
    ref = open(ref_fasta, 'r') 
    ref_lines = ref.readlines()
    ref_lines = ref_lines[1].rstrip().upper()
    amplicon = re.sub('>', '', ref_lines[0])
    ref.close()

    if not default:
        try:
            if any([number<0 for number in seq_range]):
                default=True
            else:
                seq_range = (min(seq_range), max(seq_range))
        except KeyError:
            default=True
    if default:
        amp_length=len(ref_lines)
        center=amp_length//2
        if amp_length*0.33>40:
            seq_range=((center-40),(center+40))
        else:
            seq_range=((center-int(amp_length*0.33)),(center+int(amp_length*0.33)))

    zero_site = []
    cleave_sites=''
    
    try:
        cleave_sites = parser_dict['alignment_plot']['abs_cleave_sites']
        separated=re.findall('\(.*?\)', cleave_sites)
        if len(separated)==0:
            separated = [cleave_sites]
        separated = [re.sub(r'([()])', '', i) for i in separated]
        separated = [re.split(' |,',i) for i in separated]      
        try:
            for j in separated[ref_num]:
                zero_site.append(int(j))
        except ValueError:
            log.warning('Improperly formatted absolute cleave site values')
            zero_site=[]
    except KeyError:
        pass
    if zero_site==[]:
        try:
            rel_cleave_sites = parser_dict['alignment_plot']['rel_cleave_sites'].split(',')
            rel_cleave_site=int(rel_cleave_sites[ref_num])
            if len(sgRNA)==0:
                log.warning('No guide RNA sequences given!')
                log.warning('Proceding without annotating cut sites')
            else:
                for j in sgRNA:
                    try:
                        if rel_cleave_site>0:
                            zero_site.append(ref_lines.index(j, 0) + rel_cleave_site)
                        elif rel_cleave_site<0:
                            zero_site.append(ref_lines.index(j, 0) + len(j) + rel_cleave_site)
                        else:
                            zero_site.append(ref_lines.index(j, 0))
                            log.warning('0 is an uninformative value for "rel_cleave_sites"')
                            log.warning("Assuming that cleavage occurs exactly at 5' end of guide... ")
                    except (ValueError, TypeError):
                        log.error(j + ' not found in reference sequence!')
   
        except (KeyError, ValueError):
                log.warning('Missing or incorrect absolute and relative cut options in config file')
                log.warning('Proceding without annotating cut sites')
    return(zero_site, seq_range, ref_lines, default)

def pre_alignment_plot(out_dir, in_file, in_file_type, parser_dict, seq_range, default, zero_site, ref_lines, target, log):
    """sets up all values for and feeds into the alignent plot """
    if in_file_type == 'topSeqs':
        track, mut_range=parse_topseqs(in_file)
        if len(track) == 0:
            log.error('No sequences above coverage threshold!')
            return
        ordered_track = sorted(list(six.iteritems(track)), key=lambda x: x[1], reverse=True)
        outfile = os.path.join(out_dir, target, 'mut', target + '_topSeqs_alignments.svg')
        target = target + ' topSeqs'
    else:
        track = parse_sample(in_file, parser_dict['mutation_calling']['summary_threshold'])
        ordered_track = sorted(list(six.iteritems(track)), key=lambda x: x[1], reverse=True)
        outfile = re.sub('.xls', '.svg', in_file)
        target = target + ' ' + os.path.basename(re.sub('.xls', '', in_file))
    ref_lines_sub = ref_lines[seq_range[0]:seq_range[1]]
    mut_list=[]
    cov_range_len = []
    for i in ordered_track:
        mut,seq= i[0].split('\t')
        mut=mut.rstrip()
        seq=seq.rstrip()
        if in_file_type == 'topSeqs':
            if mut_range[i[0]][1]==mut_range[i[0]][0]:
                cov_range = mut_range[i[0]][0]
            else:
                cov_range = "-".join(mut_range[i[0]])
            cov_range_len.append(len(cov_range))
        else:
            cov_range = str(i[1])
            cov_range_len.append(len(cov_range))
        if mut =='WT':
            if not (parser_dict['mutation_calling']['restrict_alignment'] and not default):
                seq = seq[seq_range[0]:seq_range[1]]
            mut_list.append({'muts':'WT', 'seq':seq, 'cov':cov_range+'%','ins':[]})
        else:
            seq, mut, ins = parse_mut(mut, seq)
            if mut == '':
                mut = 'snv'
            if not (parser_dict['mutation_calling']['restrict_alignment'] and not default):
                seq = seq[seq_range[0]:seq_range[1]]
            mut_list.append({'muts':mut, 'seq':seq, 'cov':cov_range+'%', 'ins':ins})
    visualize(zero_site, ref_lines_sub, mut_list, target, outfile, max(cov_range_len), seq_range[0], seq_range[1])

