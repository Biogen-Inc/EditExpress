import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import OrderedDict
import re
import six
#from configparser import SafeConfigParser
def pie_chart(track, total, outfile, colors):
    """creates pie chart using matplotlib"""
    track = list(six.iteritems(track))
    names = []
    sizes = []
    for key, value  in track:
        names.append(key)
        sizes.append(int(value))
     
    sizes = [float(x)/sum(sizes) for x in sizes]
    labels = ['{0}: {1:1.2f}%'.format(n,m*100) for n,m in zip(names, sizes)]
    fig, ax = plt.subplots()
    ax.spines['left'].set_position('zero')
    wedges, text = ax.pie(sizes, shadow=False, startangle=75, colors=colors)

    for w in range(len(wedges)):
        if sizes[w]>0:
            wedges[w].set_linewidth(0.5)
            wedges[w].set_edgecolor('black')
    ax.axis('equal')
    plt.legend(loc = 'upper left', labels=labels, bbox_to_anchor = (0.778,1), fontsize = 9)
    plt.subplots_adjust(left=0, bottom=0)
    plt.savefig(outfile)
    plt.close()
def parse_sample_level(infile):
    """grabs mutation and frame information from sample level mutation table"""
    with open(infile, 'r') as mut_file:
        next(mut_file)
        mut_track = OrderedDict([('unmodified',0), ('substitution',0), ('deletion',0), ('insertion',0), ('mixed',0)])
        frame_track = OrderedDict([('unmodified',0), ('no frameshift',0), ('frameshift',0)])
        for line in mut_file:
            parse=line.split('\t')
            mut = parse[1]
            frame = parse[2]
            if frame=='no':
                frame='no frameshift'
            elif frame=='yes':
                frame='frameshift'
            if mut == 'WT' or mut=='exp_WT':
                frame = 'unmodified'
                mut = 'unmodified'
            elif 'D' in mut and not any(x in mut for x in ['A', 'G', 'C', 'T']) and not 'I' in mut:
                mut = 'deletion'
            elif 'I' in mut and not any(x in mut for x in ['A', 'G', 'C', 'T']) and not 'D' in mut:
                mut = 'insertion'
            elif any(x in mut for x in ['A', 'G', 'C', 'T']) and not 'D' in mut and not 'I' in mut:
                mut = 'substitution'
            else:
                mut = 'mixed'
            n_hits = int(parse[3])
            total_mapped = int(parse[4])
            cov = float(parse[5])
            #if cov>=0.5: 
            mut_track[mut] += int(n_hits)
            frame_track[frame] += int(n_hits)
        pie_chart(mut_track, total_mapped, re.sub('.mut.xls', '.mut_piechart.png', infile), ['green',  'dodgerblue','gold','purple', 'orangered'])
        pie_chart(frame_track, total_mapped, re.sub('.mut.xls', '.frame_piechart.png', infile), ['cornflowerblue','orange', 'limegreen'])
                
if __name__ == '__main__': 
    
    in_file = sys.argv[1]
    parse_sample_level(in_file)

