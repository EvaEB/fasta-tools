import sys
import matplotlib.pyplot as plt
import numpy as np

def get_consensus(seqs):
    consensus = ''
    seqs = seqs.values()
    for i in range(len(seqs[0])):
        letters = [j[i] for j in seqs]
        consensus+= max(set(letters), key=letters.count)
    return consensus


def read_fasta(filename,consensus=-1):
    '''
    reads in a fasta file as dictionary with sequence identifiers as keys and
    sequences as values
    '''
    seqs = {}
    counter = 0
    with open(fasta_file) as f:
        for line in f.readlines():
            if '>' in line:
                if counter == consensus:
                    seqs['consensus'] = ''
                    current = 'consensus'
                else:
                    current = line.strip()[1:]
                    seqs[current] = ''
                counter+=1
            else:
                seqs[current]+=line.strip()
    if consensus == -1:
        seqs['consensus'] = get_consensus(seqs)
    return seqs

def skyline(seqs=None,filename=None,consensus=-1):
    '''
    creates a highlighter plot (as on hiv.lanl.gov) of already aligned sequences
    '''
    colors = {'A': '#bf6c60', 'T': '#99cc33', 'G': '#2a3326', 'C': '#36b8d9', '-': '#792080'}
    if seqs is None:
        seqs = read_fasta(filename,consensus)
    to_plot = seqs.keys()
    seq_len = len(seqs[seqs.keys()[0]])
    current = 'consensus'
    pos = 0
    mins=None
    while len(to_plot) > 0:
        plt.plot([0,seq_len],[pos,pos],color='grey')
        plt.text(seq_len,pos,current,va='center',ha='left')
        changes = [p for p,a,b in zip(range(seq_len),seqs['consensus'], seqs[current]) if a != b]
        to_plot.remove(current)
        for change in changes:
            plt.scatter(change,pos,color=colors[seqs[current][change]])
        mins = None
        for i in to_plot:
            changes = sum(1 for a, b in zip(seqs[current], seqs[i]) if a != b)
            if changes == 0:
                next_current = i
                mins = 0
                break
            elif mins is None or changes < mins:
                mins = changes
                next_current = i
        current = next_current
        pos -= 1

    plt.axis('off')
    plt.show()


if __name__ == '__main__':
    fasta_file = sys.argv[1]
    try:
        action = sys.argv[2]
    except IndexError:
        print 'no action given, assuming skyline'
        action = 'skyline'

    seqs = read_fasta(fasta_file,consensus=0)

    if action == 'skyline':
        skyline(seqs)
