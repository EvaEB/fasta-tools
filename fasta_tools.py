import sys
import matplotlib.pyplot as plt
import numpy as np

def read_fasta(filename):
    '''
    reads in a fasta file as dictionary with sequence identifiers as keys and
    sequences as values
    '''
    seqs = {}
    with open(fasta_file) as f:
        for line in f.readlines():
            if '>' in line:
                current = line.strip()[1:]
                seqs[current] = ''
            else:
                seqs[current]+=line.strip()

    return seqs

def skyline(seqs=None,filename=None):
    '''
    creates a highlighter plot (as on hiv.lanl.gov) of already aligned sequences
    '''
    colors = {'A': '#bf6c60', 'T': '#99cc33', 'G': '#2a3326', 'C': '#36b8d9', '-': '#792080'}
    if seqs is None:
        seqs = read_fasta(filename)
    seq_len = len(seqs.values()[0])
    first = seqs.keys()[0]
    for i,idx in enumerate(seqs):
        plt.plot([0,seq_len],[i,i],color='grey')
        changes = np.where(np.array(list(seqs[idx])) != np.array(list(seqs[first])))
        for change in changes[0]:
            plt.scatter(change,i,color=colors[seqs[idx][change]])
    plt.axis('off')
    plt.show()


if __name__ == '__main__':
    fasta_file = sys.argv[1]
    seqs = read_fasta(fasta_file)
    skyline(seqs)
