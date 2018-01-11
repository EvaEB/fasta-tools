import sys
import matplotlib.pyplot as plt
import numpy as np
from ete3 import Tree

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

def njTree(seqs=None,filename=None):
    '''creates and prints/plots a simple nj tree'''
    if seqs is None:
        seqs = read_fasta(filename,consensus)

    seq_len = len(seqs.values()[0])
    del seqs['consensus']

    diff = np.ones((len(seqs),len(seqs)))*(seq_len+1)
    for i in range(len(seqs)):
        for j in range(i+1,len(seqs)):
            diff[i,j] = sum( seqs.values()[i][k] != seqs.values()[j][k] for k in range(seq_len) )

    print diff
    newick = ''
    names = seqs.keys()

    while True:
        min_value = np.min(diff)
        if min_value == seq_len+1:
            break
        min_location = np.where(np.min(diff)==diff)
        min_location = zip(min_location[0],min_location[1])

        for i in min_location:
            if names[i[0]]!='' and names[i[1]]!='':
                names[i[0]] = '{},{}'.format(names[i[0]], names[i[1]])
                names[i[1]] = ''
                diff[i[1]:,i[0]] = (diff[i[1]:,i[1]]+diff[i[1]:,i[0]])/2
                diff[:,i[1]] = seq_len+1

                diff[i[0],i[1]:] = (diff[i[1],i[1]:]+diff[i[0],i[1]:])/2
                diff[i[1],:] = seq_len+1

        names = [add_brackets(i) for i in names]
        print diff
        print names
        raw_input()
    newick = '('+[i for i in names if i != ''][0]+');'
    #print newick
    t = Tree(newick)
    #print t
    t.show()

def add_brackets(txt):
    if ',' in txt:
        if txt[0] == '(' and txt[-1] == ')':
            counter = 0
            for i in range(len(txt)):
                if txt[i] == '(':
                    counter+=1
                elif txt[i] == ')':
                    counter-=1
                if counter == 0:
                    break
            if i+1 == len(txt):
                return txt
            else:
                return '({})'.format(txt)
        else:
            return '({})'.format(txt)
    else:
        return txt
if __name__ == '__main__':
    fasta_file = sys.argv[1]
    try:
        action = sys.argv[2]
    except IndexError:
        print 'no action given, assuming skyline'
        action = 'skyline'

    seqs = read_fasta(fasta_file,consensus=0)

    if action == 'skyline':
        skyline(seqs,consensus = 0)
    if action == 'tree':
        njTree(seqs)
