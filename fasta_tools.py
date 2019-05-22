from __future__ import print_function,division
import os
import sys
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from ete3 import Tree, NodeStyle, faces
from tqdm import tqdm
from copy import deepcopy
import itertools

def get_consensus(seqs):
    consensus = ''
    if type(seqs) is dict:
        seqs = list(seqs.values())
    for i in range(len(seqs[0])):
        letters = [j[i] for j in seqs]
        consensus+= max(set(letters), key=letters.count)
    return consensus


def translate_seq(seq,incomplete=False):
    codontable = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}

    IUPAC = {'A':'A','T':'T','G':'G','C':'C',
             'R':'AG','Y':'CT','M':'AC','K':'GT','S':'GC','W':'AT',
             'H':'ATC','B':'GTC','V':'GAC','D':'GAT','N':'ATGC'}

    if incomplete:
        seq = seq[:len(seq)-(len(seq)%3)]
    codons = [seq[i:i+3] for i in range(0,len(seq),3)]
    try:
        return ''.join([codontable[i] for i in codons])
    except KeyError:
        result = ''
        for codon in codons:
            try:
                result+=codontable[codon]
            except KeyError:
                print(codon)
                codons = [[j for j in IUPAC[i]] for i in codon]
                pos_codons = [codontable[''.join(i)] for i in itertools.product(*codons)]
                if len(set(pos_codons))>1:
                    result+='X'
                else:
                    result+=pos_codons[0]
        return result

def difference_seqs(seq1,seq2):
    '''returns the positions where seq1 and seq2 differ'''
    if len(seq1) != len(seq2):
        raise ValueError("sequences must be of same length. length seq1:{}, seq2:{}".format(len(seq1),len(seq2)))
    IUPAC = {'A':['A'],
             'T':['T'],
             'G':['G'],
             'C':['C'],
             'R':['A','G'],
             'Y':['C','T'],
             'M':['A','C'],
             'K':['G','T'],
             'S':['G','C'],
             'H':['A','T','C'],
             'W':['A','T'],
             'V':['G','A','C'],
             'B':['G','T','C'],
             'D':['G','A','T'],
             'N':['A','T','G','C'],
             '-':['-']}


    differs = []
    for i in range(len(seq1)):
        if len(set(IUPAC[seq1[i]]).intersection(IUPAC[seq2[i]]))==0:
            differs.append(i)

    return differs



def read_fasta(filename,consensus=-1):
    '''
    reads in a fasta file as dictionary with sequence identifiers as keys and
    sequences as values
    '''
    seqs = {}
    counter = 0
    with open(filename) as f:
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

def seqs_from_fasta_string(string,consensus=-1):
    seqs = {}
    counter=0
    lines = string.split('\n')
    for line in lines:
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

def highlighter(seqs=None,filename=None,fasta_string=None,consensus=-1,show=True,outname='test.png'):
    '''
    creates a highlighter plot (as on hiv.lanl.gov) of already aligned sequences
    '''
    colors = {'A': '#bf6c60', 'T': '#99cc33', 'G': '#2a3326', 'C': '#36b8d9', '-': '#792080'}
    if seqs is None:
        if fasta_string is None:
            seqs = read_fasta(filename,consensus)
        else:
            seqs = seqs_from_fasta_string(fasta_string)
    to_plot = list(seqs.keys())
    seq_len = len(seqs[list(seqs.keys())[0]])
    current = 'consensus'
    pos = 0
    mins=None
    while len(to_plot) > 0:
        plt.plot([0,seq_len],[pos,pos],color='grey')
        plt.text(seq_len,pos,current,va='center',ha='left')
        changes = difference_seqs(seqs['consensus'], seqs[current])

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
    if show:
        plt.show()
    else:
        plt.savefig(outname,dpi=100 )
        plt.close()

def njTree(seqs=None,filename=None,consensus=-1,fasta_string=None,show=True):
    '''creates and prints/plots a simple nj tree'''
    if seqs is None:
        if fasta_string is None:
            seqs = read_fasta(filename,consensus)
        else:
            seqs = seqs_from_fasta_string(fasta_string,consensus=consensus)
    else:
        seqs = deepcopy(seqs)
    seq_len = len(list(seqs.values())[0])
    try:
        del seqs['consensus']
    except KeyError:
        pass

    diff = np.ones((len(seqs),len(seqs)))*(seq_len+1)
    for i in range(len(seqs)):
        for j in range(i+1,len(seqs)):
            diff[i,j] = sum(list(seqs.values())[i][k] != list(seqs.values())[j][k] for k in range(seq_len) )

    #print diff
    newick = ''
    names = list(seqs.keys())

    while True:
        min_value = np.min(diff)
        if min_value == seq_len+1:
            break
        min_location = np.where(np.min(diff)==diff)
        min_location = zip(min_location[0],min_location[1])

        counter=0
        for i in min_location:
            if names[i[0]]!='' and names[i[1]]!='':
                counter+=1
                names[i[0]] = '{}:{},{}'.format(names[i[0]],min_value, names[i[1]],min_value)
                #print names[i[0]]
                names[i[1]] = ''
                diff[i[1]:,i[0]] = (diff[i[1]:,i[1]]+diff[i[1]:,i[0]])/2
                diff[:,i[1]] = seq_len+1

                diff[i[0],i[1]:] = (diff[i[1],i[1]:]+diff[i[0],i[1]:])/2
                diff[i[1],:] = seq_len+1

        names = [add_brackets(i,':{}'.format(min_value)) for i in names]
        # print names
        # raw_input()
    newick = '('+[i for i in names if i != ''][0]+');'
    #print newick
    t = Tree(newick)
    #print t
    if show:
        t.show()
    return newick



def highlighter_one(seq,consensus,name='test.png',show=False):
    plt.figure(figsize=[8,1])
    colors = {'A': '#bf6c60', 'T': '#99cc33', 'G': '#2a3326', 'C': '#36b8d9', '-': '#792080'}
    plt.plot([0,len(seq)],[1,1],color='black',zorder=-1)
    changed = [i for i in range(len(seq)) if seq[i]!=consensus[i]]
    for i in changed:
        plt.scatter(i,1,color=colors[seq[i]],s=100)
    plt.axis('off')
    plt.tight_layout()
    if not show:
        plt.savefig(name,dpi=60 )
        plt.close()
    else:
        plt.show()


def advancedTree(seqs, colored=['Pop1', 'Pop2']):
    nw = njTree(seqs,show=False)
    colors = {name.replace(' ', ''):'#'+''.join(['{:02x}'.format(int(c*256)) for c in list(cm.jet(i))[:3]]) for i,name in zip(np.linspace(0,1,len(colored)),colored)}
    tree = Tree(nw)
    files = []
    for node in tree.traverse():
        nst = NodeStyle()
        nst['size'] = 0
        for col in colored:
            if col in node.name:
                nst['size'] = 10
                nst['fgcolor'] = colors[col]
            node.set_style(nst)
        if node.is_leaf():
            highlighter_one(seqs[node.name],seqs['consensus'],node.name+'.png')
            files.append(node.name+'.png')
            node.add_face(faces.ImgFace(node.name+'.png'),0,position='aligned')
    tree.show()
    for f in files:
        os.remove(f)

def diversity_index(seqs):
    '''returns the average number of differences between all the sequences'''
    seqs = [seqs[i] for i in seqs if i!='consensus']

    divers = []
    for i in tqdm(range(len(seqs) )):
        for j in range(i+1,len(seqs)):
            difs = difference_seqs(seqs[i], seqs[j])
            divers.append(len(difs))
    return sum(divers)/len(divers)

def add_brackets(txt,extra=''):
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
                return '({}{})'.format(txt,extra)
        else:
            return '({}{})'.format(txt,extra)
    else:
        return txt

if __name__ == '__main__':
    fasta_file = sys.argv[1]
    try:
        action = sys.argv[2]
    except IndexError:
        print('no action given, assuming highlighter')
        action = 'highlighter'


    seqs = read_fasta(fasta_file,consensus=-1)

    if action == 'highlighter':
        fig = plt.figure(figsize=[10,20])
        highlighter(seqs,consensus = 0,show=True,outname=fasta_file.split('.')[0]+'.pdf')
    if action == 'tree':
        njTree(seqs)
    if action == 'advancedTree':
        advancedTree(seqs,colored=['pop1', 'pop2'])
    if action == 'diversityIndex':
        print(diversity_index(seqs))
