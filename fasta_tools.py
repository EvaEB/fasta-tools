import sys

def read_fasta(filename):
    seqs = {}
    with open(fasta_file) as f:
        for line in f.readlines():
            if '>' in line:
                current = line.strip()[1:]
                seqs[current] = ''
            else:
                seqs[current]+=line.strip()

    return seqs


if __name__ == '__main__':
    fasta_file = sys.argv[1]
    seqs = read_fasta(fasta_file)
    print seqs
