from Bio import SeqIO

seq_fname = 'data/MERS-CoV.fasta'
anno_fname = 'data/MERS-CoV.csv'

annos = {}
with open(seq_fname, 'r') as infile:
    for seq in SeqIO.parse(infile, 'fasta'):
        annos[seq.id] = seq.description.split('|')

with open(anno_fname, 'w') as outfile:
    for fields in annos.values():
        outfile.write(','.join(fields)+'\n')
