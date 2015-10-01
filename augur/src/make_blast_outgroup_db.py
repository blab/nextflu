from Bio import SeqIO
import os

fname =  'source-data/vaccrefMix_HAanno_nuc_Sep2015.fa'
ofname = 'source-data/outgroups.fasta'

with open(ofname, 'w') as ofile:
    for seq in SeqIO.parse(fname, 'fasta'):
        seq.seq = seq.seq.translate()
        seq.id = seq.name = '|'.join(seq.description.split()[1].split('|')[:2]).replace(' ','')
        seq.description=''
        SeqIO.write(seq, ofile, 'fasta')

os.system('makeblastdb -in '+ofname+ ' -dbtype prot ')
