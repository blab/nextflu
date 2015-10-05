from Bio import SeqIO
import os

fname =  'source-data/outgroups_nucleotides.fasta'
ofname = 'source-data/outgroups.fasta'

ref_viruses = []
with open(ofname, 'w') as ofile:
    for seq in SeqIO.parse(fname, 'fasta'):
        seq.seq = seq.seq.translate()
        ref_viruses.append(seq.description.split()[1].split('|')[1].replace(' ',''))
        seq.id = seq.name = '|'.join(seq.description.split()[1].split('|')[:2]).replace(' ','')
        seq.description=''
        SeqIO.write(seq, ofile, 'fasta')

os.system('makeblastdb -in '+ofname+ ' -dbtype prot ')

with open('ref_strain_names.txt', 'w') as ofile:
    for strain in set(ref_viruses):
        ofile.write(strain+'\n')
