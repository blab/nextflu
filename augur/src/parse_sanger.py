from Bio import SeqIO

infname = 'data/ig_ebo_145_27Jun15.fas'
informat = {'lab':0, 'strain':1, 'country':2, 'region':3, 'date':4}
outformat = ['virus', 'accession', 'strain', 'lab', 'country', 'region','location','smth', 'date']
outfname = 'data/sanger_27Jun15.fas'

generic = {'virus':'EBOV', 'accession':'', 'lab':'', 'country':'', 'region':'', 'location':'', 'smth':'' }

with open(outfname, 'w') as outfile:
	for seq in SeqIO.parse(infname, 'fasta'):
		print seq.name
		fields = seq.name.rstrip('*').split('|')
		params = {key:fields[informat[key]] for key in informat}
		new_name = '|'.join([params[key] if key in params else generic[key]
			for key in outformat])

		seq.name=new_name
		seq.id=new_name
		seq.description=''
		SeqIO.write(seq, outfile, 'fasta')
