from Bio import SeqIO

infname = 'data/ig_ebo_145_29Jun15_correct.fas'
informat = {'lab':0, 'strain':1, 'country':2, 'region':3, 'date':4}
outformat = ['virus', 'accession', 'strain', 'lab', 'country', 'region','location','smth', 'date']
outfname = 'data/sanger_27Jun15.fas'
existing_file_name = 'data/ebola_908.fasta'

generic = {'virus':'EBOV', 'accession':'', 'lab':'', 'country':'', 'region':'', 'location':'', 'smth':'' }

existing = []
for seq in SeqIO.parse(existing_file_name, 'fasta'):
	existing.append(seq.name.split('|')[outformat.index('strain')].upper())

skipped=0
written=0
with open(outfname, 'w') as outfile:
	for seq in SeqIO.parse(infname, 'fasta'):
		fields = seq.name.rstrip('*').split('|')
		params = {key:fields[informat[key]] for key in informat}
		new_name = '|'.join([params[key] if key in params else generic[key]
			for key in outformat])

		seq.name=new_name
		seq.id=new_name
		seq.description=''
		if fields[informat['strain']].upper() in existing:
			print 'skipped:',seq.name
			skipped+=1
		else:
			print 'written:',seq.name			
			SeqIO.write(seq, outfile, 'fasta')
			written+=1

print "written",written,"new sequences"
print "skipped",skipped,"sequences"