from Bio import Entrez, SeqIO
from StringIO import StringIO
from datetime import datetime

Entrez.email = "richard.neher@tuebingen.mpg.de"     # Always tell NCBI who you are
handle = Entrez.esearch(db="nucleotide", term='"Middle East respiratory syndrome coronavirus"[Organism] OR MERS-CoV[All Fields]', retmax=1000)
record = Entrez.read(handle)

translation = {
    "Homo sapiens":"human",
    "Saudi Arabia":"KSA",
    "United Arab Emirates":"UEA",
    "Oman":"OMN",
    "Jordan":"JOR",
}

def fix_country(x):
    if x in translation:
        return translation[x]
    else:
        return x

def fix_host(x):
    if "homo" in x.lower():
        return "Human"
    elif "camel" in x.lower():
        return "Camel"
    else:
        return x

def get_date(x):
    try:
        return datetime.strptime(x.upper(),'%d-%b-%Y').isoformat(sep=' ').split()[0]
    except:
        try:
            return datetime.strptime(x.upper(),'%b-%Y').isoformat(sep=' ').split()[0]
            try:
                return datetime.strptime(x.upper(),'%Y').isoformat(sep=' ').split()[0]
            except:
                return x
        except:
            return x


with open('data/genbank_mers_seqs.fasta', 'w') as mersout:
    for genbank_id in ['KT006149']: #record["IdList"]+['KR011263']+['KJ477102']:
        handle = Entrez.efetch(db="nucleotide", id=genbank_id, rettype="gb")
        seq= SeqIO.read(StringIO(handle.read()), format = 'genbank')
        if len(seq)>5000:
            src= seq.features[0].qualifiers
            fasta_fields = []

            if 'strain' in src:
                fasta_fields.append(src['strain'][0])
            elif 'isolate' in src:
                fasta_fields.append(src['isolate'][0])
            else:
                fasta_fields.append('')

            fasta_fields.append(genbank_id)

            if 'collection_date' in src:
                fasta_fields.append(get_date(src['collection_date'][0]))
            else:
                fasta_fields.append('')

            if 'host' in src:
                fasta_fields.append(fix_host(src['host'][0]))
            else:
                fasta_fields.append('')

            if 'country' in src:
                fasta_fields.append(fix_country(src['country'][0]))
            else:
                fasta_fields.append('')


            if 'location' in src:
                fasta_fields.append(src['location'][0])
            else:
                fasta_fields.append('')

            if 'lab' in src:
                fasta_fields.append(src['lab'][0])
            else:
                fasta_fields.append('')

            fasta_header = '|'.join(fasta_fields)  
            print fasta_header
            seq.id = fasta_header
            SeqIO.write(seq, mersout, format='fasta')

