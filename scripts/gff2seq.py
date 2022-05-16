import sys
import re
import time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from xopen import xopen as open

in_ref_fa = sys.argv[1]
in_gtf = sys.argv[2]

out_fa = sys.argv[3]
try: type = sys.argv[4]
except: type = 'CDS'

d_ref_fa = {}
for record in SeqIO.parse(open(in_ref_fa),'fasta'):
	d_ref_fa[record.id] = record.seq # Seq

d_transcript = {}
d_strand = {}
d_gene_id = {}
d_chr = {}
i = 0
for line in open(in_gtf,'r'):
	i += 1
	if not line.strip(): 
		continue
	if re.compile(r'#').match(line):
		continue
	temp = line.strip().split('\t')
	if len(temp) < 9:
		continue
	if temp[2] != type:
		continue
	chr = temp[0]
	start = int(temp[3])-1
	end = int(temp[4])
	try: gene_id = re.compile(r'Parent=(\S+?)($|;)').search(temp[8]).groups()[0]
	except AttributeError: gene_id = re.compile(r'ID=(\S+?)($|;)').search(temp[8]).groups()[0]
	d_chr[gene_id] = chr
	strand = temp[6]

	if gene_id in d_strand:
		if d_strand[gene_id] != strand:
			print('Error: Id %s has conflict strand' % gene_id)
		else:
			pass
	else:
		d_strand[gene_id] = strand

	if  chr in d_ref_fa:
		try:
			d_transcript[gene_id].append((start,end))
		except:
			d_transcript[gene_id] = [(start,end)]
	else:
		print('Error: Seq %s is not in genomes' % chr)

def merge_seqs(seqs):
	for seq in seqs:
		try:
			myseq += seq
		except:
			myseq = seq
	return myseq

f = open(out_fa,'w')
for key,values in sorted(d_transcript.items()):
	myseqs = []
	for value in sorted(values):
		myseqs.append(d_ref_fa[d_chr[key]][value[0]:value[1]])
	myseq = merge_seqs(myseqs)
	if d_strand[key] == '+':
		myseq = str(myseq)
	elif d_strand[key] == '-':
		myseq = str(myseq.reverse_complement())
	else:
		print('Error: strand is not either "+" or "-".')
	print('>%s\n%s' % (key,myseq), file=f)

f.close()


