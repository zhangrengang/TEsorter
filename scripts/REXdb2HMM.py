'''python REXdb2HMM.py protein_database_viridiplantae_v3.0.fasta'''
import sys,os
from Bio import SeqIO
from RunCmdsMP import run_cmd
class REXdb:
	def __init__(self, dbfas, tmpdir='/tmp'):
		self.dbfas = dbfas
		self.tmpdir = tmpdir
	def BinSeqsByClade(self):
		d_seq_bins = {}
		for rc in SeqIO.parse(self.dbfas, 'fasta'):
			clade = rc.description.split('#')[1]
			try: d_seq_bins[clade] += [rc]
			except KeyError: d_seq_bins[clade] = [rc]
		return d_seq_bins
	def pipeline(self):
		for clade, records in self.BinSeqsByClade().iteritems():
			clade0 = clade.replace('/', '__').replace(':', '--')
			cladeSeqs = '{}/{}.fa'.format(self.tmpdir, clade0)
			with open(cladeSeqs, 'w') as f:
				for rc in records:
					SeqIO.write(rc, f, 'fasta')
			alnSeqs = cladeSeqs + '.aln'
			cmd = 'mafft --auto {} > {} 2> /dev/null'.format(cladeSeqs, alnSeqs)
#			cmd = 'mafft --auto {} 2> /dev/null| prepareAlign | mafft --auto - > {} 2> /dev/null'.format(
				cladeSeqs, alnSeqs)
			run_cmd(cmd)
			alnHMM = cladeSeqs + '.hmm'
			cmd = 'hmmbuild -n {} {} {} > /dev/null'.format(clade, alnHMM, alnSeqs)
			run_cmd(cmd)
			
def main(dbfas=sys.argv[1]):
	REXdb(dbfas).pipeline()

if __name__ == '__main__':
	main()
				
