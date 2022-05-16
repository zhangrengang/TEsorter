import sys
from collections import OrderedDict
from .RunCmdsMP import run_cmd, logger

BLASType = {
    'qseqid': str,
    'sseqid': str,
    'pident': float,
    'length': int,
    'mismatch': int,
    'gapopen': int,
    'qstart': int,
    'qend': int,
    'sstart': int,
    'send': int,
    'evalue': float,
    'bitscore': float,
    'qlen': int,
    'slen': int,
    'qcovs': float,
    'qcovhsp': float,
    'sstrand': str,
	}

def blast(db_seq, qry_seq, seqtype='nucl', db_name=None, blast_out=None, blast_outfmt='6', blast_opts='', ncpu=4):
	if seqtype == 'nucl':
		blast_app = 'blastn'
	elif seqtype == 'prot':
		blast_app = 'blastp'
	else:
		raise ValueError('Unknown molecule type "{}" for blast'.format(seqtype))

	if db_name is None:
		db_name = db_seq
	if blast_out is None:
		blast_out = qry_seq + '.blastout'

	cmd = 'makeblastdb -in {} -dbtype {} -out {}'.format(db_seq, seqtype, db_name)
	run_cmd(cmd, logger=logger, fail_exit=True)

	cmd = '{} -query {} -db {} -out {} -outfmt {} -num_threads {} {}'.format(
			blast_app, qry_seq, db_name, blast_out, blast_outfmt, ncpu, blast_opts)
#	cmd += " " + blast_opts
	run_cmd(cmd, logger=logger)
	return blast_out

class BlastOut(object):
	def __init__(self, blast_out, outfmt=None):
		self.blast_out = blast_out
		if outfmt is not None:
			outfmt = outfmt.strip(''''"''')
		if outfmt is None:
			self.outfmt = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'.split()
		elif outfmt[0] == '6':
			self.outfmt = outfmt.split()[1:]
		else:
			raise ValueError('Only support for blast outfmt 6 = tabular, but {} input'.format(outfmt))
	def __iter__(self):
		return self.parse()
	def parse(self):
		for line in open(self.blast_out):
			values = line.strip().split('\t')
			yield BlastOutRecord(self.outfmt, values)
	def filter_besthit(self, fout=sys.stdout):
		d_best_hit = OrderedDict()
		for rc in self.parse():
			if rc.qseqid in d_best_hit:
				if rc.bitscore > d_best_hit[rc.qseqid].bitscore:
					d_best_hit[rc.qseqid] = rc
			else:
				d_best_hit[rc.qseqid] = rc
		if fout is None:
			return d_best_hit
		for qseqid, rc in d_best_hit.items():
			rc.write(fout)
		return d_best_hit

class BlastOutRecord(object):
	def __init__(self, outfmt, values):
		self.values = values
		for key, value in zip(outfmt, values):
			setattr(self, key, BLASType[key](value))
	def write(self, fout=sys.stdout):
		print('\t'.join(self.values), file=fout)
	@property
	def scov(self):
		return 1.0*(abs(self.send - self.sstart) +1)/self.slen
