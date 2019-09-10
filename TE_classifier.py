#!/bin/env python
# coding: utf-8
'''# Author: zrg1989@qq.com
'''
import sys
import os
import re
import glob
import argparse
from math import log
from collections import Counter, OrderedDict
from Bio import SeqIO

bindir = os.path.dirname(os.path.realpath(__file__))
sys.path = [bindir + '/bin'] + sys.path
from RunCmdsMP import pp_run
from split_records import split_fastx_by_chunk_num

__version__ = '0.1'

DB = {
	'gydb' : bindir + '/database/GyDB2.hmm',
	'rexdb': bindir + '/database/REXdb_protein_database_viridiplantae_v3.0_plus_metazoa_v3.hmm',
	'rexdb-plant': bindir + '/database/REXdb_protein_database_viridiplantae_v3.0.hmm',
	'rexdb-metazoa': bindir + '/database/REXdb_protein_database_metazoa_v3.hmm',
	}
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
	
def pipeline(args):
	if not os.path.exists(args.tmp_dir):
		os.makedirs(args.tmp_dir)
		
	# search against DB and parse
	gff, geneSeq = LTRlibAnn(
			ltrlib = args.sequence, 
			hmmdb = args.hmm_database, 
			seqtype = args.seq_type,
			prefix = args.prefix,
			processors = args.processors,
			tmpdir = args.tmp_dir,
			mincov = args.min_coverage,
			maxeval = args.max_evalue,
			)
			
	# classify	
	classify_out = args.prefix + '.cls'
	fc = open(classify_out, 'w')
	d_class = OrderedDict()
	for rc in Classifier(gff, db=args.hmm_database, fout=fc):
		d_class[rc.ltrid] = rc
	fc.close()
	if not args.disable_pass2:
		
	# pass-2 classify
class BlastClassifier(object):
	def __int__(self, ltrlib, d_class, seqtype='nucl', ncpu=4):
		
	
def classify_by_blast(db_seq, qry_seq, blast_out=None, seqtype='nucl', ncpu=4, 
					  min_identtity=80, min_coverge=80, min_length=80):
	blast_outfmt = '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs qcovhsp sstrand'
	blast_out = blast(db_seq, qry_seq, seqtype=seqtype, blast_out=blast_out, blast_outfmt=blast_outfmt, ncpu=ncpu)
	d_best_hit = BlastOut(blast_out, blast_outfmt).filter_besthit(fout=None)
	for qseqid, rc in d_best_hit.iteritems():
		if not (rc.pindent >= min_identtity and rc.qcovs >= min_coverge and rc.qlen >= min_length):
			del d_best_hit[qseqid]
	d_class = OrderedDict([(qseqid, rc.sseqid) for qseqid, rc in d_best_hit.iteritems()])
	return d_class
	
def blast(db_seq, qry_seq, seqtype='nucl', blast_out=None, blast_outfmt=None, ncpu=4):
	if seqtype == 'nucl':
		blast = 'blastn'
	elif seqtype == 'prot':
		blast = 'blastp'
	else:
		raise ValueError('Unknown molecule type "{}" for blast'.format(seqtype))
	if blast_out is None:
		blast_out = qry_seq + '.blastout'
	if blast_outfmt is None:
		blast_outfmt = '6'
	cmd = 'makeblastdb -in {} -dbtype {}'.format(db_seq, seqtype)
	os.system(cmd)
		
	cmd = '{} -query {} -db {} -out {} -outfmt {} -num_threads {}'.format(qry_seq, db_seq, blast_out, blast_outfmt, ncpu)
	os.system(cmd)
	return blast_out

class BlastOut(object):
	def __init__(self, blast_out, outfmt=None):
		self.blast_out = blast_out
		if outfmt is None:
			self.outfmt = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'.split()
		elif outfmt[0] == '6':
			self.outfmt = outfmt.split()[1:]
		else:
			raise ValueError('Only support for blast outfmt 6 = tabular')
	def __iter__(self):
		return self.parse()
	def parse(self):
		for line in open(self.blast_out):
			values = line.strip().split('\t')
			yield BlastOutRecord(self.outfmt, values)
	def filter_besthit(self, fout=None):
		d_best_hit = OrderedDict()
		for rc in self.parse():
			if rc.qseqid in d_best_hit:
				if rc.bitscore > d_best_hit[rc.qseqid]:
					d_best_hit[rc.qseqid] = rc
			else:
				d_best_hit[rc.qseqid] = rc
		if fout is None:
			return d_best_hit
		for qseqid, rc in d_best_hit.iteritems():
			rc.write(fout)
class BlastOutRecord(object):
	def __init__(self, outfmt, values):
		self.values = values
		for key, value in zip(outfmt, values):
			setattr(self, key, BLASType[key](value))
	def write(self, fout=sys.stdout):
		print >> fout, '\t'.join(self.values)
class Classifier(object):
	def __init__(self, gff, db='rexdb', fout=sys.stdout): # gff is sorted
		self.gff = gff
		self.db = db
		self.fout = fout
		if self.db.startswith('rexdb'):
			self.markers = {'GAG', 'PROT', 'INT', 'RT', 'RH'}
		elif self.db == 'gydb':
			self.markers = {'GAG', 'AP', 'INT', 'RT', 'RNaseH', 'ENV'}
	def __iter__(self):
		return self.classify()
	def parse(self):
		record = []
		last_lid = ''
		for line in open(self.gff):
			if line.startswith('#'):
				continue
			line = LTRgffLine(line)
			lid = line.ltrid
			if record and  not lid == last_lid:
				yield record
				record = []
			record.append(line)
			last_lid = lid
		yield record
	def classify(self, ):
		line = ['#TE', 'Order', 'Superfamily', 'Clade', 'Code', 'Strand', 'hmmmatchs']
		print >> self.fout, '\t'.join(line)
		for rc in self.parse():
			rc_flt = rc #[line for line in rc if line.gene in self.markers]
			#if len(rc_flt) == 0:
			#	rc_flt = [line for line in rc]
			strands = [line.strand for line in rc_flt]
			if len(set(strands)) >1 :
				strand = '?'
			else:
				strand = strands[0]
			if strand == '-':
				rc_flt.reverse()
				rc.reverse()
			lid = rc_flt[0].ltrid
			domains = ' '.join(['{}|{}'.format(line.gene, line.clade)  for line in rc])
			genes  = [line.gene  for line in rc_flt]
			clades = [line.clade for line in rc_flt]
			names = [line.name for line in rc_flt]
			if self.db.startswith('rexdb'):
				order, superfamily, max_clade, coding = self.identify_rexdb(genes, names)
			elif self.db == 'gydb':
				order, superfamily, max_clade, coding = self.identify(genes, clades)
			line = [lid, order, superfamily, max_clade, coding, strand, domains]
			print >> self.fout, '\t'.join(line)
			self.ltrid, self.order, self.superfamily, self.clade, self.code, self.strand, self.domains = line
			yield self
	def identify_rexdb(self, genes, clades):
		perfect_structure = {
            ('LTR', 'Copia'): ['Ty1-GAG', 'Ty1-PROT', 'Ty1-INT', 'Ty1-RT', 'Ty1-RH'],
            ('LTR', 'Gypsy'): ['Ty3-GAG', 'Ty3-PROT', 'Ty3-RT', 'Ty3-RH', 'Ty3-INT'],
			}
		clade_count = Counter(clades)
		max_clade = max(clade_count, key=lambda x: clade_count[x])
		order, superfamily = self._parse_rexdb(max_clade)
		if len(clade_count) == 1:
			max_clade = max_clade.split('/')[-1]
		elif len(clade_count) > 1:
			max_clade = 'mixture'
			superfamlies = [self._parse_rexdb(clade)[1] for clade in clades]
			if len(Counter(superfamlies)) > 1:
				superfamily = 'mixture'
				orders = [self._parse_rexdb(clade)[0] for clade in clades]
				if len(Counter(orders)) > 1:
					order = 'mixture'
		try:
			ordered_genes = perfect_structure[(order, superfamily)]
			my_genes = [gene for gene in genes if gene in set(ordered_genes)]
			if ordered_genes == my_genes:
				coding = 'cmpl' # completed gene structure
			else:
				coding = 'lost'
		except KeyError:
			coding = 'unknown'
		return order, superfamily, max_clade, coding
	def _parse_rexdb(self, clade): # full clade name
		if clade.startswith('Class_I/LTR/Ty1_copia'):
			order, superfamily = 'LTR', 'Copia'
		elif clade.startswith('Class_I/LTR/Ty3_gypsy'):
			order, superfamily = 'LTR', 'Gypsy'
		elif clade.startswith('Class_I/'): # LINE, pararetrovirus, Penelope, DIRS
			order, superfamily = clade.split('/')[1], 'unknown'
		elif clade.startswith('Class_II/'):
			try: order, superfamily = clade.split('/')[2:4]
			except ValueError: order, superfamily = clade.split('/')[2], 'unknown'
		return order, superfamily
	def identify(self, genes, clades):
		perfect_structure = {
			('LTR', 'Copia')         : ['GAG', 'AP', 'INT', 'RT', 'RNaseH'],
			('LTR', 'Gypsy')         : ['GAG', 'AP', 'RT', 'RNaseH', 'INT'],
			('LTR', 'Pao')           : ['GAG', 'AP', 'RT', 'RNaseH', 'INT'],
			('LTR', 'Retroviridae')  : ['GAG', 'AP', 'RT', 'RNaseH', 'INT', 'ENV'],
			('LTR', 'Caulimoviridae'): ['GAG', 'AP', 'RT', 'RNaseH', ],
			}
		d_map = self.clade_map
		clade_count = Counter(clades)
		max_clade = max(clade_count, key=lambda x: clade_count[x])
		try: (order, superfamily) = d_map[max_clade]
		except KeyError: 
			(order, superfamily) = ('Unknown', 'unknown')
			print >>sys.stderr, 'unknown clade: {}'.format(max_clade)
		try:
			ordered_genes = perfect_structure[(order, superfamily)]
			my_genes = [gene for gene in genes if gene in set(ordered_genes)]
			if ordered_genes == my_genes:
				coding = 'cmpl' # completed gene structure and the same order
			else:
				coding = 'lost'
		except KeyError:
			coding = 'unknown'
		return order, superfamily, max_clade, coding
	def replace_annotation(self, rawseq, fout=sys.stdout, idmap = None):
		d_class = self.classification
		i = 0
		for rc in SeqIO.parse(rawseq, 'fasta'):
			if idmap is None:
				intid = rc.id.split('#')[0]
			else:
				try: intid = idmap[rc.id.split('#')[0]]
				except KeyError as e:	# this should be rare
					print >>sys.stderr, '[Warning] skipped', e
			if intid in d_class:
				neword, newfam = d_class[intid]
				re_org = self.re_orgnize(rc.id, neword, newfam)
				if re_org:
					i += 1
					rc.id = re_org
#					print >>sys.stderr, rc.description, len(rc.seq)
			SeqIO.write(rc, fout, 'fasta')
		print >> sys.stderr, i, 'sequences re-classified'
	def re_orgnize(self, rawid, neword, newfam):
		rawid, rawcls = rawid.split('#')
		try: raword, rawfam = rawcls.split('/')[:2]
		except ValueError: raword, rawfam = rawcls, 'unknown'
		if raword.lower() == 'unknown' or (rawfam.lower() == 'unknown' and raword == neword):
			return '{}#{}/{}'.format(rawid, neword, newfam)
		else:
			return False
	@property
	def classification(self):
		return {rc.ltrid.split('#')[0]: (rc.order, rc.superfamily) for rc in self.classify()}
	@property
	def clade_map(self):
		return {rc.clade: (rc.order, rc.superfamily) for rc in CladeInfo()}
			
class CladeInfo():
	def __init__(self, infile=DB['gydb']+'.info'):
		self.infile = infile
	def __iter__(self):
		return self.parse()
	def parse(self):
		i = 0
		for line in open(self.infile):
			i += 1
			temp = line.strip().split('\t')
			if i == 1:
				title = temp
				continue
			self.dict = dict(zip(title, temp))
			if self.dict['Clade'] == 'NA':
				self.clade = self.dict['Cluster_or_genus']
			else:
				self.clade = self.dict['Clade']
			if self.clade == '17.6': # exception
				self.clade = '17_6'
			self.superfamily = self.dict['Family'].split('/')[-1]
			if self.superfamily == 'Retroviridae':	# deltaretroviridae gammaretroviridae
				self.clade = self.dict['Cluster_or_genus'].replace('virus', 'viridae')
			if self.superfamily == 'Retrovirus':	# an exception
				self.superfamily = 'Retroviridae'
			self.order = 'LTR' if self.dict['System'] in {'LTR_retroelements', 'LTR_Retroelements', 'LTR_retroid_elements'} else self.dict['System']
			self.clade = self.clade.replace('-', '_') # A-clade V-clade C-clade
			
			yield self
			self.clade = self.clade.lower()
			yield self
		self.order, self.superfamily, self.clade, self.dict = ['LTR', 'Copia', 'ty1/copia', {}]  # AP_ty1/copia
		yield self
		for clade, order in zip(['retroelement', 'shadow', 'all'], ['LTR', 'Unknown', 'Unknown']): # CHR
			self.order, self.superfamily, self.clade, self.dict = [order, 'unknown', clade, {}]  # CHR_retroelement
			yield self	
				
class GffLine(object):
	def __init__(self, line):
		temp = line.strip().split('\t')
		self.chr, self.source, self.type, self.start, self.end, self.score, self.strand, self.frame, self.attributes = temp
		self.start, self.end = int(self.start), int(self.end)
		try: self.score = float(self.score)
		except: pass
		try: self.frame = int(self.frame)
		except: pass
		self.attributes = self.parse(self.attributes)
	def parse(self, attributes):
		return dict(kv.split('=') for kv in attributes.split(';'))
class LTRgffLine(GffLine):
	def __init__(self, line):
		super(LTRgffLine, self).__init__(line)
		self.gene = self.attributes['gene']
		self.clade = self.attributes['clade']
		self.ltrid = '|'.join(self.attributes['ID'].split('|')[:-1])
		self.name = self.attributes['ID'].split('|')[-1].split(':')[0]

class HmmScan():
	def __init__(self, hmmout, hmmfmt='domtbl'):
		self.hmmout = hmmout
		self.hmmfmt = hmmfmt
	def __iter__(self):
		return self.parse()
	def parse(self):
		for line in open(self.hmmout):
			if line.startswith('#'):
				continue
			if self.hmmfmt == 'domtbl':
				yield HmmDomRecord(line)
class HmmDomRecord():
	def __init__(self, line):
		temp = line.strip().split()
		self.tname, self.tacc, self.tlen, self.qname, self.qacc, self.qlen, \
			self.evalue, self.score, self.bias, \
			self.domi, self.domn, self.cevalue, self.ievalue, self.domscore, self.dombias, \
			self.hmmstart, self.hmmend, self.alnstart, self.alnend, self.envstart, self.envend, \
			self.acc \
				= temp[:22]
		self.tlen, self.qlen, self.domi, self.domn, \
			self.hmmstart, self.hmmend, self.alnstart, self.alnend, self.envstart, self.envend = \
			map(int, [self.tlen, self.qlen, self.domi, self.domn, \
				self.hmmstart, self.hmmend, self.alnstart, self.alnend, self.envstart, self.envend])
		self.evalue, self.score, self.bias, self.cevalue, self.ievalue, self.domscore, self.dombias, self.acc = \
			map(float, [self.evalue, self.score, self.bias, self.cevalue, self.ievalue, self.domscore, self.dombias, self.acc])
		self.tdesc = ' '.join(temp[22:])
	@property
	def hmmcov(self):
		return round(1e2*(self.hmmend - self.hmmstart + 1) / self.tlen, 1)
def seq2dict(inSeq):
	return dict([(rc.id, rc) for rc in SeqIO.parse(inSeq, 'fasta')])
	
def parse_hmmname(hmmname, db='gydb'):
	db = db.lower()
	if db == 'gydb':
		temp = hmmname.split('_')
		gene, clade = temp[0], '_'.join(temp[1:])
	elif db.startswith('rexdb'):	# Class_I/LTR/Ty3_gypsy/chromovirus/Tekay:Ty3-RT
		gene = hmmname.split(':')[1] #.split('-')[1]
		clade = hmmname.split(':')[0].split('/')[-1]
	elif db.startswith('pfam'):
		gene = hmmname
		clade = hmmname

	return gene, clade

def hmm2best(inSeq, inHmmouts, prefix=None, db='rexdb', seqtype='nucl', mincov=20, maxeval=1e-3):
	if prefix is None:
		prefix = inSeq
	d_besthit = {}
	for inHmmout in inHmmouts:
		for rc in HmmScan(inHmmout):
			suffix = rc.qname.split('|')[-1]
			if suffix.startswith('aa') or suffix.startswith('rev_aa'):
				qid = '|'.join(rc.qname.split('|')[:-1])
			else:
				qid = rc.qname
			domain,clade = parse_hmmname(rc.tname, db=db)
			if db.startswith('rexdb'):
				cdomain = domain.split('-')[1]
				if cdomain == 'aRH':
					cdomain = 'RH'
				key = (qid, cdomain)
				if key in d_besthit:
					best_rc = d_besthit[key]
					if rc.score > best_rc.score:
						best_domain, _ = parse_hmmname(best_rc.tname, db=db)
						if domain == best_domain:
							d_besthit[key] = rc
						elif rc.envstart <= best_rc.envend and rc.envend >= best_rc.envstart: # overlap
							d_besthit[key] = rc
				else:
					d_besthit[key] = rc
			else:
				key = (qid, domain)
				if key in d_besthit:
					if rc.score > d_besthit[key].score:
						d_besthit[key] = rc
				else:
					d_besthit[key] = rc
#	print d_besthit
	d_seqs = seq2dict(inSeq)
	lines = []
	for (qid, domain), rc in d_besthit.items():
		if rc.hmmcov < mincov or rc.evalue > maxeval:
			continue
#		gid = '{}|{}|{}'.format(qid, domain, rc.tname)
		rawid = qid
#		clade = '_'.join(rc.tname.split('_')[1:])
		gene,clade = parse_hmmname(rc.tname, db=db)
		if db.startswith('rexdb'):
			domain = gene
		gid = '{}|{}'.format(qid, rc.tname)
		gseq = d_seqs[rc.qname].seq[rc.envstart-1:rc.envend]
		if seqtype == 'nucl':
			strand, frame = parse_frame(rc.qname.split('|')[-1])
			if strand == '+':
				nuc_start = rc.envstart * 3 - 2  + frame
				nuc_end = rc.envend* 3 + frame
			elif strand == '-':
				nuc_start = rc.qlen*3 - (rc.envend* 3 + frame) + 1
				nuc_end = rc.qlen*3 - (rc.envstart* 3 + frame) + 1
			else:
				nuc_start = rc.envstart
				nuc_end = rc.envend
		elif seqtype == 'prot':
			strand, frame = '+', '.'
			nuc_start, nuc_end, = rc.envstart, rc.envend
		match = re.compile(r'(\S+?):(\d+)\.\.(\d+)').match(qid)
		if match:
			qid, ltrstart, ltrend = match.groups()
			ltrstart = int(ltrstart)
			nuc_start = ltrstart + nuc_start - 1
			nuc_end = ltrstart + nuc_end -1
		attr = 'ID={};gene={};clade={};evalue={};coverage={};probability={}'.format(gid, domain, clade, rc.evalue, rc.hmmcov, rc.acc)
		gffline = [qid, 'ltrapl', 'CDS', nuc_start, nuc_end, rc.score, strand, frame, attr, rc.evalue, rc.hmmcov, rc.acc, rawid, gid, gseq]
		lines.append(gffline)
	gff, seq, tsv = '{}.gff3'.format(prefix), '{}.faa'.format(prefix), '{}.tsv'.format(prefix)
	fgff = open(gff, 'w')
	fseq = open(seq, 'w')
	ftsv = open(tsv, 'w')
	print >> ftsv, '\t'.join(['#id', 'length', 'evalue', 'coverge', 'probability'])
	for line in sorted(lines, key=lambda x: (x[0], x[-3], x[3])):
		gffline = line[:9]
		gffline = map(str, gffline)
		print >> fgff, '\t'.join(gffline)
		gid, gseq = line[-2:]
		gdesc = line[8]
		print >> fseq, '>{} {}\n{}'.format(gid, gdesc, gseq)
		evalue, hmmcov, acc = line[-6:-3]
		line = [gid, len(gseq), evalue, hmmcov, acc]
		print >> ftsv, '\t'.join(map(str, line))
	fgff.close()
	fseq.close()
	return gff, seq
def translate(inSeq, prefix=None):
	if prefix is None:
		prefix = inSeq
	prog = 'perl {}/bin/Six-frame_translate.pl'.format(bindir)
	outSeq = prefix + '.aa'
	cmd = '{} {} > {}'.format(prog, inSeq, outSeq)
	os.system(cmd)
	return outSeq
def hmmscan(inSeq, hmmdb='rexdb.hmm', hmmout=None):
	if hmmout is None:
		hmmout = prefix + '.domtbl'
	cmd = 'hmmscan --notextw -E 0.01 --domE 0.01 --noali --cpu 4 --domtblout {} {} {} > /dev/null'.format(hmmout, hmmdb, inSeq)
	os.system(cmd)
	return hmmout
def hmmscan_pp(inSeq, hmmdb='rexdb.hmm', hmmout=None, tmpdir='./tmp'):
	chunk_prefix = '{}/{}'.format(tmpdir, 'chunk_aaseq')
	_, _, _, chunk_files = split_fastx_by_chunk_num(
			inSeq, prefix=chunk_prefix, chunk_num=processors, seqfmt='fasta', suffix='')
	domtbl_files = [chunk_file + '.domtbl' for chunk_file in chunk_files]
	commands = [ 
		'hmmscan --notextw -E 0.01 --domE 0.01 --noali --domtblout {} {} {}'.format(domtbl_file, hmmdb, chunk_file) \
			for chunk_file, domtbl_file in zip(chunk_files, domtbl_files)]
	jobs = pp_run(commands)
	for cmd, (stdout, stderr, status) in zip(commands, jobs):
		if not status == 0:
			print >>sys.err, "Warning: exit code {} for CMD '{}'".format(status, CMD)
	# cat files
	if hmmout is None:
		hmmout = prefix + '.domtbl'
	with open(hmmout, 'w') as f:
		for domtbl_file in domtbl_files:
			for line in open(domtbl_file):
				f.write(line)
	return hmmout

def LTRlibAnn(ltrlib, hmmdb='rexdb', seqtype='dna', prefix=None, 
			force_write_hmmscan=False,
			processors=4, tmpdir='./tmp', 
			mincov=20, maxeval=1e-3):
	if prefix is None:
		prefix = '{}.{}'.format(ltrlib, hmmdb)
	
	if seqtype == 'nucl':
		print >>sys.stderr, 'translating {} in six frames'.format(ltrlib)
		aaSeq = translate(ltrlib)
	elif seqtype == 'prot':
		aaSeq = ltrlib
	
	print >>sys.stderr, 'HMM scanning against {}'.format(DB[hmmdb])
	domtbl = aaSeq + '.domtbl'
	if not (os.path.exists(domtbl) and os.path.getsize(domtbl) >0) or force_write_hmmscan:
		if processors > 1:
			hmmscan_pp(aaSeq, hmmdb=DB[hmmdb], hmmout=domtbl, tmpdir=tmpdir):
		else:
			hmmscan(aaSeq, hmmdb=DB[hmmdb], hmmout=domtbl)
	else:
		print >>sys.stderr, 'use existed {} and skip hmmscan'.format(domtbl)
	print >>sys.stderr, 'generating gene anntations'
	gff, geneSeq = hmm2best(aaSeq, [domtbl], db=hmmdb, prefix=prefix, seqtype=seqtype)
	return gff, geneSeq
def replaceCls(ltrlib, seqtype='dna', db='rexdb'):
	gff = ltrlib + '.' + db + '.gff3'
	if not os.path.exists(gff):
		gff, geneSeq, aaSeq = LTRlibAnn(ltrlib, seqtype=seqtype, hmmdb=db)
	annout = '{}.anno'.format(gff)
	newlib = '{}.reclassified'.format(ltrlib)
	fann = open(annout, 'w')
	flib = open(newlib, 'w')
	Classifier(gff, fout=fann).replace_annotation(ltrlib, fout=flib)
	fann.close()
	flib.close()
def parse_frame(string):
	if string.startswith('rev'):
		strand = '-'
	elif string.startswith('aa'):
		strand = '+'
	else:
		return '.', '.' #None,None
	frame = int(string[-1]) -1
	return strand, frame
					
def main():
	subcmd = sys.argv[1]
	if subcmd == 'InsertionTimePlot':
		genome = sys.argv[2]
		try: type = sys.argv[3]
		except IndexError: type = 'intact'
		try: mu = float(sys.argv[4])
		except IndexError: mu=1.3e-8
		InsertionTimePlot(genome, type, mu=mu)
	elif subcmd == 'LTRlibAnn': # hmmscan + HmmBest
		ltrlib = sys.argv[2]	# input is LTR library (fasta)
		try: 
			hmmdb = sys.argv[3] # rexdb, gydb, pfam, etc.
			try: seqtype = sys.argv[4]
			except IndexError: seqtype = 'dna'
			LTRlibAnn(ltrlib, hmmdb=hmmdb, seqtype=seqtype)
		except IndexError:
			LTRlibAnn(ltrlib)
	elif subcmd == 'HmmBest':
		inSeq = sys.argv[2] # aa seq # input: LTR library (translated protein)
		prefix = inSeq
		inHmmouts = sys.argv[3:]     # input: hmmscan output (inSeq search against hmmdb)
		hmm2best(inSeq, inHmmouts, prefix)
	elif subcmd == 'Classifier':
		gff = sys.argv[2]	# input: gff3 output by LTRlibAnn or HmmBest
		try: db = sys.argv[3]	# rexdb or gydb
		except IndexError: db = 'rexdb'
		for line in Classifier(gff, db=db):
			continue
	elif subcmd == 'replaceCls':	# LTRlibAnn + Classifier
		ltrlib = sys.argv[2]	    # input: LTR library (nucl fasta) 
		replaceCls(ltrlib)
	elif subcmd == 'replaceClsLR':
		genome = sys.argv[2]		# input: genome input for LTR_retriever pipeline
		Retriever(genome).re_classify()
	else:
		raise ValueError('Unknown command: {}'.format(subcmd))
		

	
def Args():
	parser = argparse.ArgumentParser(version=__version__)
	parser.add_argument("-s","--sequence", action="store",type=str,
					dest="input", required=True, 
					help="input TE sequences in fasta format [required]")
	parser.add_argument("-db","--hmm-database", action="store",type=str,
					dest="database name", default='rexdb', choices=DB.keys(),  
					help="the database used [default=%(default)s]")
	parser.add_argument("-st","--seq-type", action="store",type=str,
					dest="sequence type", default='nucl', choices=['nucl', 'prot'],  
					help="'nucl' for DNA and 'prot' for protein [default=%(default)s]")
	parser.add_argument("-pre", "--prefix", action="store",
					dest="prefix", default=None, type=str,
					help="output prefix [default='-s']")
	parser.add_argument("-fw", "--force-write-hmmscan", action="store_ture",
					dest="write", default=False, 
					help="if False, will use existed hmmscan outfile and skip hmmscan [default=%(default)s]")
	parser.add_argument("-p", "--processors", action="store",
					dest="processors", default=4, type=int,
					help="processors to use [default=%(default)s]")
	parser.add_argument("-tmp", "--tmp-dir", action="store",
					dest="tmp", default='./tmp', type=str,
					help="directory for temporary files [default=%(default)s]")
	parser.add_argument("-cov", "--min-coverage", action="store",
					dest="covearge", default=20, type=float,
					help="mininum covearge for protein domains in HMMScan output.[default=%(default)s]")
	parser.add_argument("-eval", "--max-evalue", action="store",
					dest="evalue", default=1e-3, type=float,
					help="maxinum E-value for protein domains in HMMScan output.[default=%(default)s]")				
	parser.add_argument("-dp2", "--disable-pass2", action="store_true",
					dest="disable", default=False, 
					help="do not to classify the unclassified sequences [default=%(default)s]")
	parser.add_argument("-nolib", "--no-library", action="store_true",
					dest="disable", default=False, 
					help="do not generate a library for RepeatMasker [default=%(default)s]")
	args = parser.parse_args()
	if args.prefix is None:
		args.prefix = '{}.{}'.format(args.sequence, args.hmm_database)
	return args
	
if __name__ == '__main__':
	main()
