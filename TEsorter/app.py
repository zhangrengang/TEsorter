#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Author: Zhang, Ren-Gang and Wang, Zhao-Xuan, Jacques Dainat
'''
import sys
import os
import re
import shutil
import glob
import argparse
import subprocess
import itertools
from collections import Counter, OrderedDict
from Bio import SeqIO
import logging
logging.basicConfig(level = logging.INFO,format = '%(asctime)s -%(levelname)s- %(message)s')
logger = logging.getLogger(__name__)

bindir = os.path.dirname(os.path.realpath(__file__))
sys.path = [bindir + '/bin'] + sys.path

import uuid
uid = uuid.uuid1()
default_tmpdir = './tmp-{}'.format(uid)

from .modules.translate_seq import six_frame_translate
# for multi-processing HMMScan
from .modules.RunCmdsMP import run_cmd, pp_run, pool_func
from .modules.split_records import split_fastx_by_chunk_num, cut_seqs
from .modules.small_tools import tr_numeric
#from .modules.small_tools import open_file as open
from xopen import xopen as open

# for pass-2 blast classifying
from .modules.get_record import get_records
from .modules.Blast import blast, BlastOut
from .version import __version__

DB = {
	'gydb' : bindir + '/database/GyDB2.hmm',
	'rexdb': bindir + '/database/REXdb_protein_database_viridiplantae_v4.0_plus_metazoa_v3.1.hmm',
    'rexdb-plant': bindir + '/database/REXdb_protein_database_viridiplantae_v4.0.hmm',
    'rexdb-metazoa': bindir + '/database/REXdb_protein_database_metazoa_v3.1.hmm',

	'rexdb-v3': bindir + '/database/REXdb_protein_database_viridiplantae_v3.0_plus_metazoa_v3.hmm',
	'rexdb-plantv3': bindir + '/database/REXdb_protein_database_viridiplantae_v3.0.hmm',
	'rexdb-metazoav3': bindir + '/database/REXdb_protein_database_metazoa_v3.hmm',
#	'rexdb-tir': bindir + '/database/REXdb_v3_TIR.hmm',
	'rexdb-pnas': bindir + '/database/Yuan_and_Wessler.PNAS.TIR.hmm',
	'rexdb-line': bindir + '/database/Kapitonov_et_al.GENE.LINE.hmm',
	'sine': bindir + '/database/AnnoSINE.hmm',
	}
	
ORDERS = ['LTR', 'pararetrovirus', 'DIRS', 'Penelope', 'LINE', 'SINE', 
		  'TIR', 'Helitron', 'Maverick', 'mixture', 'Unknown', 'Total']


def Args():
	parser = argparse.ArgumentParser(
			prog='TEsorter',
			description='lineage-level classification of transposable elements using conserved protein domains.',
			)
	parser.add_argument('-v', '--version', action='version',
					version='%(prog)s {version}'.format(version=__version__))
	parser.add_argument("sequence", action="store",type=str,
					help="input TE/LTR or genome sequences in fasta format [required]")
	parser.add_argument("-db","--hmm-database", action="store",type=str,
					default='rexdb', choices=list(DB.keys()),
					help="the database name used [default=%(default)s]")
	parser.add_argument("--db-hmm", action="store",type=str,
					default=None,
					help="the database HMM file used (prior to `-db`) [default=%(default)s]")
	parser.add_argument("-st","--seq-type", action="store",type=str,
					default='nucl', choices=['nucl', 'prot'],
					help="'nucl' for DNA or 'prot' for protein [default=%(default)s]")
	parser.add_argument("-pre", "--prefix", action="store",
					default=None, type=str,
					help="output prefix [default='{-s}.{-db}']")
	parser.add_argument("-fw", "--force-write-hmmscan", action="store_true",
					default=False,
					help="if False, will use the existed hmmscan outfile and skip hmmscan [default=%(default)s]")
	parser.add_argument("-p", "--processors", action="store",
					default=4, type=int,
					help="processors to use [default=%(default)s]")
	parser.add_argument("-tmp", "--tmp-dir", action="store",
					default=default_tmpdir, type=str,
					help="directory for temporary files [default=%(default)s]")
	parser.add_argument("-cov", "--min-coverage", action="store",
					default=20, type=float,
					help="mininum coverage for protein domains in HMMScan output (ranging: 0-100) [default=%(default)s]")
	parser.add_argument("-eval", "--max-evalue", action="store",
					default=1e-3, type=float,
					help="maxinum E-value for protein domains in HMMScan output (ranging: 0-10) [default=%(default)s]")
	parser.add_argument("-prob", "--min-probability", action="store",
					default=0.5, type=float,
					help="mininum posterior probability for protein domains in HMMScan output (ranging: 0-1) [default=%(default)s]")	
	parser.add_argument("-score", "--min-score", action="store",
                    default=0.1, type=float,
                    help="mininum score for protein domains in HMMScan output (ranging: 0-2) [default=%(default)s]")

	parser.add_argument("-mask", action="store", type=str, nargs='+',
					default=None, choices=['soft', 'hard'],
					help="output masked sequences (soft: masking with lowercase; hard: masking with N) [default=%(default)s]")
	parser.add_argument("-nocln", "--no-cleanup", action="store_true",
					default=False,
					help="do not clean up the temporary directory [default=%(default)s]")
	parser.add_argument("-cite", "--citation", action="store_true",
                    default=False,
                    help="print the citation and exit [default=%(default)s]")

	group_element = parser.add_argument_group('ELEMENT mode (default)', 
					'Input TE/LTR sequences to classify them into clade-level.')
	group_element.add_argument("-dp2", "--disable-pass2", action="store_true",
					default=False,
					help="do not further classify the unclassified sequences [default=%(default)s for `nucl`, True for `prot`]")
	group_element.add_argument("-rule", "--pass2-rule", action="store",
					default='80-80-80', type=str,
					help="classifying rule [identity-coverage-length] in pass-2 based on similarity [default=%(default)s]")
	group_element.add_argument("-nolib", "--no-library", action="store_true",
					default=False,
					help="do not generate a library file for RepeatMasker [default=%(default)s]")
	group_element.add_argument("-norc", "--no-reverse", action="store_true",
					default=False,
					help="do not reverse complement sequences if they are detected in minus strand [default=%(default)s]")

	group_genome = parser.add_argument_group('GENOME mode', 
					'Input genome sequences to detect TE protein domains throughout whole genome.')
	group_genome.add_argument("-genome", action="store_true",
					default=False,
					help="input is genome sequences [default=%(default)s]")
	group_genome.add_argument("-win_size",  action="store",
					default=int(270e3), type=int,
					help="window size of chunking genome sequences [default=%(default)s]")
	group_genome.add_argument("-win_ovl",  action="store",
					default=int(30e3), type=int,
					help="overlap size of windows [default=%(default)s]")
					
	args = parser.parse_args()
#	if args.prefix is None:
#		args.prefix = '{}.{}'.format(os.path.basename(args.sequence), args.hmm_database)

	if args.citation:
		print_citation()
		sys.exit()

	if args.seq_type == 'prot':
		args.disable_pass2 = True
		args.no_reverse = True
	if not args.disable_pass2:
		for key, par in zip(['p2_identity', 'p2_coverage', 'p2_length'], args.pass2_rule.split('-')):
			setattr(args, key, float(par))
#	if args.hmm_database:
#		args.seq_type = 'prot'
	return args

def print_citation():
	print('''If you use the TEsorter tool, please cite:
Zhang RG, Li GL, Wang XL et. al. TEsorter: an accurate and fast method to classify LTR retrotransposons in plant genomes [J]. Horticulture Research, 2022, 9: uhac017 [https://doi.org/10.1093/hr/uhac017]

If you use the REXdb database ('-db rexdb/rexdb-plant/rexdb-metazoa'), please cite:
Neumann P, Novák P, Hoštáková N et. al. Systematic survey of plant LTR-retrotransposons elucidates phylogenetic relationships of their polyprotein domains and provides a reference for element classification [J]. Mobile DNA, 2019, 10: 1 https://doi.org/10.1186/s13100-018-0144-1

If you use the GyDB database ('-db gydb'), please cite:
Llorens C, Futami R, Covelli L et. al. The Gypsy Database (GyDB) of mobile genetic elements: release 2.0 [J]. Nucleic Acids Research, 2011, 39: 70–74 [https://doi.org/10.1093/nar/gkq1061]

If you use the AnnoSINE database ('-db sine'), please cite:
Li Y, Jiang N, Sun Y. AnnoSINE: a short interspersed nuclear elements annotation tool for plant genomes [J]. Plant Physiology, 2022, 188: 955–970 [http://doi.org/10.1093/plphys/kiab524]

If you use the LINE/RT database ('-db rexdb-line'), please cite:
Kapitonov VV, Tempel S, Jurka J. Simple and fast classification of non-LTR retrotransposons based on phylogeny of their RT domain protein sequences [J]. Gene, 2009, 448: 207–213 [http://doi.org/10.1016/j.gene.2009.07.019]

If you use the DNA/TIR database ('-db rexdb-pnas'), please cite:
Yuan YW, Wessler SR. The catalytic domain of all eukaryotic cut-and-paste transposase superfamilies [J]. Proceedings of the National Academy of Sciences, 2011, 108: 7884–7889 [http://doi.org/10.1073/pnas.1104208108]
''')

def check_db(full_path):
#	folder_path = os.path.dirname(full_path)
#	if folder_path == '':
#		folder_path = '.'
	data_file = os.path.basename(full_path)
#	logger.info( 'db path: '+folder_path )
	logger.info( 'db file: '+ full_path )

	if not os.path.exists(full_path):
		logger.error( 'db file: '+full_path+' does not exist!' )
		sys.exit()
	else:
#		for root, dirs, files in os.walk(folder_path):
		if os.path.exists(full_path+".h3i"):
			logger.info(data_file+'\tOK')
		else:
			logger.info( 'db '+data_file+' not yet ready, building db!' )
			command = "hmmpress -f "+full_path
			#Execute the command
			process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
			stdout, stderr = process.communicate()
			logger.info(stdout.decode('utf-8'))

def pipeline(args):
	logger.info('Command: {}'.format(' '.join(sys.argv)))
	logger.info('Version: {}'.format(__version__))

	logger.info( 'VARS: {}'.format(vars(args)))
	logger.info( 'checking dependencies:' )
	db_name = args.hmm_database
	db_file = DB[args.hmm_database]
	if args.db_hmm is not None:
		db_file = args.db_hmm
		db_name = os.path.splitext(os.path.basename(db_file))[0]
#	if args.db_name is not None:
#		db_name = args.db_name
	if args.prefix is None:
		args.prefix = '{}.{}'.format(os.path.basename(args.sequence), db_name)
	Dependency().check_hmmer(db=db_file)
	if not args.disable_pass2:
		Dependency().check_blast()
	if not os.path.exists(args.tmp_dir):
		os.makedirs(args.tmp_dir)

	logger.info( 'check database '+ db_file)
	check_db(db_file)

	gap = 'X' if args.seq_type == 'prot' else 'N' 
	kargs = dict(hmmdb = db_file,
            db_name = db_name,
            prefix = args.prefix,
            force_write_hmmscan = args.force_write_hmmscan,
            processors = args.processors,
            tmpdir = args.tmp_dir,
            mincov = args.min_coverage,
            maxeval = args.max_evalue,
            minprob = args.min_probability,
			minscore = args.min_score
            )

	if args.genome:
		logger.info( 'Start identifying pipeline (GENOME mode)' )
		#print(open(args.sequence))
		#print([(rc.id, len(rc.seq)) for rc in SeqIO.parse(open(args.sequence), 'fasta')])
		seq_type = 'nucl'
		gff, geneSeq = genomeAnn(genome=args.sequence, 
			window_size=args.win_size, window_ovl=args.win_ovl, 
			seqtype = seq_type,
			**kargs
			)
		mask_gff3(args.sequence, gff, args.prefix, types=args.mask, gap=gap)
		cleanup(args)
		logger.info( 'Pipeline done.' )
		return	# genome mode stop at here
	logger.info( 'Start classifying pipeline (ELEMENT mode)' )
	lens = [len(rc.seq) for rc in SeqIO.parse(open(args.sequence), 'fasta')]
	if max(lens) > 1e6:
		logger.warn('the longest sequence is {} bp, so it may be genome. If so, \
please switch to the GENOME mode by specifiy `-genome`')
	seq_num = len(lens)
	logger.info('total {} sequences'.format(seq_num))
	# search against DB and parse
	seq_type = 'prot' if db_name == 'sine' else args.seq_type
	gff, geneSeq = LTRlibAnn(
			ltrlib = args.sequence,
			seqtype = seq_type,
			**kargs
			)

	mask_gff3(args.sequence, gff, args.prefix, types=args.mask, gap=gap)
	# classify
	classify_out = args.prefix + '.cls.tsv'
	fc = open(classify_out, 'w')
	d_class = OrderedDict()
	for rc in Classifier(gff, db=db_name, fout=fc):
		d_class[rc.id] = rc
	fc.close()
	classfied_num = len(d_class)
	logger.info('{} sequences classified by HMM'.format(classfied_num))
	logger.info('see protein domain sequences in `{}` and annotation gff3 file in `{}`'.format(geneSeq, gff))

	# pass-2 classify
	if classfied_num == 0 and not args.disable_pass2:
			logger.warning('skipping pass-2 classification for zero classification in step-1')
			args.disable_pass2 = True
	if not args.disable_pass2:
		logger.info('classifying the unclassified sequences by searching against the classified ones')
		classified_seq = '{}/pass1_classified.fa'.format(args.tmp_dir)
		unclassified_seq = '{}/pass1_unclassified.fa'.format(args.tmp_dir)
		get_records(args.sequence, classified_seq, list(d_class.keys()), type='fasta', process='get', format_id=format_gff_id)
		get_records(args.sequence, unclassified_seq, list(d_class.keys()), type='fasta', process='remove',format_id=format_gff_id)

		logger.info('using the {} rule'.format(args.pass2_rule))
		d_class2 = classify_by_blast(classified_seq, unclassified_seq,
						seqtype=args.seq_type, ncpu=args.processors,
						min_identtity=args.p2_identity, min_coverge=args.p2_coverage, min_length=args.p2_length,
						)
		fc = open(classify_out, 'a')
		for unclfed_id, clfed_id in list(d_class2.items()):
			clfed = d_class[clfed_id]
			order, superfamily, clade = clfed.order, clfed.superfamily, 'unknown'
			line = [unclfed_id, order, superfamily, clade, 'none', '?', 'none']
			print('\t'.join(line), file=fc)
			# update
			d_class[unclfed_id] = CommonClassification(*line)
		fc.close()
		logger.info('{} sequences classified in pass 2'.format(len(d_class2)))
		logger.info('total {} sequences classified.'.format(len(d_class)))
	logger.info('see classified sequences in `{}`'.format(classify_out))

	# output library
	if not args.no_library:
		out_lib = args.prefix + '.cls.lib'
		logger.info( 'writing library for RepeatMasker in `{}`'.format(out_lib) )
		fout = open(out_lib, 'w')
		for rc in SeqIO.parse(open(args.sequence), 'fasta'):
			if rc.id in d_class:
				cl = d_class[rc.id]
				strand = cl.strand
				cl = fmt_cls(cl.order, cl.superfamily, cl.clade)
				if not args.no_reverse and strand == '-':
					rc.seq = rc.seq.reverse_complement()
			else:
				cl = 'Unknown'
			rc.id = rc.id.split('#')[0] + '#' + cl
			SeqIO.write(rc, fout, 'fasta')
		fout.close()

	# rename id of protein domains
	pep_lib = args.prefix + '.cls.pep'
	logger.info( 'writing classified protein domains in `{}`'.format(pep_lib) )
	fout = open(pep_lib, 'w')
	for rc in SeqIO.parse(geneSeq, 'fasta'):
		raw_id = '|'.join(rc.id.split('|')[:-1])
		#assert raw_id in d_class
		cl = d_class[raw_id]
		cl = fmt_cls(cl.order, cl.superfamily, cl.clade)
		d_desc = dict([pair.split('=', 1)for pair in rc.description.split()[-1].split(';')])
		gene, clade = d_desc['gene'], d_desc['clade']
		new_id = '{}#{}#{}|{}'.format(raw_id.split('#')[0], cl, gene, clade)
		rc.id = new_id
		SeqIO.write(rc, fout, 'fasta')
	fout.close()
	logger.info('Summary of classifications:')
	summary(d_class)
	cleanup(args)
	logger.info( 'Pipeline done.' )
def mask_gff3(inSeq, inRM, outPrefix, types=['hard'], **kargs):
	if not types:
		return
	from .modules.RepeatMasker import mask
	for type in types:
		outSeqfile = '{}.{}masked'.format(outPrefix, type)
		soft = True if type == 'soft' else False
		logger.info( '{}-masking `{}`; output `{}`'.format(type, inSeq, outSeqfile) )
		with open(outSeqfile, 'w') as outSeq:
			masked, total = mask(inSeq, inRM, outSeq, soft=soft, **kargs)
	logger.info('{} / {} ({:.2%}) masked'.format(masked, total, masked/total))
def cleanup(args):
	# clean up
	if not args.no_cleanup:
		logger.info( 'cleaning the temporary directory {}'.format(args.tmp_dir) )
		shutil.rmtree(args.tmp_dir)

def summary(d_class):
	d_sum = {}
	for sid, clf in d_class.items():
		key = (clf.order, clf.superfamily)
		d_sum[key] = [0, 0, [], 0] # #seqs, #seqs in clades, #clades, #full domains
	for sid, clf in d_class.items():
		key = (clf.order, clf.superfamily)
		d_sum[key][0] += 1
		if clf.clade not in {'unknown', 'mixture'}:
			d_sum[key][1] += 1
			d_sum[key][2] += [clf.clade]
		if clf.completed == 'yes':
			d_sum[key][3] += 1
	
	template = '{:<16}{:<16}{:>15}{:>15}{:>15}{:>15}'
	line = ['Order', 'Superfamily', '# of Sequences', '# of Clade Sequences', '# of Clades', '# of full Domains']

	print(template.format(*line), file=sys.stdout)
	for (order, superfamliy), summary in \
			sorted(list(d_sum.items()), key=lambda x: (ORDERS.index(x[0][0]), x[0][1])):
		line = [order, superfamliy, summary[0], summary[1], len(set(summary[2])), summary[3]]
		line = list(map(str, line))
		line = template.format(*line)
		print(line, file=sys.stdout)

def fmt_cls(*args):
	values = []
	for arg in args:
		if arg == 'unknown' or arg in set(values):
			continue
		values += [arg]
	return '/'.join(values)

class CommonClassification(object):
	def __init__(self, id=None, order=None, superfamily=None,
				clade=None, completed=None, strand=None, domains=None):
		self.id = id
		self.order = order
		self.superfamily = superfamily
		self.clade = clade
		self.completed = completed
		self.strand = strand
		self.domains = domains
class CommonClassifications:
	def __init__(self, clsfile):
		self.clsfile = clsfile
	def __iter__(self):
		return self._parse()
	def _parse(self):
		for i, line in enumerate(open(self.clsfile)):
			if i == 0:
				continue
			line = line.strip().split('\t')
			yield CommonClassification(*line)

def classify_by_blast(db_seq, qry_seq, blast_out=None, seqtype='nucl', ncpu=4,
					  min_identtity=80, min_coverge=80, min_length=80):
	'''pass-2 classify'''
	if os.path.getsize(db_seq) == 0:
		return {}

	blast_outfmt = "'6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs qcovhsp sstrand'"
	blast_out = blast(db_seq, qry_seq, seqtype=seqtype, blast_out=blast_out, blast_outfmt=blast_outfmt, ncpu=ncpu)
	with open(blast_out+'.best', 'w') as fb:
		d_best_hit = BlastOut(blast_out, blast_outfmt).filter_besthit(fout=fb)
		d_best_hit_copy = d_best_hit.copy()
	for qseqid, rc in d_best_hit_copy.items():
		if not (rc.pident >= min_identtity and rc.qcovs >= min_coverge and rc.length >= min_length):
			del d_best_hit[qseqid]
	d_class = OrderedDict([(qseqid, rc.sseqid) for qseqid, rc in d_best_hit.items()])
	return d_class


class Classifier(object):
	def __init__(self, gff=None, db='rexdb', fout=sys.stdout): # gff is sorted
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
		line = ['#TE', 'Order', 'Superfamily', 'Clade', 'Complete', 'Strand', 'Domains']
		print('\t'.join(line), file=self.fout)
		for rc in self.parse():
			rc_flt = rc
			strands = [line.strand for line in rc_flt]
			strands_num = len(set(strands))
			if strands_num >1 :
				strand = '?'
			elif strands_num == 1:
				strand = strands[0]
			else:
				continue
			if strand == '-':
				rc_flt.reverse()
			lid = rc_flt[0].ltrid
			domains = ' '.join(['{}|{}'.format(line.gene, line.clade)  for line in rc])
			order, superfamily, max_clade, coding = self.classify_element(rc_flt)
			line = [lid, order, superfamily, max_clade, coding, strand, domains]
			print('\t'.join(line), file=self.fout)
			yield CommonClassification(*line)
	def classify_element(self, lines):
		genes  = [line.gene  for line in lines]
		clades = [line.clade for line in lines]
		names = [line.name for line in lines]
		if self.db.startswith('rexdb'):
			order, superfamily, max_clade, coding = self.identify_rexdb(genes, names)
		elif self.db == 'gydb':
			order, superfamily, max_clade, coding = self.identify_gydb(genes, clades)
		elif self.db == 'sine':
			order, superfamily, max_clade, coding = 'SINE', 'unknown','unknown','unknown'
		else:
			order, superfamily, max_clade, coding = 'Unknown', 'unknown','unknown','unknown'
		return order, superfamily, max_clade, coding
	def identify_rexdb(self, genes, clades):
		perfect_structure = {
#            ('LTR', 'Copia'): ['Ty1-GAG', 'Ty1-PROT', 'Ty1-INT', 'Ty1-RT', 'Ty1-RH'],
#            ('LTR', 'Gypsy'): ['Ty3-GAG', 'Ty3-PROT', 'Ty3-RT', 'Ty3-RH', 'Ty3-INT'],
			('LTR', 'Copia'): ['GAG', 'PROT', 'INT', 'RT', 'RH'],
			('LTR', 'Gypsy'): ['GAG', 'PROT', 'RT', 'RH', 'INT'],
			('LTR', 'Bel-Pao'): ['GAG', 'PROT', 'RT', 'RH', 'INT'],
			}
		clade_count = Counter(clades)
		counts = list(clade_count.values())
		max_clade = max(clade_count, key=lambda x: clade_count[x])
		order, superfamily = self._parse_rexdb(max_clade)
		if len(clade_count) == 1 or (clade_count[max_clade] > 1 and counts[0] > counts[1]):
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
				coding = 'yes' # completed gene structure
			else:
				coding = 'no'
		except KeyError:
			coding = 'unknown'
		if superfamily not in {'Copia', 'Gypsy'}:
			max_clade = 'unknown'
		if max_clade.startswith('Ty'): # Ty3_gypsy, Ty1_copia, Ty1-outgroup in metazoa_v3
			max_clade = 'unknown'
		return order, superfamily, max_clade, coding
	def _parse_rexdb(self, clade): # full clade name
		if clade.startswith('Class_I/LTR/Ty1_copia'):
			order, superfamily = 'LTR', 'Copia'
		elif clade.startswith('Class_I/LTR/Ty3_gypsy'):
			order, superfamily = 'LTR', 'Gypsy'
		elif clade.startswith('Class_I/LTR/'): # LTR/Bel-Pao, LTR/Retrovirus
			order, superfamily = clade.split('/')[1:3]
		elif clade.startswith('Class_I/'): # LINE, pararetrovirus, Penelope, DIRS
			try: order, superfamily = clade.split('/')[1:3]
			except ValueError: order, superfamily = clade.split('/')[1], 'unknown'
		elif clade.startswith('Class_II/'): # TIR/hAT, Helitro, Maverick
			try: order, superfamily = clade.split('/')[2:4]
			except ValueError: order, superfamily = clade.split('/')[2], 'unknown'
		elif clade.startswith('NA'): # "NA:Retrovirus-RH"
			order, superfamily = 'LTR', 'Retrovirus'
		else:	# not get it
			logger.warning( 'unknown clade: {}'.format(max_clade) )
		return order, superfamily
	def identify_gydb(self, genes, clades):
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
			logger.warning( 'unknown clade: {}'.format(max_clade) )
		if len(clade_count) == 1 or clade_count[max_clade] > 1:
			max_clade = max_clade
		elif len(clade_count) > 1:
			max_clade = 'mixture'
			superfamlies = [d_map.get(clade, [None,None])[1] for clade in clades]
			if len(Counter(superfamlies)) > 1:
				superfamily = 'mixture'
				orders = [d_map.get(clade, [None,None])[0] for clade in clades]
				if len(Counter(orders)) > 1:
					order = 'mixture'

		try:
			ordered_genes = perfect_structure[(order, superfamily)]
			my_genes = [gene for gene in genes if gene in set(ordered_genes)]
			if ordered_genes == my_genes:
				coding = 'yes' # completed gene structure and the same order
			else:
				coding = 'no'
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
					logger.warn( 'skipped KeyError: {}'.format(e) )
			if intid in d_class:
				neword, newfam = d_class[intid]
				re_org = self.re_orgnize(rc.id, neword, newfam)
				if re_org:
					i += 1
					rc.id = re_org
			SeqIO.write(rc, fout, 'fasta')
		logger.info( 'sequences re-classified'.format(i) )
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
		self.clade_map = {
			'Ty_(Pseudovirus)': 'pseudovirus',
			'Cer2-3': 'cer2-3',
			'412/Mdg1': '412_mdg1',
			'TF1-2': 'TF',
			'Micropia/Mdg3': 'micropia_mdg3',
			'CoDi-I': 'codi_I',
			'CoDi-II': 'codi_II',
			'17.6': '17_6',
			}

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
			self.dict = dict(list(zip(title, temp)))
			if self.dict['Clade'] == 'NA':
				self.clade = self.dict['Cluster_or_genus']
			else:
				self.clade = self.dict['Clade']
			self.superfamily = self.dict['Family'].split('/')[-1]
			if self.superfamily == 'Retroviridae':	# deltaretroviridae gammaretroviridae
				self.clade = self.dict['Cluster_or_genus'].replace('virus', 'viridae')
			if self.superfamily == 'Retrovirus':	# an exception
				self.superfamily = 'Retroviridae'
			self.order = 'LTR' if self.dict['System'] in {'LTR_retroelements', 'LTR_Retroelements', 'LTR_retroid_elements'} else self.dict['System']
			yield self
			if self.clade in self.clade_map:
				self.clade = self.clade_map[self.clade]
				yield self
				if self.clade == '412_mdg1':
					self.clade = '412-mdg1'  # 412-mdg1 and 412_mdg1
					yield self

			self.clade = self.clade.replace('-', '_') # A-clade V-clade C-clade
			yield self
			self.clade = self.clade.lower()
			yield self
		self.order, self.superfamily, self.clade, self.dict = ['LTR', 'Copia', 'ty1/copia', {}]  # AP_ty1/copia
		yield self
		order_map = { 	# some unknown clade
			'retroelement': 'LTR',
			'retroviridae': 'LTR',
			'B-type_betaretroviridae': 'LTR',
			'D-type_betaretroviridae': 'LTR',
			'caulimoviruses': 'LTR',
			'caulimoviridae_dom2': 'LTR',
			'errantiviridae': 'LTR',
			'retropepsins': 'LTR',
			'VPX_retroviridae': 'LTR',
			'cog5550': 'Unknown',
			'ddi': 'Unknown',
			'dtg_ilg_template': 'Unknown',
			'saspase': 'Unknown',
			'GIN1': 'Unknown',
			'shadow': 'Unknown',
			'all': 'Unknown',
			'pepsins_A1a': 'Unknown',
			'pepsins_A1b': 'Unknown',
			}
		for clade, order in list(order_map.items()):
			self.order, self.superfamily, self.clade, self.dict = [order, 'unknown', clade, {}]  # CHR_retroelement
			yield self

class GffLine(object):
	def __init__(self, line):
		if isinstance(line, str):
			temp = line.strip().split('\t')
		else:
			temp = line
		if len(temp) == 8:
			temp = list(temp) + ['']
		self.chr, self.source, self.type, self.start, self.end, self.score, self.strand, \
			self.frame, self.attributes = temp
		self.start, self.end = int(self.start), int(self.end)
		try: self.score = float(self.score)
		except: pass
		try: self.frame = int(self.frame)
		except: pass
		try: self.attributes = self.parse(self.attributes)
		except: pass
	def parse(self, attributes):
		return dict(kv.split('=', 1) for kv in attributes.split(';') if kv)

class LTRgffLine(GffLine):
	def __init__(self, line):
		super(LTRgffLine, self).__init__(line)
		self.gene = self.attributes['gene']
		self.clade = self.attributes['clade']
		self.ltrid = '|'.join(self.attributes['ID'].split('|')[:-1])
		self.name = self.attributes['ID'].split('|')[-1].split(':')[0]

class HmmScan(object):
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

class HmmDomRecord(object):
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
			list(map(tr_numeric, [self.tlen, self.qlen, self.domi, self.domn, \
					self.hmmstart, self.hmmend, self.alnstart, self.alnend, self.envstart, self.envend]))
		self.evalue, self.score, self.bias, self.cevalue, self.ievalue, self.domscore, self.dombias, self.acc = \
			list(map(tr_numeric, [self.evalue, self.score, self.bias, self.cevalue, self.ievalue, self.domscore, self.dombias, self.acc]))
		self.tdesc = ' '.join(temp[22:])
		
	@property
	def hmmcov(self):
		return round(1e2*(self.hmmend - self.hmmstart + 1) / self.tlen, 1)



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
	elif db.startswith('sine'):
		gene = 'SINE' #hmmname
		clade = 'SINE'
	else:
		gene, clade = hmmname, hmmname
	return gene, clade
class HmmCluster(object):
	def __int__(self, hmmout, seqtype = 'nucl'): # only for nucl
		self.hmmout = hmmout
		self.seqtype = seqtype

	def cluster(self):
		d_cluster = OrderedDict()
		for rc in HmmScan(self.hmmout):
			suffix = rc.qname.split('|')[-1]
			qid = '|'.join(rc.qname.split('|')[:-1])
			strand, frame = parse_frame(rc.qname.split('|')[-1])
			rc.qid, rc.strand, rc.frame = qid, strand, frame
			key = (rc.tname, qid, strand)
			try: d_cluster[key] += [rc]
			except KeyError: d_cluster[key] = [rc]

		return d_cluster

	def extand(self, records, max_mal=5):
		if len(records) == 1:
			return records
		records = sorted(records, key=lambda x:x.hmmstart)
		best_idx, best_rc = self.maxscore(records)
		new_rcs = [best_rc]
		for rc in records[:best_idx][::-1]:		# <<- left extand
			right_rc = new_rcs[0]
			mal_pos = right_rc.hmmstart - rc.hmmend
			if abs(mal_pos) <= max_mal:
				if mal_pos <= 0:
					diff = 1-mal_pos
					rc.hmmend -= diff
					rc.alnend -= diff
				new_rcs = [rc] + new_rcs
		for rc in records[best_idx:]:			# ->> right extand
			left_rc = new_rcs[-1]
			mal_pos = rc.hmmstart - left_rc.hmmend
			if abs(mal_pos) <= max_mal:
				if mal_pos <= 0:
					diff = 1-mal_pos
					rc.hmmstart += diff
					rc.alnstart += diff
				new_rcs += [rc]
		return HmmClusterRecord(new_rcs)

	def maxscore(self, records):
		for i, rc in enumerate(records):
			if i == 0:
				best = (i, rc)
				continue
			if rc.score > best[1].score:
				best = (i, rc)
		return best

class HmmClusterRecord(object):
	def __init__(self, records):
		self.records = records
		self.score = sum([rc.score for rc in records])
		self.hmmstart = records[0].hmmstart
		self.hmmend = records[-1].hmmend
		self.alnstart = records[0].alnstart
		self.alnend = records[-1].alnend
		self.tlen = records[0].tlen
		self.hmmcov = round(1e2*(self.hmmend - self.hmmstart + 1) / self.tlen, 1)
		self.evalue = multi(*[rc.evalue for rc in records])

def multi(*n):
	result = 1
	for i in n:
		result = result * i
	return result

def group_resolve_overlaps(lines):
	'''assume multiple chromsomes'''
	resolved_lines = []
	for chrom, items in itertools.groupby(lines, key=lambda x:x[0]):
		logger.info('resolving overlaps in {}'.format(chrom))
		resolved_lines += resolve_overlaps(list(items))
	return resolved_lines
def overlap(self, other):
	# gff line: gffline = [qid, 'TEsorter', 'CDS', nuc_start, nuc_end, rc.score, strand, frame, attr
	ovl = max(0, min(self[4], other[4]) - max(self[3], other[3]))
	return 100*ovl/(min((self[4]-self[3]+1), (other[4]-other[3]+1)))
	
def resolve_overlaps(lines, max_ovl=20, ):
	'''assume only one chromsome'''
	last_line = None
	discards = []
	ie, io = 0, 0
	for line in sorted(lines, key=lambda x:x[3]):
		discard = None
	#	print(last_line, line)
		if last_line:
			if line == last_line:	# equal
				ie += 1
				line_pair = [last_line, line]	# retain, discard
			else:
				if overlap(line, last_line) > max_ovl:
					io += 1
					if line[5] > last_line[5]:	# score is prior
						line_pair = [line, last_line]
					else: 
						line_pair = [last_line, line]
				else:	# no overlap or too short overlap
					last_line = line
					continue
			
			retain, discard = line_pair

			discards += [discard]

		if not last_line or discard != line:
			last_line = line
#	logger.info('discard {} equal and {} overlapped hits; {} in total'.format(ie, io, ie+io))
	return sorted(set(lines) - set(discards), key=lambda x:x[3])

def _hmm2best(inHmmouts, db='rexdb', seqtype='nucl', genome=False):
	'''best HMM hit based on score'''
	d_besthit = {}
	for inHmmout in inHmmouts:
		for rc in HmmScan(inHmmout):
			suffix = rc.qname.split('|')[-1]
			if seqtype == 'nucl' and (suffix.startswith('aa') or suffix.startswith('rev_aa')):
				qid = '|'.join(rc.qname.split('|')[:-1])
			else:
				qid = rc.qname
			domain,clade = parse_hmmname(rc.tname, db=db)
			key = (qid,)
			if genome:
				key += (rc.envstart, rc.envend)
			# normlize score
			rc.score = round(rc.domscore / rc.tlen, 2)
			rc.evalue = rc.ievalue
				
			if db.startswith('rexdb'):
				cdomain = domain.split('-')[1]
				if cdomain == 'aRH' and not genome:
					cdomain = 'RH'
				if cdomain == 'TPase' and not genome:
					cdomain = 'INT'
				key += (cdomain,)
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
			else: # gydb
				key += (domain,)
				if key in d_besthit:
					if rc.score > d_besthit[key].score:
						d_besthit[key] = rc
				else:
					d_besthit[key] = rc
	return d_besthit
	
def seqs2dict(inSeqs):
	d = {}
	for inSeq in inSeqs:
		for rc in SeqIO.parse(inSeq, 'fasta'):
			d[rc.id] = rc
	return d
	#return dict([(rc.id, rc) for rc in SeqIO.parse(inSeq, 'fasta') for inSeq in inSeqs])

def format_gff_id(id):
#	return id.replace(';', '_').replace('=', '_')
	return re.compile(r'[;=\|]').sub("_", id)

def hmm2best(inSeqs, inHmmouts, nucl_len=None, prefix=None, db='rexdb', seqtype='nucl', 
			mincov=20, maxeval=1e-3, minprob=0.6, minscore=0.5, genome=False):
	if prefix is None:
		prefix = inSeqs[0]
	if nucl_len is None and seqtype=='nucl':
		raise ValueError('Sequences length must provide for `{}` sequences'.format(seqtype))
	d_besthit = _hmm2best(inHmmouts, db=db, seqtype=seqtype, genome=genome)
	
	d_seqs = seqs2dict(inSeqs)
	lines = []
	for key, rc in list(d_besthit.items()):
		if genome:
			qid, s, e, domain = key
		else:
			qid, domain = key
		if rc.hmmcov < mincov or rc.evalue > maxeval or rc.acc < minprob or rc.score < minscore:
			continue
		rawid = qid
		gene,clade = parse_hmmname(rc.tname, db=db)
		if db.startswith('rexdb'):
			domain = gene.split('-')[1]
		gid = '{}|{}'.format(format_gff_id(qid), rc.tname)
		#gid = '{}|{}'.format(qid, rc.tname)
		try: gseq = d_seqs[rc.qname].seq[rc.envstart-1:rc.envend]
		except KeyError as e:
			raise KeyError('{}\nIt seems that the HMM domtbl file is not consistent with the sequence file. Please retry with `-fw`.')
		gseq = str(gseq)
		if seqtype == 'nucl':
			strand, frame = parse_frame(rc.qname.split('|')[-1])
			if strand == '+':
				nuc_start = (rc.envstart-1) * 3  + frame + 1
				nuc_end = rc.envend * 3 + frame
			elif strand == '-':
				nucl_length = nucl_len[qid]
				nuc_start = nucl_length - (rc.envend * 3 + frame) + 1
				nuc_end = nucl_length - ((rc.envstart-1) * 3 + frame)
			else:	# not translated
				nuc_start = rc.envstart
				nuc_end = rc.envend
		elif seqtype == 'prot':
			strand, frame = '+', '.'
			nuc_start, nuc_end, = rc.envstart, rc.envend
		match = re.compile(r'(\S+?):(\d+)[\.\-]+(\d+)').match(qid)
		if match:
			qid, linestart, ltrend = match.groups()
			linestart = int(linestart)
			nuc_start = linestart + nuc_start - 1
			nuc_end = linestart + nuc_end -1
		_add = ''
		gffline = (qid, 'TEsorter', 'CDS', nuc_start, nuc_end, rc.score, strand, frame, )
		if genome:
			gid = '{}:{}-{}|{}'.format(qid, nuc_start, nuc_end, rc.tname)
			element = LTRgffLine(gffline + ({'ID':gid, 'gene':domain, 'clade':clade},))
			order, superfamily, max_clade, coding = Classifier(db=db).classify_element([element])
#			if order == 'Unknown':
#				logger.warn('unknown element: {}, is excluded'.format(gid))
#				continue
			cls = fmt_cls(order, superfamily, max_clade)
			nstop = list(gseq).count('*')
			match = '{} {} {}'.format(rc.tname, rc.hmmstart, rc.hmmend)
			_add = 'Classification={};Target={};nstop={};'.format(cls, match, nstop)
		name = '{}-{}'.format(clade, domain)
		attr = 'ID={};Name={};{}gene={};clade={};coverage={};evalue={};probability={}'.format(
				gid, name, _add, domain, clade, rc.hmmcov, rc.evalue,  rc.acc)
		gffline = gffline + (attr, rc.evalue, rc.hmmcov, rc.acc, rawid, gid, gseq)
		lines.append(gffline)
	
	gff, seq, tsv = '{}.dom.gff3'.format(prefix), '{}.dom.faa'.format(prefix), '{}.dom.tsv'.format(prefix)
	
	# with open(gff+'.debug', 'w') as f:
		# for line in sorted(lines, key=lambda x: (x[0], x[-3], x[3])):
			# gffline = line[:9]
			# gffline = list(map(str, gffline))
			# print('\t'.join(gffline), file=f)
			
	if genome:
		lines = group_resolve_overlaps(sorted(lines, key=lambda x: x[0]))
	else:
		lines = sorted(lines, key=lambda x: (x[0], x[-3], x[3]))
	
	fgff = open(gff, 'w')
	fseq = open(seq, 'w')
	ftsv = open(tsv, 'w')
	print('\t'.join(['#id', 'length', 'evalue', 'coverge', 'probability', 'score']), file=ftsv)
	for line in lines:
		gffline = line[:9]
		gffline = list(map(str, gffline))
		print('\t'.join(gffline), file=fgff)
		gid, gseq = line[-2:]
		gdesc = line[8]
		print('>{} {}\n{}'.format(gid, gdesc, gseq), file=fseq)
		evalue, hmmcov, acc = line[-6:-3]
		score = gffline[5]
		line = [gid, len(gseq), evalue, hmmcov, acc, score]
		print('\t'.join(map(str, line)), file=ftsv)
	fgff.close()
	fseq.close()
	return gff, seq
def summary_genome(gff, fout=sys.stdout):
	last_chr, last_end = '', 0
	d_stats = {}
	for line in open(gff):
		line = GffLine(line)
		cls = line.attributes['Classification']
		cls = tuple(cls.split('/'))
		assert len(cls) <=3
		assert cls[0] in set(ORDERS), 'Unknown order: {}'.format(cls)
		if line.end < last_end:
			continue
		elif line.start < last_end:
			start = last_end+1
		else:
			start = line.start
		_len = line.end - start + 1
		xcls = []
		if len(cls) ==3:
			xcls += [cls[:1], cls[:2]]
		elif len(cls) ==2:
			xcls += [cls[:1]]
		xcls += [cls] + [('Total', )]
		for cls in xcls:
			try: d_stats[cls] += [_len]
			except KeyError: d_stats[cls] = [_len]
	line = ['Order', 'Superfamily', 'Clade', 'Number', 'Total_length', 'Mean_length']
	fout.write('\t'.join(line)+'\n')
	for cls, lens in sorted(d_stats.items(), key=lambda x:(ORDERS.index(x[0][0]), x[0])):
		cls = ['']*(len(cls)-1) + [cls[-1]] + ['']* (3-len(cls))
		n, total = len(lens), sum(lens)
		mean = round(total/n, 1)
		line = cls + [n, total, mean]
		line = map(str, line)
		fout.write('\t'.join(line)+'\n')
def translate(inSeq, prefix=None, overwrite=True):
	if prefix is None:
		prefix = inSeq
	outSeq = prefix + '.aa'
	overwrite = not (os.path.exists(outSeq) and os.path.getsize(outSeq) >0) or overwrite
	if not overwrite:
		logger.info( 'use existed non-empty `{}` and skip translating'.format(outSeq) )
		return outSeq
	logger.info( 'translating `{}` in six frames'.format(inSeq) )
	with open(outSeq, 'w') as fp:
		six_frame_translate(inSeq, fp)
	return outSeq

def translate_pp(inSeq, prefix=None, tmpdir='./tmp', processors=4):
	if prefix is None:
		prefix = inSeq
	outSeq = prefix + '.aa'
	chunk_prefix = '{}/{}'.format(tmpdir, 'chunk_nuclseq')
	_, _, _, chunk_files = split_fastx_by_chunk_num(
			inSeq, prefix=chunk_prefix, chunk_num=processors, seqfmt='fasta', suffix='')

def hmmscan(inSeq, hmmdb='rexdb.hmm', hmmout=None, ncpu=4, bin='hmmscan'):
	if hmmout is None:
		hmmout = prefix + '.domtbl'
	cmd = '{} --nobias --notextw --noali --cpu {} --domtblout {} {} {} > /dev/null'.format(bin, 
			ncpu, hmmout, hmmdb, inSeq)
	run_cmd(cmd, logger=logger)
	return hmmout
def _translate(arg):
	inSeq, overwrite = arg
	return translate(inSeq, overwrite=overwrite)

def hmmscan_pp(inSeq, hmmdb='rexdb.hmm', hmmout=None, tmpdir='./tmp', processors=4, 
			bin='hmmscan',seqtype='nucl', force_write_hmmscan=False):
	if hmmout is None:
		hmmout = prefix + '.domtbl'
	overwrite = not (os.path.exists(hmmout) and os.path.getsize(hmmout) >0) or force_write_hmmscan
	
	chunk_prefix = '{}/{}'.format(tmpdir, 'chunk')
	if processors > 1:
		chunk_num = processors*2
		_, _, _, chunk_files = split_fastx_by_chunk_num(
			inSeq, prefix=chunk_prefix, chunk_num=chunk_num, seqfmt='fasta', suffix='')
		chunk_files = [chunk_file for chunk_file in chunk_files if os.path.getsize(chunk_file)>0]
	else:
		chunk_files = [inSeq]
	if seqtype == 'nucl':	#translate
		iterable = ((chunk, overwrite) for chunk in chunk_files)
		chunk_files = list(pool_func(_translate, iterable, processors=processors))
	
	if not overwrite:
		logger.info( 'use existed non-empty `{}` and skip hmmscan'.format(hmmout) )
		return chunk_files
	
	domtbl_files = [chunk_file + '.domtbl' for chunk_file in chunk_files]
	cmds = [
		'{} --nobias --notextw --noali --domtblout {} {} {} > /dev/null'.format(
			bin, domtbl_file, hmmdb, chunk_file) \
			for chunk_file, domtbl_file in zip(chunk_files, domtbl_files)]
	jobs = pp_run(cmds, processors=processors)
	for cmd, (stdout, stderr, status) in zip(cmds, jobs):
		if not status == 0:
			logger.warning( "exit code {} for CMD '{}'".format(status, cmd) )
			logger.warning('\n\tSTDOUT:\n{0}\n\tSTDERR:\n{1}\n\n'.format(stdout, stderr))
	# cat files
	
	with open(hmmout, 'w') as f:
		for domtbl_file in domtbl_files:
			for line in open(domtbl_file):
				f.write(line)
	return chunk_files
def genomeAnn(genome, tmpdir='./tmp', seqfmt='fasta',window_size=1e6, window_ovl=1e5, **kargs):
	cutSeq = '{}/cut.{}'.format(tmpdir, seqfmt)
	with open(cutSeq, 'w') as f:
		cut_seqs(genome, f, window_size=window_size, window_ovl=window_ovl, seqfmt=seqfmt)
	gff, geneSeq = LTRlibAnn(ltrlib=cutSeq, genome=True, tmpdir=tmpdir, **kargs)
	logger.info('Summary of classifications:')
	summary_genome(gff, fout=sys.stdout)
	return gff, geneSeq
	
def LTRlibAnn(ltrlib, hmmdb='rexdb', db_name='rexdb',  seqtype='nucl', prefix=None,
			force_write_hmmscan=False, genome=False, 
			processors=4, tmpdir='./tmp',
			mincov=20, maxeval=1e-3, minprob=0.5, minscore=0.5):
	if prefix is None:
		prefix = '{}.{}'.format(ltrlib, hmmdb)
	bin = 'hmmscan'
	domtbl = prefix + '.domtbl'
	if seqtype == 'nucl' :
		d_nucl_len = dict([(rc.id, len(rc.seq)) for rc in SeqIO.parse(open(ltrlib), 'fasta')])
	elif seqtype == 'prot':
		d_nucl_len = None

	logger.info( 'HMM scanning against `{}`'.format(hmmdb) )
	
	chunk_files = hmmscan_pp(ltrlib, hmmdb=hmmdb, hmmout=domtbl, tmpdir=tmpdir, 
			processors=processors, bin=bin, seqtype=seqtype, force_write_hmmscan=force_write_hmmscan)
	logger.info( 'generating gene anntations' )
	gff, geneSeq = hmm2best(chunk_files, [domtbl], db=db_name, nucl_len=d_nucl_len, genome=genome,
				prefix=prefix, seqtype=seqtype, mincov=mincov, maxeval=maxeval, minprob=minprob, minscore=minscore)
	return gff, geneSeq

def replaceCls(ltrlib, seqtype='nucl', db='rexdb'):
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

def parse_frame(string):	# frame=0-2
	if string.startswith('rev'):
		strand = '-'
	elif string.startswith('aa'):
		strand = '+'
	else:
		return '.', '.' #None,None
	frame = int(string[-1]) -1
	return strand, frame

class Dependency(object):
	def __init__(self,):
		pass

	def check(self):
		pass

	def check_hmmer(self, db, program='hmmscan'):
		dp_version = self.get_hmm_version(db)[:3]
		if self.check_presence(program):
			version0 = self.check_hmmer_verion(program)
			version = version0[:3]
			if version >= dp_version:
				logger.info('hmmer\t{}\tOK'.format(version0))
			elif version < dp_version:
				logger.warning('hmmer version {} is too low. Please update to {} from http://hmmer.org/download.html'.format(version, dp_version))
			else:
				logger.info('hmmer version {} is too high. You may use the version {}. However, I update the database first.'.format(version, dp_version))
				self.update_hmmer(db)
		else:
			logger.error('hmmer>={} not found'.format(dp_version))

	def get_hmm_version(self, db):
		line = open(db).readline()
		version = re.compile(r'HMMER\S+ \[([\w\.]+)').search(line).groups()[0]
		return version

	def update_hmmer(self, db):
		from small_tools import backup_file
		bk_db, db = backup_file(db)
		for suffix in ['.h3f', '.h3i', '.h3m', '.h3p']:
			backup_file(bk_db + suffix)
		cmd = 'hmmconvert {} > {}'.format(bk_db, db)
		out, err, status0 = run_cmd(cmd, logger=logger)
		cmd = 'hmmpress {}'.format(db)
		out, err, status1 = run_cmd(cmd, logger=logger)
		if status0 + status1 == 0:
			logger.info('HMM converted. it will continue')
		else:
			logger.error('HMM failed to convert. exit')
			sys.exit(1)

	def check_blast(self, program='blastn'):
		if self.check_presence(program):
			version = self.check_blast_version(program)
			logger.info('{}\t{}\tOK'.format(program, version))
		else:
			logger.error('{} not found'.format(program))

	def check_presence(self, program):
		cmd = 'which {}'.format(program)
		out, err, status = run_cmd(cmd)
		if status == 0:
			return True
		else:
			return False

	def check_hmmer_verion(self, program):
		cmd = '{} -h'.format(program)
		out, err, status = run_cmd(cmd)
		version = re.compile(r'HMMER (\S+)').search(out.decode('utf-8')).groups()[0]
		return version

	def check_blast_version(self, program):
		cmd = '{} -version'.format(program)
		out, err, status = run_cmd(cmd)
		version = re.compile(r'blast\S* ([\d\.\+]+)').search(out.decode('utf-8')).groups()[0]
		return version

def main():
	pipeline(Args())

if __name__ == '__main__':
	main()
