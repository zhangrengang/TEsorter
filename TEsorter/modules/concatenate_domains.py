#!/usr/bin/env python
import os
import sys
import re
import random
import uuid
from Bio import SeqIO
from .RunCmdsMP import run_cmd, logger

def concat_domains(inSeq, domains, outSeq=None, tmpdir='/tmp',
		targets=None, unique=False, prefix=None, raw=False,
		order=None, superfamily=None, clade=None, subsample=None,
		format_id=False):
	'''targets: target TE id set
unique: remove duplicated sequences
prefix: provide unique file prefix, otherwise, use uuid
raw: use raw id, excluding classification
order, superfamily, clade: specify order, superfamily, clade to be included
format_id: format seq id for compatibility with iqtree'''

	d_domain = {domain: [] for domain in domains}
	if prefix is None:
		uid = uuid.uuid1()
		prefix = '{}/{}'.format(tmpdir, uid)
	
	
	# intersect
	d_seqs = {}
	d_idmap = {}
	for rc in SeqIO.parse(inSeq, 'fasta'):
		domain, *_ = re.compile(r'gene=([^;\s]+)').search(rc.description).groups()
		if domain not in d_domain:
			continue
			
		rid, cid, did = re.compile(r'^(\S+)#(\S+)#(\S+)$').match(rc.id).groups()
		if targets is not None and rid not in targets:
			continue
			
		_order, *_superfamily = cid.split('/')
		_superfamily = _superfamily[0] if _superfamily else None
		_clade, *_ = re.compile(r'clade=([^;\s]+)').search(rc.description).groups()
		if order and order != _order:
			continue
		elif superfamily and superfamily != _superfamily:
			continue
		elif clade and clade != _clade:
			continue
			
		if raw:
			raw_id = rid
		else:
			raw_id = rid + '#' + cid

		new_id = format_id_for_iqtree(raw_id) if format_id else raw_id
		d_idmap[raw_id] = new_id
		
		rc.id = new_id
		d_seqs[(raw_id, domain)] = rc
		d_domain[domain] += [raw_id]

	i = 0
	for domain, rawids in d_domain.items():
		logger.info('{} sequences contain {} domains'.format(len(rawids), domain))
		i += 1
		if i == 1:
			intersect = set(rawids)
			continue
		intersect = intersect & set(rawids)
	nseq = len(intersect)
	logger.info('{} sequences contain all {} domains'.format(len(intersect), domains))
	

	if isinstance(subsample, int) and subsample>0 and nseq > subsample:
		logger.info('Subsample {} / {} ({:.2%})'.format(subsample, nseq, subsample/nseq))
		intersect = random.sample(intersect, subsample)
	
	d_idmap = {raw_id: d_idmap[raw_id] for raw_id in intersect}
	
	# write fasta
	files = []
	for i, domain in enumerate(domains):
		outfile = '{}.{}.fa'.format(prefix, domain)
		fout = open(outfile, 'w')
		for raw_id in intersect:
			rc = d_seqs[(raw_id, domain)]
			SeqIO.write(rc, fout, 'fasta')
		fout.close()
		files += [outfile]

	# align
	alnfiles = []
	for seqfile in files:
		alnfile = seqfile + '.alignment'
		cmd = 'mafft --auto {} > {}'.format(seqfile, alnfile)
#		os.system(cmd)
		run_cmd(cmd, log=True)
		alnfiles += [alnfile]
	
	# concatenate
	outSeqfile = prefix + '.aln'
	if outSeq is None:
		outSeq = open(outSeq, 'w')
	
	catAln(alnfiles, outSeq, unique=unique)
	
	try:
		outSeq.close()
		return outSeqfile, d_idmap
	except: return None, d_idmap
	

def catAln(inALNs, outALN, unique=False):
	d_seqs = {}
	lens = []
	
	for inALN in inALNs:
		for rc in SeqIO.parse(inALN, 'fasta'):
			sp = rc.id
			seq = str(rc.seq)
			
			try: d_seqs[sp] += [seq]
			except KeyError: d_seqs[sp] = [seq]
		try: _len = len(seq)
		except UnboundLocalError: _len = 0
		lens += [_len]
	description = 'genes:{} sites:{} blocks:{}'.format(len(lens), sum(lens), lens)
	
	xseqs = set([])
	i,j = 0,0 
	for sp, seqs in list(d_seqs.items()):
		i += 1
		seqs = ''.join(seqs)
		if unique and seqs in xseqs:
			j += 1
			continue
		if unique:
			xseqs.add(seqs)
		print('>{} {}\n{}'.format(sp, description, seqs), file=outALN)
	if unique and i > 0:
		logger.info('{} ({:.1%}) unique alignments retained'.format(i-j, 1-j/i))
	if i == 0:
		logger.warn('0 sequences')

def format_id_for_iqtree(id):
    return re.compile(r'[^\w]+').sub('_', id)

def main():
	inSeq = sys.argv[1]
	domains = sys.argv[2:]
	concat_domains(inSeq=inSeq, domains=domains, outSeq=sys.stdout)

if __name__ == '__main__':
	main()
