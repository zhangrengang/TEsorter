import sys,os
import re
from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter

class RMOutRecord():
	def __init__(self, line):
		d_class = dict([
			('LTR', 'ClassI'),
			('LINE', 'ClassI'),
			('SINE', 'ClassI'),
			('Retroposon', 'ClassI'),
			('DNA', 'ClassII'),
			('RC', 'ClassII'),
			('Unknown', 'Others'),
			('Satellite', 'Others'),
			('Simple_repeat', 'Others'),
			('Low_complexity', 'Others'),
			('rRNA', 'Others'),
			('snRNA', 'Others'),
			])

		temp = line.strip().split()
		title = ['score', 'perc_div', 'perc_del', 'perc_ins',
				 'query_id', 'query_begin', 'query_end', 'query_left', 'strand',
				 'repeat_family', 'super_class', 'repeat_begin', 'repeat_end', 'repeat_left', 'repeat_id',
				]
		convert = [int, float, float, float,
				   str, int, int, str, str,
				   str, str, str, int, str, str ]
		assert len(temp) == len(title) or (len(temp)-1 == len(title) and temp[-1] == "*") or (len(temp) == len(title)-1)
#		try: assert len(temp) == len(title) or (len(temp)-1 == len(title) and temp[-1] == "*")
#		except AssertionError: print >> sys.stderr, temp, '\n', title
		self.__dict__ = {key: func(value) for key,value,func in zip(title, temp, convert)}
		self.query_left = int(self.query_left.strip('()'))
		self.repeat_begin = int(self.repeat_begin.strip('()'))
		self.repeat_left = int(self.repeat_left.strip('()'))
		if self.strand == 'C':
			self.strand = '-'
		if temp[-1] == "*": # there is a higher-scoring match whose domain partly (<80%) includes the domain of this match
			self.overlap = True
		else:
			self.overlap = False
		self.order = self.super_class.split('/')[0]
#		self.order = self.order.rstrip('?')
		try: self.superfamily = self.super_class.split('/')[1]
		except IndexError: self.superfamily = ''
		try: self.Class = d_class[self.order.rstrip('?')]
		except KeyError: self.Class = 'Others'
		self.query_length = self.query_end + self.query_left
		self.query_match_length = self.query_begin - self.query_end + 1
		if self.strand == '-':
			self.repeat_begin, self.repeat_end, self.repeat_left = self.repeat_left, self.repeat_end, self.repeat_begin 
		self.target = '{} {} {}'.format(self.repeat_family, self.repeat_begin, self.repeat_end)
	def write(self, fout=sys.stdout):
		if self.strand == '+':
			line = [self.score, self.perc_div, self.perc_del, self.perc_ins, \
				self.query_id, self.query_begin, self.query_end, '({})'.format(self.query_left), self.strand, \
				self.repeat_family, self.super_class, \
				self.repeat_begin, self.repeat_end, '({})'.format(self.repeat_left), self.repeat_id]
		elif self.strand == '-':
			line = [self.score, self.perc_div, self.perc_del, self.perc_ins, \
				self.query_id, self.query_begin, self.query_end, '({})'.format(self.query_left), self.strand, \
				self.repeat_family, self.super_class, \
				'({})'.format(self.repeat_left), self.repeat_end, self.repeat_begin, self.repeat_id]
		if self.overlap:
			line += ['*']
		line = map(str, line)
		print >> fout, '\t'.join(line)
	def write_gff3(self, fout=sys.stdout):
		attr = 'Target={target}'.format(self.__dict__)
		line = [self.query_id, 'RepeatMasker', 'dispersed_repeat', self.query_begin, self.query_end, \
				self.score, self.strand, '.', attr]
		line = map(str, line)
		print >> fout, '\t'.join(line)
	def get_seq(self, seqRecord):
		id = '{}:{}..{}|{}#{}'.format(self.query_id, self.query_begin, self.query_end, self.repeat_family, self.super_class)
		teRecord = seqRecord[self.query_begin-1:self.query_end]
		teRecord.id = id
		teRecord.description = id
		return teRecord

class RMOutParser():
	def __init__(self, inRMout):
		self.inRMout = inRMout
	def __iter__(self):
		return self._parse()
	def _parse(self):
		for line in open(self.inRMout):
			if not re.compile(r'.*\d').match(line): # title
				continue
			yield RMOutRecord(line)
	def get_seqs(self, genome, fout=sys.stdout):
		from Bio import SeqIO
		d_seqs = SeqIO.to_dict(SeqIO.parse(genome, 'fasta'))
		for rc in self._parse():
			seqRecord = d_seqs[rc.query_id]
			SeqIO.write(rc.get_seq(seqRecord), fout, 'fasta')

def mask(inSeq, inRM, outSeq, exclude=None, suffix=None, simpleSoft=True, excludeSoft=True, soft=False,):
	if suffix is None:
		suffix = os.path.splitext(inRM)[-1].lstrip('.')
	if suffix in {'bed'}:
		kw = dict(chr_col=0, start_col=1, end_col=2, based=0)
	elif suffix in {'gff3', 'gff', 'gff2', 'gtf'}:
		kw = dict(chr_col=0, start_col=3, end_col=4, based=1)
	elif suffix in {'out'}:
		kw = dict(chr_col=4, start_col=5, end_col=6, based=1)
	else:
		raise ValueError('un-recagnized suffix {}'.format(suffix))
	d_exclude = {}
	if exclude is not None:
		d_exclude = blast2dict(exclude, based=1)
	d_simple = {}
	if simpleSoft and not soft:  # only support RM.out
		d_simple = out2dict(inRM, based=1)
		#d_exclude.update(d_simple)
		update(d_exclude, d_simple)
		if excludeSoft:
			#d_simple.update(d_exclude)
			update(d_simple, d_exclude)
	d_mask = region2dict(inRM, **kw)
	all_n = 0
	for rc in SeqIO.parse(inSeq, 'fasta'):
		if rc.id not in d_mask:
			SeqIO.write(rc, outSeq, 'fasta')
			continue
		regions = list(set(d_mask[rc.id]) - set(d_exclude.get(rc.id, [])))
		seq = list(rc.seq)
		for start, end in sorted(regions):
			if soft:
				seq[start:end] = map(lower, seq[start:end])
			else:
				seq[start:end] = ['N'] * (end-start)
		regions2 = d_simple.get(rc.id, [])
		for start, end in sorted(regions2):
			seq[start:end] = map(lower, seq[start:end])
		print >> sys.stderr, 'total {}, exclude {}, hard {}, soft {}.'.format(len(set(d_mask[rc.id])), len(set(d_exclude.get(rc.id, []))), len(regions), len(regions2))
		d_counter = Counter(list(seq))
		num_N = d_counter['n'] + d_counter['N']
		seq = ''.join(seq)
		seq = Seq(seq)
		if not len(seq) == len(rc.seq):
			print >> sys.stderr, '\tERROR', rc.id, len(seq), len(rc.seq), len(regions)
		print >> sys.stderr, rc.id, len(seq), num_N, 100.0*(len(seq)-num_N)/len(seq)
		if len(seq) == num_N:
			all_n += 1
		assert len(seq) == len(rc.seq)
		rc.seq = seq
		SeqIO.write(rc, outSeq, 'fasta')
	print >> sys.stderr, all_n, 'all N'
def update(dt, dq):
	for key, value in dq.iteritems():
		dt[key] = set(dt.get(key, [])) | set(value)
	return dt
def lower(base):
	return base.lower()

def out2dict(inRMout, based=1, exclude={'ClassI', 'ClassII'}):
	d_region = {}
	for rc in RMOutParser(inRMout):
		if rc.Class in exclude:
			continue
		CHR, START, END = rc.query_id, rc.query_begin, rc.query_end
		START = START - based
		try: d_region[CHR] += [(START, END)]
		except KeyError: d_region[CHR] = [(START, END)]
	return d_region

def blast2dict(inBlast, based=1):
	d_region = {}
	for line in open(inBlast):
		try: 
			CHR, START, END = re.compile(r'(\S+):(\d+)[^d]+(\d+)').match(line).groups()
		except TypeError: pass
		START, END = int(START)-based, int(END)
		try: d_region[CHR] += [(START, END)]
		except KeyError: d_region[CHR] = [(START, END)]
	return d_region

def region2dict(inBed, chr_col=0, start_col=1, end_col=2, based=0):
	d_bed = {}
	for line in open(inBed):
		temp = line.strip().split()
		try:
			CHR, START, END = temp[chr_col], temp[start_col], temp[end_col]
			START = int(START) - based
			END = int(END)
		except IndexError: continue
		except ValueError: continue
		try: d_bed[CHR] += [(START, END)]
		except KeyError: d_bed[CHR] = [(START, END)]
	return d_bed

def RMOut2mnd(inRMout, outMnd, binsize0=300, scale=1, mapq=20, minlen=1000):
	from itertools import combinations
	d_family = {}
	for rc in RMOutParser(inRMout):
		if rc.query_match_length < minlen:
			continue
		try: d_family[rc.repeat_family] += [rc]
		except KeyError: d_family[rc.repeat_family] = [rc]
	i = 0
	for repeat_family, records in d_family.iteritems():
		for rc1, rc2 in combinations(records, 2):
			binsize = int(binsize0 * scale)
			bins = max(1, (rc1.query_match_length+rc2.query_match_length) / 2 / binsize)
			interval1 = (rc1.query_end - rc1.query_begin) / bins
			interval2 = (rc2.query_end - rc2.query_begin) / bins
			bins1 = range(rc1.query_begin, rc1.query_end, interval1)
			bins2 = range(rc2.query_begin, rc2.query_end, interval2)
			if rc1.strand == '-':
				bins1.reverse()
			if rc2.strand == '-':
				bins2.reverse()
			str1, str2 = 0,0
			chr1, chr2 = rc1.query_id, rc2.query_id
			frag1, frag2 = 0, 1
			mapq1, mapq2 = mapq, mapq
			cigar1,sequence1,cigar2,sequence2,readname1,readname2 = ['-'] * 6
			for bin1, bin2 in zip(bins1,bins2):
				pos1, pos2 = bin1, bin2
				line = [str1, chr1, pos1, frag1, str2, chr2, pos2, frag2, mapq1,cigar1,sequence1,mapq2,cigar2,sequence2,readname1,readname2]
				line = map(str, line)
				print >> outMnd, ' '.join(line)
				i += 1
	print >>sys.stderr, '{} links'.format(i)

def main():
	subcmd = sys.argv[1]
	if subcmd == 'out2mnd':
		inRMout = sys.argv[2]
		outMnd = sys.stdout
		RMOut2mnd(inRMout, outMnd)
	elif subcmd == 'out2seqs':
		inRMout = sys.argv[2]
		genome = sys.argv[3]
		RMOutParser(inRMout).get_seqs(genome)
	elif subcmd == 'mask':
		genome = sys.argv[2]
		inRMout = sys.argv[3]
		outSeq = sys.stdout
		try: exclude = sys.argv[4]
		except IndexError: exclude = None
		mask(genome, inRMout, outSeq, exclude=exclude)
	else:
		raise ValueError('Unknown command: {}'.format(subcmd))

if __name__ =='__main__':
	main()

