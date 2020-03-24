#!/bin/env python
import sys
from Bio import SeqIO
from Bio.Data import CodonTable
def six_frame_translate(inFa, fout=sys.stdout, seqfmt='fasta', transl_table=1):
	d_length = {}
	for rc in SeqIO.parse(inFa, seqfmt):
		for seq, suffix0 in zip([rc.seq, rc.seq.reverse_complement()], ['aa', 'rev_aa']):
			for frame in range(0,3):
				nucl_seq = seq[frame:]
				try: aa_seq = translate_seq(nucl_seq, table=transl_table)
				except CodonTable.TranslationError: continue   # Codon 'XGA' is invalid
				suffix = '|{}{}'.format(suffix0, frame+1)
				print >> fout, '>{}{}\n{}'.format(rc.id, suffix, aa_seq)
		d_length[rc.id] = len(rc.seq)
	return d_length
			
def translate_seq(inSeq, **kargs):
	aa = inSeq.translate(**kargs)
	return aa
def translate_cds(inSeq, transl_table=1, **kargs):
	for key in kargs.keys():
		if not key in {'to_stop', 'stop_symbol', 'gap'}:
			del kargs[key]
	try:
		aa = translate_seq(inSeq, cds=True, table=transl_table, **kargs)
	except CodonTable.TranslationError as e:
		aa = translate_seq(inSeq, table=transl_table, **kargs)
	return aa

def main(inFa, outSeq=sys.stdout):
	for rc in SeqIO.parse(inFa, 'fasta'):
		print >> outSeq, '>{}\n{}'.format(rc.id, translate_seq(rc.seq))

if __name__ == '__main__':
	import sys
	inFa = sys.argv[1]
	if inFa == 'six_frame_translate':
		inFa = sys.argv[2]
		six_frame_translate(inFa)
	else:
		main(inFa)
