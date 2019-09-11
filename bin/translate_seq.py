#!/bin/env python
import sys
from Bio import SeqIO
def six_frame_translate(inFa, fout=sys.stdout):
	for rc in SeqIO.parse(inFa, 'fasta'):
		for seq, suffix0 in zip([rc.seq, rc.seq.reverse_complement()], ['aa', 'rev_aa']):
			for frame in range(0,3):
				cds_seq = seq[frame:]
				aa_seq = translate_seq(cds_seq)
				suffix = '|{}{}'.format(suffix0, frame+1)
				print >> fout, '>{}{}\n{}'.format(rc.id, suffix, aa_seq)
			
def translate_seq(inSeq):
	aa = inSeq.translate()
	return aa

def main(inFa, outSeq=sys.stdout):
	for rc in SeqIO.parse(inFa, 'fasta'):
		print >> outSeq, '>{}\n{}'.format(rc.id, translate_seq(rc.seq))

if __name__ == '__main__':
	import sys
	inFa = sys.argv[1]
	main(inFa)
