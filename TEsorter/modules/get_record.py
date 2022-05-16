#!/usr/bin/env python
#coding: utf-8
import sys, getopt
from io import IOBase
from Bio import SeqIO
from TEsorter.modules.small_tools import open_file as open

def usage():
	print('usage: %s [options] -i RECORD -a LIST -o SUBRECORD' %sys.argv[0])
	print('')
	print('  -i RECORD         input a RECORD FILE')
	print('  -a LIST           input a LIST FILE, one RECORD ID per line')
	print('  -o SUBRECORD      output to SUBRECORD FILE')
	print('  -t RECORD TYPE    RECORD FILE TYPE [table|fasta|fastq|hmm][default: table]')
	print('  -g STR            [get|remove] RECORD [default: get]')
	print('  -k NUM            if a table RECORD, the column NUM of RECORD ID[default: 1]')
	print('  -f NUM            if a table RECORD, retain the row NUM as header [default: 1]')
	print('  -s SEP            if a table RECORD, seperation [default: "\\t"]')
	print('  -d                remove duplicates when get table records')
	print('  -h                show this help message and exit')

def main():
	opts, args = getopt.getopt(sys.argv[1:], 'hi:a:o:t:g:k:f:s:x:d')
	input_file = ''
	output_file = ''
	in_accnos = sys.stdin
	type = 'table'
	process = 'get'
	col = 1
	head = 1
	accnos_sep = None
	sep = "\t"
	dedup = False
	for op, value in opts:
		if op == '-i':
			input_file = value
		elif op == '-o':
			output_file = value
		elif op == '-a':
			in_accnos = open(value)
		elif op == '-t':
			type = value
		elif op == '-g':
			process = value
		elif op == '-k':
			col = int(value)
		elif op == '-f':
			head = int(value)
		elif op == '-x':
			accnos_sep = value
		elif op == '-s':
			sep = value
		elif op == '-d':
			dedup = True
		elif op == '-h':
			usage()
			sys.exit()

	if type not in {'table','fasta','fastq', 'hmm', 'genbank'}:
		raise TypeError("type must be one of ['table','fasta','fastq'], unexpected '%s'" % (type,))
		usage()
		sys.exit()
	if process not in ['get','remove']:
		raise TypeError("process must be one of ['get','remove'], unexpected '%s'" % (process,))
		usage()
		sys.exit()
	get_records(input_file, output_file, in_accnos, type=type, process=process, col=col, head=head, accnos_sep=accnos_sep, sep=sep, dedup=dedup)

def _format_id(id):
	return id
def get_records(input_file, output_file, in_accnos,
				type='table', process='get', sep="\t", format_id=_format_id,
				col=1, head=1, accnos_sep=None, dedup=False):
	def get_record(d_accnos, record_id):
		if record_id in d_accnos:
			return True
		else:
			return False
	def remove_record(d_accnos, record_id):
		if record_id in d_accnos:
			return False
		else:
			return True
	if isinstance(in_accnos, IOBase):
		d_accnos = {line.strip().split(accnos_sep)[0] for line in in_accnos if line.strip()}
	else: # list
		d_accnos = set(in_accnos)

	lst_get = set([])
	f = open(output_file, 'w')
	if type == 'table':
		i = 0
		for line in open(input_file,'r'):
			i += 1
			if i == head:
				f.write(line)
			else:
				temp = line.strip().split(sep)
				record_id = format_id(temp[col-1])
				temp[col-1] = record_id
				if process == 'get':
					if get_record(d_accnos, record_id):
						if dedup and record_id in lst_get:
							continue
						f.write(line)
						lst_get.add(record_id)
					else:
						continue
				elif process == 'remove':
					if remove_record(d_accnos, record_id):
						f.write(line)
					else:
						lst_get.add(record_id)
						continue

	elif type in { 'fasta', 'fastq', 'genbank'}:
		for seq_record in SeqIO.parse(open(input_file),type):
			record_id = format_id(seq_record.id)
			seq_record.id = record_id
			if process == 'get':
				if get_record(d_accnos, record_id):
					SeqIO.write(seq_record, f, type)
					lst_get.add(record_id)
				else:
					continue
			elif process == 'remove':
				if remove_record(d_accnos, record_id):
					SeqIO.write(seq_record, f, type)
				else:
					lst_get.add(record_id)
					continue
	elif type == 'hmm':
		from HMMER import HMMParser
		for rc in HMMParser(open(input_file)):
			values = {getattr(rc, key, None) for key in ['NAME', 'ACC']}
			values = {v for v in values if v}
			if process == 'get':
				if values & d_accnos:
					rc.write(f)
					lst_get = lst_get | values
			elif process == 'remove':
				if not (values & d_accnos):
					rc.write(f)
					lst_get = lst_get | values
	f.close()

	not_get = d_accnos-set(lst_get)
	if not_get:
		for not_get_id in not_get:
			print(not_get_id)

if __name__ == '__main__':
	main()
