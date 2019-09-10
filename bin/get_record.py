#!/usr/bin/python
#coding: utf-8
#2016-8-25, Wei,Y.L. print not get id, one id per line.

from Bio import SeqIO
import sys, getopt
from small_tools import open_file as open
#open=open_file

def usage():
	print 'usage: %s [options] -i RECORD -a LIST -o SUBRECORD' %sys.argv[0]
	print ''
	print '  -i RECORD         input a RECORD FILE'
	print '  -a LIST           input a LIST FILE, one RECORD ID per line'
	print '  -o SUBRECORD      output to SUBRECORD FILE'
	print '  -t RECORD TYPE    RECORD FILE TYPE [table|fasta|fastq|hmm][default: table]'
	print '  -g STR            [get|remove] RECORD [default: get]'
	print '  -k NUM            if a table RECORD, the column NUM of RECORD ID[default: 1]'
	print '  -f NUM            if a table RECORD, retain the row NUM as header [default: 1]'
	print '  -h                show this help message and exit'

opts, args = getopt.getopt(sys.argv[1:], 'hi:a:o:t:g:k:f:')
input_file = ''
output_file = ''
in_accnos = sys.stdin
type = 'table'
process = 'get'
col = 1
head = 1

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
	elif op == '-h':
		usage()
		sys.exit()

if type not in ['table','fasta','fastq', 'hmm']:
	print "Error: type must be one of ['table','fasta','fastq'], unexpected '%s'" % (type,)
	usage()
	sys.exit()
if process not in ['get','remove']:
	print "Error: process must be one of ['get','remove'], unexpected '%s'" % (process,)
	usage()
	sys.exit()

#def overwrite_output(filename):
#        import os
#        if os.path.exists(filename):
#                os.remove(filename)
#                os.mknod(filename)
#        else:
#                os.mknod(filename)

#overwrite_output(output_file)

def get_record(d_accnos, record_id):
	if d_accnos.has_key(record_id): #record_id in d_accnos.keys():
		return True
	else:
		return False

def remove_record(d_accnos, record_id):
	if d_accnos.has_key(record_id): #record_id not in d_accnos.keys():
		return False #True
	else:
		return True #False

d_accnos = {}
for line in in_accnos:
	if not line.strip():
		continue
	temp = line.strip('\n\r').split()[0]
	d_accnos[temp] = None

#def parse_table(input_file, d_accnos, process, col, head, type):
#	i = 0
lst_get = []
f = open(output_file, 'w')
if type == 'table':
	i = 0
#	f = open(output_file, 'a')
	for line in open(input_file,'r'):
		i += 1
		if i == head:
			f.write(line)
		else:
			temp = line.strip('\n').split('\t')
			record_id = temp[col-1]
			if process == 'get':
				if get_record(d_accnos, record_id):
					f.write(line)
					lst_get.append(record_id)
				else:
					continue
			elif process == 'remove':
				if remove_record(d_accnos, record_id):
					f.write(line)
				else:
					lst_get.append(record_id)
					continue
	#f.close()

elif type == 'fasta' or type == 'fastq':
#	f = open(output_file, 'a')
	for seq_record in SeqIO.parse(open(input_file),type):
		record_id = seq_record.id
		if process == 'get':
			if get_record(d_accnos, record_id):
				SeqIO.write(seq_record, f, type)
				lst_get.append(record_id)
			else:
				continue
		elif process == 'remove':
			if remove_record(d_accnos, record_id):
#				lst_get.append(record_id)
				SeqIO.write(seq_record, f, type)
			else:
				lst_get.append(record_id)
				continue
elif type == 'hmm':
	from HMMER import HMMParser
	for rc in HMMParser(open(input_file)):
		if process == 'get':
			if rc.NAME in d_accnos or rc.ACC in d_accnos:
				rc.write(f)
				lst_get += [rc.NAME, rc.ACC]
		elif process == 'remove':
			if not (rc.NAME in d_accnos or rc.ACC in d_accnos):
				rc.write(f)
				lst_get += [rc.NAME, rc.ACC]
f.close()

not_get = set(d_accnos.keys())-set(lst_get)
if not_get:
	for not_get_id in not_get:
		print not_get_id
