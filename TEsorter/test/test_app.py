#!/usr/bin/env python

import subprocess
import os
import sys

def fill_path(file):
	path =  os.path.realpath(__file__)
	tail = path.rsplit('/',1)
	return tail[0]+"/"+file

def main():
	#
	DATA="rice6.9.5.liban"

	#Create the command
	command = "TEsorter "+fill_path(DATA)
	print("Running the following command: "+command)

	#Execute the command
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
	stdout, stderr = process.communicate()
	print (stdout.decode('utf-8'))

if __name__ == '__main__':
	main()
