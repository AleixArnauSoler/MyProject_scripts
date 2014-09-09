#!/usr/bin/python

import sys
import subprocess

print('\n#####################################################################\n')
print('This program is used for joinnig files from the folder where you are.\n')
print('(follow the instructions)\n')
print ('THAT\'s YOUR FOLDER CONTENTS:\n')
subprocess.call(['ls'])
print('\n#####################################################################')

#Opening files and adding them in a single list.
file=[]
n_files=input('\nNUMBER OF FILES TO JOIN: ')
for x in range (n_files):
	filename=raw_input('\nFILENAME TO ADD: ')
	fp = open( filename, 'r' )
	file1=fp.readlines()
	file1=file1+list('\n')
	file=file+file1

# Writing the list to the outputfile.
outputname=raw_input('\nWRITE THE OUTPUT NAMEFILE: ')
output_file=open(outputname,'w')
for line in file:
	output_file.write(line)

print('\n#####################################################################')
print('\n!!!YOUR OUTPUT FILES HAVE BEEN CREATED SUCCESSFULY!!!\n')
