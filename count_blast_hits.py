#!/usr/bin/python

# This script count number of hits from a blast file.

filename=raw_input('WRITE the name of the blast file: ')
input_filename=open(filename)
input_file=input_file.readlines()

hits_list=[]
for line in input_file:
	line=line.split()
	hits=line[0]
	hits_list.append(hits)

final_list=set(hits_list)
length_list=len(final_list)

print('\n')
print (length_list)
print('\n')
