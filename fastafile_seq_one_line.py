#!/usr/bin/python

#-------------------------------------------------------------------------------#
#																				#
# Script write by: Aleix Arnau Soler (aleix.arnau1990@gmail.com) in the RBGE. 	#
#																				#
#	25/06/2014																	#
#																				#
# This script will convert your fasta file where each sequence is written in 	#
# different lines to a fasta file where sequences are in just one line in upper #
# letters. 																		#
#-------------------------------------------------------------------------------#

# Opening original fasta file

fasta_filename=raw_input('Write the fasta_file name: ')
transcriptome = open(fasta_filename)
transcriptome = transcriptome.readlines()

# opening output file
output_filename=raw_input('Write the output_file name: ')
output_transcriptome=open(output_filename,'w')

# main loop
count=0
for line in transcriptome:
	if count>0:
		if line.startswith('>'):
			output_transcriptome.write('\n'+line)
		else:
			sequence=line.upper()
			sequence=sequence.rstrip('\n')
			output_transcriptome.write(sequence)
	else:
		output_transcriptome.write(line)
	count+=1

output_transcriptome.close()
