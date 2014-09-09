#!/usr/bin/python

import time
import subprocess

#-------------------------------------------------------------------------------#
#																				#
# Script write by: Aleix Arnau Soler (aleix.arnau1990@gmail.com) in the RBGE. 	#
#																				#
#	07/05/2014																	#
#																				#
# This script has to be run in the same folder where there are the input files 	#
#																				#
#-------------------------------------------------------------------------------#

######################
#	 FUNCTIONS: 	 #
######################

# This function returs a dictionary with each read and their sequence from the fasta file.
def reads_with_seq(fasta_file):
	reads_file = open(fasta_file)
	reads=reads_file.readlines()
	reads_dictionary={}
	for line in reads:
		if line.startswith('>'):
			read_seq=''
			read_name=line.rstrip('\n')
			read_name=read_name.lstrip('>')
		else:
			read_seq = read_seq + str(line)
			read_seq = read_seq.rstrip()
			read_seq = read_seq.upper()
		reads_dictionary[read_name]=read_seq
	reads_file.close()
	return reads_dictionary

# This function provide a list with headers from all reads.
def list_all_reads(filename):
	filename = open(filename)
	file=filename.readlines()
	list_headers=[]
	for line in file:
		if line.startswith('>'):
			header=line.lstrip('>')
			header=header.rstrip('\n')
			list_headers.append(header)
	filename.close()
	return list_headers

# This functions provide a list with headers from all reads in clusters (NO singlets)
def list_reads_no_singlets(filename):
	filename=open(filename)
	file=filename.readlines()
	list_headers=[]
	for header in file:
		header=header.rstrip('\n')
		list_headers.append(header)
	filename.close()
	return list_headers

# This function return fasta format
def format_fasta(name,sequence):
	fasta_string='>{header}\n{sequence}\n'.format(header=name,sequence=sequence)
	return fasta_string

# This function create a fasta file with the singlet reads
def create_fasta_file(outputfilename,list_reads,dictionary):
	outputfile=(outputfilename)
	output_file=open(outputfile,'wt')
	for read in list_reads:
		output_file.write(format_fasta(read,dictionary[read]))

# This function returs a dictionary with all headers of singlet reads with their read ID.
def headers_singlets_index(index_file,list_singlets):
	headers_file = open(index_file)
	headers=headers_file.readlines()
	headers_dictionary={}
	for line in headers:	
		headers=line.split('\t')
		header=headers[0]
		index=headers [1].rstrip('\n')
		if index in list_singlets:
			headers_dictionary[index]=header
	headers_file.close()
	return headers_dictionary

# This function provide a list with the  singlet reads (all reads - non-singlet reads = singlet reads)
def list_singlets(list_all_reads,list_reads_no_singlets):
	list_singlets=[]
	for read in list_all_reads:
		if read not in list_reads_no_singlets:
			list_singlets.append(read)
	return list_singlets

# This function create a fasta file with the singlet reads with the original header
def create_fasta_file_original_headers(outputfilename,list_reads,dictionary_reads,dictionary_headers):
	outputfile=(outputfilename)
	output_file=open(outputfile,'wt')
	for read in list_reads:
		output_file.write(format_fasta(dictionary_headers[read],dictionary_reads[read]))

def translate_seconds(seconds):
	seg_ingresados = seconds
	for tmp in range (1,seg_ingresados+1):
		#print "Segundos..." + str(tmp)
		min = tmp/60
		sec= tmp%60
		hour = 0
		if min > 60:
			#vamos por las horas
			hour =min/60
			min= min%60
	if hour == 0:
		return str(min) + " minuts and "+str(sec)+" seconds"
	else:
		return str(hour)+" hours, "+str(min)+" minuts and "+str(sec)+"seconds"
	time.sleep(1)

###############################		END OF FUNCTIONS		#######################################################



print('\n#####################################################################\n')
subprocess.call(['ls'])
print('\n#####################################################################')

reads_fastafile=raw_input('\nWRITE THE NAME OF THE FASTA FILE WITH ALL READS: ')
no_singlet_reads=raw_input('\nWRITE THE NAME OF THE TEXT FILE WITH THE LIST OF NON-SINGLET READS: ')
headers_file=raw_input('\nWRITE THE NAME OF THE FILE WITH ALL THE HEADERS AND THEIR IDs: ')
ID=raw_input('\nWRITE AN IDENTIFIER FOR YOUR OUTPUTS: ')

print('\n#####################################################################')

# Tiempo de inicio de ejecucion.
inicio = time.time()

# Dictionary of all reads and their sequences from the fasta file
print('\nCreating Dict_reads_sequneces.....')
dict_reads_sequences=reads_with_seq(reads_fastafile)

# List with all read headers:
print('\nCreating list_all_reads.....')
list_all_reads = list_all_reads(reads_fastafile)

# List with non-singlets read headers.
print('\nCreating list_reads_no_singlets.....')
list_reads_no_singlets = list_reads_no_singlets(no_singlet_reads)

# List with all singlet reads.
print('\nCreating list_singlets.....')
list_singlets=list_singlets(list_all_reads,list_reads_no_singlets)
					# list_singlets=[read for read in list_all_reads if read not in list_reads_no_singlets]

# Dictionary with all index and their equivalent header from singlet reads.
print('\nCreating dict_singlets_headers.....')
dict_singlets_headers=headers_singlets_index(headers_file,list_singlets)

# Create a fasta file of singlet reads with original headers
print('\nCreating the fasta file with original headers.....')
create_fasta_file_original_headers(ID+'_singlets.fasta',list_singlets,dict_reads_sequences,dict_singlets_headers)

# To create a fasta file with the index headers.
		# Create a fasta file of singlet reads
		# print('\nCreating the fasta file.....\n')
		# create_fasta_file(ID+'_singlets.fasta',list_singlets,dict_reads_sequences)

# Tiempo de fin de ejecucion.
fin = time.time()

# Tiempo de ejecucion.
tiempo_total = fin - inicio
time_required=translate_seconds(int(tiempo_total))
print('\nThis script has required '+time_required+' to finish.\n')

print('There are '+str(len(list_singlets))+' singlet reads.\n')


