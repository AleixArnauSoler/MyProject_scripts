#!/usr/bin/python

import subprocess

#-------------------------------------------------------------------------------#
#																				#
# Script write by: Aleix Arnau Soler (aleix.arnau1990@gmail.com) in the RBGE. 	#
#																				#
#	23/06/2014																	#
#																				#
# This script has to be run in the same folder where there are the input files 	#
#																				#
#-------------------------------------------------------------------------------#

######################
#	 FUNCTIONS: 	 #
######################
# This function returs a dictionary with each singlet and their sequence.
def singlets_with_seq(fnafile):
	singlet_file = open(fnafile)
	singlet=singlet_file.readlines()
	singlet_dictionary={}
	for line in singlet:
		if line.startswith('>'):
			singlet_seq=''
			singlet_name=line.lstrip('>')
			singlet_name=singlet_name.rstrip('\n')
		else:
			singlet_seq = singlet_seq + str(line)
			singlet_seq = singlet_seq.rstrip()
			singlet_seq = singlet_seq.upper()
		singlet_dictionary[singlet_name]=singlet_seq
	singlet_file.close()
	return singlet_dictionary

# This function returs a dictionary with each transcript and their sequence.
def isotigs_with_seq(fnafile):
	isotigs_file = open(fnafile)
	isotigs=isotigs_file.readlines()
	isotig_dictionary={}
	for line in isotigs:
		if line.startswith('>'):
			isotig_seq=''
			isotig_split_line=line.split(' ')
			isotig_name=isotig_split_line[0]
			isotig_name=isotig_name.lstrip('>')
			isotig_name=isotig_name.rstrip('\n')
		else:
			isotig_seq = isotig_seq + str(line)
			isotig_seq = isotig_seq.rstrip()
			isotig_seq = isotig_seq.upper()
		isotig_dictionary[isotig_name]=isotig_seq
	isotigs_file.close()
	return isotig_dictionary



# This function returns a list with all singlets that are in the blast file (which match isotigs).
def singlets_match_transcripts(blastname):
	blast_file=open(blastname)
	hits = blast_file.readlines()
	singlets_list=[]
	for line in hits:
		line_split=line.split('\t')
		singlet=line_split[1]
		singlets_list.append(singlet)
	blast_file.close()
	return singlets_list





# This function returns a list with all transcripts that are in the blast file (which match singlets).
def transcripts_match_singlets(blastname):
	blast_file=open(blastname)
	hits = blast_file.readlines()
	isotig_list=[]
	for line in hits:
		line_split=line.split('\t')
		isotig=line_split[0]
		isotig_list.append(isotig)
	blast_file.close()
	return isotig_list

# This function returns a list with all transcripts from the transcriptome
def isotigs_name(fnafile):
	isotigs_file = open(fnafile)
	isotigs=isotigs_file.readlines()
	isotig_list=[]
	for line in isotigs:
		if line.startswith('>'):
			isotig_split_line=line.split(' ')
			isotig_name=isotig_split_line[0]
			isotig_name=isotig_name.lstrip('>')
			isotig_list.append(isotig_name)	
	isotigs_file.close()
	return isotig_list

# This function return isotigs that don't are in the list of match
def filter_match(list1,list2):
	isotigs=[]
	for isotig1 in list1:
		if isotig1 not in list2:
			isotigs.append(isotig1)
	return isotigs

# This function return fasta format
def format_fasta(name,sequence):
	fasta_string='>{header}\n{sequence}\n'.format(header=name,sequence=sequence)
	return fasta_string

# This function create a fasta file with the isotigs which there are in the list using the dictionary
def fasta(outputfilename,listisotigs,dictionary):
	outputfile=(outputfilename)
	output_file=open(outputfile,'wt')
	for isotig in listisotigs:
		output_file.write(format_fasta(isotig,dictionary[isotig]))




print('\n#####################################################################')
print('\nHi! This script will create:\n   * An output file with all transcripts that match to the singlets from unmapped reads\n   * An output file with all transcripts that match to the singlets from unmapped reads.\n')
print('You\'ll have to introduce the two input filenames and the two output filenames. So, look at the files you have in that folder to write correctly their names:\n')

subprocess.call(['ls'])

print('\n#####################################################################')

transcriptome_file =raw_input('\nWRITE THE TRANSCRIPTOME\'S NAMEFILE (the file that comes from Newbler): ')
singlets_file =raw_input('\nWRITE THE SINGLETS\'S NAMEFILE: ')
blast_file =raw_input('\nWRITE THE NAMEFILE OF THE BLAST FILE: ') 
output_singlets=raw_input('\nWRITE THE BANE YOU WANT TO GIVE TO THE FASTA FILE WITH THE SINGLETS WHICH MATCH ISOTIGS: ')
output_match=raw_input('\nWRITE THE NAME YOU WANT TO GIVE TO THE FILE WITH THE ISOTIGS THAT MATCH SINGLETS: ')
output_no_match=raw_input('\nWRITE THE NAME YOU WANT TO GIVE TO THE FILE WITH THE ISOTIGS THAT DON\'T MATCH SINGLETS: ')


#Dictionary of all singlets with their sequence
dictionary_singlets=singlets_with_seq(singlets_file)
# List of all singlets which match with blastn to the transcripts
list_singlets_match=singlets_match_transcripts(blast_file)


# Dictionary of all isotigs with their sequence
dictionary_isotigs =isotigs_with_seq(transcriptome_file)
# List of all isotigs from transcriptome
list_isotigs_name=isotigs_name(transcriptome_file)
# List of all transcripts which match with blastn to the singlets
list_isotigs_match=transcripts_match_singlets(blast_file)

list_transcripts_match=sorted(set(list_isotigs_match))
print(len(list_transcripts_match))


# List of all trancripts which don't match with blastn to the singlets
list_isotigs_not_match=filter_match(list_isotigs_name,list_transcripts_match)


# Create fasta of singlets that match transcripts
fasta(output_singlets,list_singlets_match,dictionary_singlets)
# Create fasta of isotigs that match singlets
fasta(output_match,list_transcripts_match,dictionary_isotigs)
# Create fasta of isotigs that don't match singlets
fasta(output_no_match,list_isotigs_not_match,dictionary_isotigs)

print('\n#####################################################################')
print('\n!!!YOUR OUTPUT FILES HAVE BEEN CREATED SUCCESSFULY!!!\n')
print('\nLOOK AT THEM:\n')

subprocess.call(['ls'])

