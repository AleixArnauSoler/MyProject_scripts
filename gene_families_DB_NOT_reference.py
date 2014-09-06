#!/usr/bin/python

import subprocess
import re
import collections
import time
import datetime

#################################################################################
#																				#
# Script write by: Aleix Arnau Soler (aleix.arnau1990@gmail.com) in the RBGE. 	#
#																				#
#									27/05/2014									#
#																				#
#################################################################################

####################################################################################################################
# 							>>>		FUNCTIONS:		<<<
####################################################################################################################

# This fuction provide the reverse complement of a sequence.
def rc(seq): 
	basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', '-': '-'} 
	useq = seq.upper()
	sort = str(type(useq))
	letters = list(useq) 
	completters = [basecomplement[base] for base in letters] 
	sort = str(type(completters))
	rc_list =  completters[::-1]
	sort = str(type(rc_list))
	rc_seq = "".join(rc_list)
	sort = str(type(rc_list))
	return rc_seq

# This function provide longest ORF from a sequence.
def codons(seq):
	try:
		ORF_list = []
		codon_list = []
		for i in range(len(seq)):
			codon = [seq[j:j+3] for j in range(i, len(seq), 3)]
			current_seq	 = []
			for i in codon:
				if i != 'TAG' and i != 'TAA' and i != 'TGA':
					current_seq.append(i)
				else:
					joined_seq = ''.join(current_seq)
					ORF_list.append(joined_seq)
					del(current_seq[:])
					break
		return(max(ORF_list, key = len))
	except ValueError:
		print('I dont think there are any sequences here')

# This function returs a list wit all isotigs from a fasta file.
def list_isotig(fasta_file):
	isotigs_file = open(fasta_file)
	isotigs=isotigs_file.readlines()
	list_isotigs=[]
	for line in isotigs:
		if line.startswith('>'):
			isotig_seq=''
			isotig_name=line.rstrip('\n')
			isotig_name=isotig_name.lstrip('>')
			list_isotigs.append(isotig_name)
	isotigs_file.close()
	return list_isotigs

# This function returs a dictionary with each isotig and their sequence from the fasta file.
def isotigs_with_seq(fasta_file):
	isotigs_file = open(fasta_file)
	isotigs=isotigs_file.readlines()
	isotig_dictionary={}
	for line in isotigs:
		if line.startswith('>'):
			isotig_seq=''
			isotig_name=line.rstrip('\n')
			isotig_name=isotig_name.lstrip('>')
		else:
			isotig_seq = isotig_seq + str(line)
			isotig_seq = isotig_seq.rstrip()
			isotig_seq = isotig_seq.upper()
		isotig_dictionary[isotig_name]=isotig_seq
	isotigs_file.close()
	return isotig_dictionary

# This function return fasta format
def format_fasta(name,sequence):
	fasta_string='>{header}\n{sequence}\n'.format(header=name,sequence=sequence)
	return fasta_string

# This fuction returns a dictionary with each sequence and their header from a fastafile.
def dictionary_seq(filename):
	dictionary = {}
	for line in filename:
		if line.startswith(">"):
			C_seq = ''
			C_split_line = line.split(' ')
			C_name = C_split_line[0]
			C_name = C_name.rstrip()
			C_name = C_name.lstrip('>')
		else:
			C_seq = C_seq + str(line)
			C_seq = C_seq.rstrip()
			C_seq = C_seq.upper()
		dictionary[C_name] = C_seq
	return dictionary

# This function translate the seconds required to run the script in hours, minuts and seconds.
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

####################################################################################################################
# 							>>>		END OF FUNCTIONS:		<<<
####################################################################################################################							

################################		PROGRAM			############################################################

print('\n#####################################################################')
print('\n Hi! This script will create a DB with clusters & gene family information about the transcriptomes you introduce:\n')
print('\n It will need that you introduce the next information to be run:\n\n')
print('   * The number of species you want to analize.')
print('   * An identifier for each specie.\n   * The path of each specie\'s transcriptome.')
print('   * A name to identify the output file gene families DB.\n')
print(' REMEMBER:\n If you misspell some of the inputs required you will must start once again.')
print(' If some file is in the working directory you only need to write the name of the file (is not necessary all path).')
print('\n######################################################################')

																		


# Temporal output file.
all_transcriptomes=open('all_transcriptomes.fasta','w')

# Introduce the number of species to analize.
number_species=(raw_input('\nHow many species do you want to introduce? '))

# Opening species's data.
dictionary_seq_dictionaries={}
identifiers_list=[]

# Introduce paths for each specie.
for i in range(int(number_species)):
	# Identifier for the specie.
	identifier=raw_input('\nWrite the identifier for the specie '+str(i+1)+' (the same used to identify the transcriptome): ')
	# Creat a list with all identifiers.
	identifiers_list.append(identifier)
	
	# Opening transcriptome
	path=raw_input('\nWrite the path of the fasta_file of the specie '+identifier+': ')
	print (' opening the transcriptome of specie '+str(i+1)+'...')

	# List of isotigs.	
	list_isotigs=list_isotig(path)
	# Dictionary with each isotig and its sequence.
	dictionari_transcriptome=isotigs_with_seq(path)

	# Dictionary with all sequence dictionaries.
	filename = open(path)
	dictionary_seq_dictionaries['{0}'.format(identifier)]=dictionary_seq(filename)
	filename.close()

	# Editing output fasta.
	for isotig in list_isotigs:
		all_transcriptomes.write(format_fasta(identifier+'_'+isotig,dictionari_transcriptome[isotig]))
	
	
# Create a tuple from the list of identifiers (for further analysis)
identifiers_tuple=tuple(identifiers_list)

# Data and hour.
now=datetime.datetime.now()
Datetime=now.strftime("Date: %Y-%m-%d\nHour: %H:%M")

# Create the output file.
final_output = raw_input('\nWrite the output filename: ')
final_output = open(final_output, 'w')

# Create the singlet isotigs output_file.
singlets_output = open('singlets_'+str(identifiers_list), 'w')

# Create a file for possible errors.
errors_output= open ('errors_'+str(identifiers_tuple), 'w')

# Time of ejecution.
time_beginning = time.time()


####################################################	COMAND LINE:	- BLAST & MCL - 	##################################################################

# Create a BlastDB with all transcriptomes.
subprocess.call(['makeblastdb -in all_transcriptomes.fasta -dbtype nucl'],shell=True)

# Blast all transcriptomes against the DB.
subprocess.call(['blastn -query all_transcriptomes.fasta -db all_transcriptomes.fasta -out blast_gene_families.abc -outfmt "6 qseqid sseqid evalue"'],shell=True)


# MCL
subprocess.call(['mcxload -abc blast_gene_families.abc --stream-mirror --stream-neg-log10 -stream-tf "ceil(200)" -o blast_gene_families.mci -write-tab blast_gene_families.tab'],shell=True)
			# CHECK PARAMETERS

subprocess.call(['mcl blast_gene_families.mci -I 3 -use-tab blast_gene_families.tab'],shell=True)

##############################################################################################################################################################



# Open the dump file..
mcl_file = ('out.blast_gene_families.mci.I30')
print ('\n opening '+mcl_file +' file...\n')
f = open(mcl_file)
dump = f.readlines()
print ('Analazing... wait...\n')



########################################	ANALYSIS	################################################################################

# Counting clusters and clusters with copy numbers.
cluster_sum=0
cluster_with_copy=0
for line in dump:
	cluster_sum+=1
	if any(re.match(item, line) for item in identifiers_list):
		item = line.rstrip('\n')
		clusters = item.split('\t')
		try:
			clusters[1]
			cluster_with_copy+=1
		except IndexError:
			pass

# Editing final output file.
final_output.write('This is a gene families DB from '+str(int(number_species)+1)+' species.\n\n')
final_output.write(Datetime)
final_output.write('\n\nSpecie IDs:')
for ID in identifiers_list:
	final_output.write(' '+ID+'\n')					
final_output.write('\nTotal number of CLUSTERS analysed: '+str(cluster_sum))
final_output.write('\nTotal number of CLUSTERS with COPY NUMBER (counting singlets): '+str(cluster_with_copy)+'\n\n')

# Headers.
for i in range(int(number_species)):
	final_output.write('Number of isotigs from specie ' + str(i+1) + '.\t')
final_output.write('Minimum Kn/ks value.\tMaximum Kn/ks value.\tRange Kn/ks value.\t')
final_output.write('Number consensus bases\tNumber total bases aligned\tIsotig names from the cluster.\n')

# Main Loop
dictionary_copy_dictionaries={}

cluster_sum2=0
for line in dump:
	
	consensus_bases = 0
	total_bases = 0
	knks_list=[]
	knks_list_no1 = []
	whole_name = []

	dictionary_copy_dictionaries={}

	cluster_sum2+=1
	if any(re.match(item, line) for item in identifiers_list):
		file = open('dump.orfs', 'w')
		item = line.rstrip('\n')
		clusters = item.split('\t')														

		try:
			clusters[1]

			consensus_bases = 0
			total_bases = 0
			knks_list=[]
			knks_list_no1 = []
			whole_name = []

			dictionary_copy_dictionaries={}
			copy_total=0


			for ID in identifiers_list:
				copy=0
				for isotig in clusters:
					if isotig.startswith(ID):
						isotig=isotig.lstrip(ID).lstrip('_')					# CONCORDEN LS IDs?
						file.write('>' + ID + isotig + '\n' + dictionary_seq_dictionaries[ID][isotig] + '\n')	
						# List with all isotigs in the cluster
						whole_name.append(ID + isotig)
						copy+=1
						copy_total+=1
						dictionary_copy_dictionaries['{0}'.format(ID)]= copy
				
			file.close()
				
			print('\nCLUSTER '+str(cluster_sum2))
			print('total copy number is ' + str(copy_total))
			print('\nthe sequence names are: ' + str(clusters))
			print('\n')		
							
			longest = open('dump.orfs')
			outfile = open('dump.longest.orfs', 'w')
			# Dictionary with the sequences from dump.orfs
			dictionary_longests=dictionary_seq(longest)		
			longest.close()
				
			if int(len(dictionary_longests)) == 0:
				break

			for isotig in dictionary_longests:		
				seq = dictionary_longests[isotig]			

				if int(len(list(seq)))>200:			
					FandR = []			
					Fwd = codons(seq)			
					Rev = codons(rc(seq))			
					FandR.append(Fwd)
					FandR.append(Rev)

					try:				
						outfile.write('>' + isotig + '\n' + max(FandR, key = len) +'\n')
					except TypeError:
						print('no sequences')
						
			outfile.close()


			####################################################	COMAND LINE 	####################################################################
				
			try:
				subprocess.call(["python dna2pep-1.1/dna2pep.py dump.longest.orfs --fasta dump.orfs.pep"], shell=True)
				# in clusters which have cucumber in them - it doesnt recognise the cucumber sequences?
							
				subprocess.call(["linsi --thread 8 --quiet dump.orfs.pep > dump.orfs.pep.align"], shell=True)
					
				subprocess.call(["python RevTrans-1.4/revtrans.py dump.longest.orfs dump.orfs.pep.align dump.orfs.revtrans"], shell=True)
						

				subprocess.call(["linsi --thread 8 dump.longest.orfs > dump.nuc1"], shell=True)	
					
			##########################################################################################################################################


				f = open('dump.orfs.revtrans')
				file = open('dump.nuc', 'w')
				aln = f.read()
				aln = aln.replace('isotig', '')
				for ID in identifiers_list:
					aln = aln.replace(ID, ID[0])
				file.write(aln)
				f.close()
				file.close()	
					
					
				####################################################	R ANALYSIS 	####################################################################
					
				# this bit is to count number of bases that are aligned, using consensus with 0.6 threshold
					
				subprocess.call(["Rscript consensus.R"], shell=True)
						
				f = open('alignment_consensus')
				consensus = f.readlines()
				
				for i in consensus:
					total_bases+=1
					i=i.rstrip('\n')
					if i == 'NA':
						pass
					elif i == '-':
						pass
					else:
						consensus_bases+=1

				##########################################################################################################################################
					

				####################################################	COMAND LINE 	##################################################################

				subprocess.call(["/usr/genome/Pipeline/EMBOSS-6.6.0/emboss/seqret -sequence dump.nuc -outseq dump.revtrans.phylip -osformat2 phylip"], shell=True)
				
				subprocess.call(["/usr/genome/Pipeline/standard-RAxML-master/raxmlHPC-SSE3 -p 100 -s dump.revtrans.phylip -n dump.revtrans.raxml -m PROTCATWAG"], shell=True)

				subprocess.call(['codeml'], shell=True) 

				# Maybe better create some command prompts in /usr/bin/ with these executable file (ex:codeml.exe)
			
			except IndexError:
				knks_list_no1 = ['N/A']					
				#########################################################################################################################################
			try:

				f = open('paml_out')
				paml_file = f.read()
					
				x = paml_file.split('comparison.)\n\n')
				q = (x[1])		
					
				q = re.sub(r'\([^)]*\)', '', x[1])
					
				q = q.replace('\n', ' ')
				q = q.replace('  ', ' ')
				q = q.replace('(', '')
				q = q.replace(')', '')
				q = q.replace('NOTE: -1 means that NG86 is inapplicable.', '')
				
				q = q.split()
				knks_list = [float(item) for item in q if not item.startswith(identifiers_tuple)]

				
				for i in knks_list:
					if i != '-1.0000':
						knks_list_no1.append(i)
			except IOError:
				knks_list_no1 = ['N/A']
		
		except IndexError:
			singlets_output.write(line)

	# else:
	# 	# CLUSTERS with only resference transcriptome isotigs.
	# 	print('CLUSTER '+str(cluster_sum2))
	# 	print('NO COPY NUMBERS')
	# 	item = line.rstrip('\n')
	# 	clusters = item.split('\t')	
		
	# 	try:
	# 		clusters[1]

	# 		consensus_bases = ('N/A')
	# 		total_bases = ('N/A')
	# 		knks_list_no1 = []
	# 		whole_name=[]

	# 		dictionary_copy_dictionaries={}
			
																												
	# 	except IndexError:
	# 		singlets_output.write(line)
	# 	except KeyError:
	# 		errors_output.write(line)

	# EDITING OUTPUT


	
	try:
		clusters[1]


		# Number of isotigs from the other species.
		for i in identifiers_list:
			if i in dictionary_copy_dictionaries:
				final_output.write(str(dictionary_copy_dictionaries[i]) + '\t')
			else:
				final_output.write(str(0) + '\t')		

	
		try:
			# Minimum, maximum and avarage Kn/Ks values.
			final_output.write(str(min(knks_list_no1)))
			final_output.write('\t')
		except ValueError:
			final_output.write('N/A' + '\t')
		except TypeError:
			final_output.write('N/A' + '\t')
		except NameError:
			final_output.write('N/A' + '\t')		
	
		try:
			final_output.write(str(max(knks_list_no1)))
			final_output.write('\t')
		except ValueError:
			final_output.write('N/A' + '\t')
		except TypeError:
			final_output.write('N/A' + '\t')
		except NameError:
			final_output.write('N/A' + '\t')	
	
		try:
			final_output.write(str(max(knks_list_no1) - min(knks_list_no1)))
			final_output.write('\t')		
		except ValueError:
			final_output.write('N/A' + '\t')
		except TypeError:
			final_output.write('N/A' + '\t')
		except NameError:
			final_output.write('N/A' + '\t')	
					
		# Number of consensus bases on the alignment (same base in different isotigs).
		try:
			final_output.write(str(consensus_bases) + '\t')
		except ValueError:
			final_output.write(('N/A') + '\t')
		except NameError:
			final_output.write('N/A' + '\t')
	
		# Number of total bases aligned.
		try:
			final_output.write(str(total_bases) + '\t')
		except ValueError:
			final_output.write(('N/A') + '\t')
		except NameError:
			final_output.write('N/A' + '\t')
	
		# Name of the isotigs from the species analised.
		for isotig in whole_name:
			final_output.write(isotig + '\t')
		final_output.write('\n')

	except IndexError:
		pass

	try:
		subprocess.call(["rm *dump.revtrans.raxml*"], shell=True)
	except :
		pass
	try:
		subprocess.call(["rm paml_out"], shell=True)
	except :
		pass

#########################################################################################################################################
# 										END OF ANALYSIS
#########################################################################################################################################

# Time finished.
time_end = time.time()

# Time required.
total_time = time_end - time_beginning
time_required=translate_seconds(int(total_time))

print('#####################################################################')
print('\n!!!IT HAS FINISHED SUCCESSFULY!!!')
print('\n * Total number of CLUSTERS: '+str(cluster_sum))
print(' * Total number of CLUSTERS with COPY NUMBER: '+str(cluster_with_copy))
print('\n This script has required '+time_required+' to finish.\n')
