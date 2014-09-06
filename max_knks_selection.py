#!/usr/bin/python

filename=raw_input("\nWrite the cluster DB filename:")
DB_file=open(filename)
DB_file=DB_file.readlines()

knks_value=input("\nSelect a threshold for the max knks value:")

output_name=raw_input("\nWrite an output name:")
output=open(output_name,"w")

count=0

for line in DB_file:
	line=line.split("\t")

	try:
		value=line[5]					# Need to be adabted to correct column depending on the DB
		if float(value)>float(knks_value):
			print line
			for i in range(len(line)-1):
							
				output.write(line[i]+"\t")
			output.write("\n")

	except IndexError:
		pass
	except ValueError:
		pass
