#!/usr/bin/env python


def getAttributes(text):
	attrs = {}
	fields = [x.strip() for x in text.split(';')]
	for field in fields:
		if field!="":
			values = [x.strip() for x in field.split('"')[0:-1]]
		if len(values) > 0:
			attrs[values[0].lower()] = values[1]
	return attrs

# inputFile = open("/home/ldelafuente/squanti/myNeuralData/refseq_ensembl_miRNA_smallRNAs_badAligned_filtered.genePred", "r")

# transcripts = []
# for i in inputFile:
# 	transcripts.append(i.split("\t")[0])

# inputFile.close()

# print(len(transcripts))

# ## reference 1

# inputFile = open("/home/ldelafuente/squanti/myNeuralData/Refseq_global_PATH1_ok.gtf", "r")
# output = open("/home/ldelafuente/squanti/myNeuralData/Refseq_global_PATH1_ok_filtered.gtf", "w")

# for line in inputFile:
# 	fields = line.split("\t")
# 	attrs = getAttributes(fields[8])
# 	pbid = attrs["transcript_id"]	
# 	if pbid in transcripts:
# 		output.write(line)	

# inputFile.close()
# output.close()

# ## reference 2

# inputFile = open("/home/ldelafuente/squanti/myNeuralData/ENSEMBL_global_PATH1_ok.gtf", "r")
# output = open("/home/ldelafuente/squanti/myNeuralData/ENSEMBL_global_PATH1_ok_filtered.gtf", "w")

# for line in inputFile:
# 	fields = line.split("\t")
# 	attrs = getAttributes(fields[8])
# 	pbid = attrs["transcript_id"]	
# 	if pbid in transcripts:
# 		output.write(line)	

# inputFile.close()
# output.close()


### alias 

inputFile = open("/home/ldelafuente/squanti/myNeuralData/ENSEMBL_global_PATH1_ok_filtered.gtf", "r")
aliasFile = open("/home/ldelafuente/squanti/ensembl_NCBI_genes.txt", "r")
output = open("/home/ldelafuente/squanti/myNeuralData/ENSEMBL_global_PATH1_ok_filtered_NCBIaliases.gtf", "w")
ensembl_ncbi = open("/home/ldelafuente/squanti/aliases_changed.txt", "w")
aliases = {}

for i in aliasFile:
	genes = i.split("\t")
	aliases[genes[2]] = genes[1]
aliasFile.close()




for line in inputFile:
	fields = line.split("\t")
	attrs = getAttributes(fields[8])
	pbid = attrs["gene_name"]	
	if pbid in aliases:
		new_geneName = aliases[pbid]	
		line.replace(pbid, new_geneName) 
		print("found")
		ensembl_ncbi.write(pbid+"\t"+new_geneName+"\n")
	output.write(line)	

inputFile.close()
output.close()
ensembl_ncbi.close()





