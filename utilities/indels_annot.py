#!/usr/bin/env python

import re
import pysam
import os.path
import sys
import itertools
from optparse import OptionParser
from collections import defaultdict


def indels(samFile):

	sam = pysam.AlignmentFile(samFile, "r")
	out_file = samFile.split(".")[0]+"_indels.txt"
	output = open(out_file, "w")
 	output.write("#isoform\tindelStart\tindelEnd\tnt\tnearJunction\tjunctionStart\tjunctionEnd\tindelType\n")

	read_names = []
	indelsJunc = []
	indelsTotal = {}


	for read in sam.fetch():
	 	if read.is_unmapped:
	 		continue
		cigarLine = read.cigar;

		#print read.query_name
		# print(read.query_name)
		# print(cigarLine)
		# print(read.pos)
		
		## reading splice junctions and storing information
		pos = read.pos + 1 
		spliceSites = []

		for (cigarType,cigarLength) in cigarLine:
			if (cigarType not in [1,4,5,8,7]):#["I","S","H","X","="]
				pos_end = pos + cigarLength - 1
				#print(str(pos)+"\t"+str(pos_end)+"\t")
				if (cigarType == 3):  #skip
					spliceSites.append([pos, pos_end])	
				pos = pos_end + 1

		#print(spliceSites)
		## reading indels, comparing with splice junctions and writting information
		pos = read.pos + 1

		for (cigarType,cigarLength) in cigarLine:

			if (cigarType not in [1,4,5,8,7]):#["I","S","H","X","="]
				pos_end = pos + cigarLength -1

			if(cigarType == 1): 
					pos_indel = pos
					length_indel = cigarLength
					pos_end_indel = pos
					log = False
					spliceSitesNearIndel = []
					name = str(read.query_name).split("|")[0]

					
					# indels in the sequence
					if name in indelsTotal:
						indelsTotal[name] = indelsTotal[name] + 1
					else:
						indelsTotal[name] = 1

					# indels near spliceSties
					for i in spliceSites:
						#print i
						if any(abs(pos_indel-e)<10 for e in i) or any(abs(pos_end_indel-e)<10 for e in i):
							spliceSitesNearIndel.append(i)
							log = True
					#print(spliceSitesNearIndel)
					#print(str(log))
					if len(spliceSitesNearIndel)==0:
						output.write(str(name)+"\t"+str(pos_indel)+"\t"+str(pos_end_indel)+"\t"+str(length_indel)+"\t"+str(log)+"\t"+"NA"+"\t"+"NA"+"\tinsertion\n")

					else:	
						for j in spliceSitesNearIndel:
							output.write(str(name)+"\t"+str(pos_indel)+"\t"+str(pos_end_indel)+"\t"+str(length_indel)+"\t"+str(log)+"\t"+str(j[0])+"\t"+str(j[1])+"\tinsertion\n")
							indelsJunc.append(("-").join([name,str(j[0]),str(j[1])]))

			if(cigarType == 2): 
					pos_indel = pos
					length_indel = cigarLength
					pos_end_indel = pos_end
					log = False
					spliceSitesNearIndel = []
					name = str(read.query_name).split("|")[0]

					# indels in the sequence
					if name in indelsTotal:
						indelsTotal[name] = indelsTotal[name] + 1
					else:
						indelsTotal[name] = 1

					# indels near spliceSties
					for i in spliceSites:
						#print i
						if any(abs(pos_indel-e)<10 for e in i) or any(abs(pos_end_indel-e)<10 for e in i):
							spliceSitesNearIndel.append(i)
							log = True
					#print(spliceSitesNearIndel)
					#print(str(log))
					# if 	len(spliceSitesNearIndel)>1:
					# 	print(read.query_name)
					# 	print("wola")
					# 	print(spliceSitesNearIndel)
					if len(spliceSitesNearIndel)==0:
						output.write(str(name)+"\t"+str(pos_indel)+"\t"+str(pos_end_indel)+"\t"+str(length_indel)+"\t"+str(log)+"\t"+"NA"+"\t"+"NA"+"\tdeletion\n")

					else:	
						for j in spliceSitesNearIndel:
							output.write(str(name)+"\t"+str(pos_indel)+"\t"+str(pos_end_indel)+"\t"+str(length_indel)+"\t"+str(log)+"\t"+str(j[0])+"\t"+str(j[1])+"\tdeletion\n")
							indelsJunc.append(("-").join([name,str(j[0]),str(j[1])]))

			if (cigarType not in [1,4,5,8,7]):#["I","S","H","X","="]
				pos = pos_end + 1


	sam.close()
	output.close()
	return(indelsJunc, indelsTotal)


if __name__ == "__main__":
    import sys
    indels(sys.argv[1], sys.argv[2])