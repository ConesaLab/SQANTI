#!/usr/bin/env python
# SQANTI: Structural and Quality Annotation of Novel Transcript Isoforms
# Authors: Lorena de la Fuente, Hector del Risco, Cecile Pereira and Manuel Tardaguila

from collections import defaultdict
import getopt
import sys
import os.path, os
import timeit
import subprocess
import argparse
utilitiesPath =  os.path.dirname(os.path.realpath(__file__))+"/utilities/" 
sys.path.insert(0, utilitiesPath)
from rt_switching import rts
from indels_annot import indels
import psutil
import gc
import distutils.spawn

class myQueryTranscripts:

	def __init__(self, transcript, tss_diff, tts_diff, exon_number, length, str_class, subtype=None, genes=None, transcripts=None, chrom=None, strand=None, bite = "FALSE", RT_switching = "FALSE", canonical="canonical", min_cov = "NA",  min_cov_pos = "NA", min_samp_cov="NA", sd = "NA", FL = "NA", nIndels = "NA", nIndelsJunc = "NA", proteinID=None, ORFlen="NA", CDS_start="NA", CDS_end="NA", isoExp = "NA", geneExp = "NA" , coding = "non_coding", refLen = "NA", refExons = "NA", FSM_class = None, percAdownTTS = None):
		
		self.transcript  = transcript
		self.tss_diff    = tss_diff 	
		self.tts_diff    = tts_diff 	
		self.genes 		 = [] 			
		self.transcripts = [] 			
		self.exon_number = exon_number 
		self.length      = length
		self.str_class   = str_class  	# structural classification of the isoform.
		self.chrom       = chrom
		self.strand 	 = strand
		self.subtype 	 = subtype
		self.RT_switching= RT_switching
		self.canonical   = canonical
		self.min_samp_cov = min_samp_cov
		self.min_cov     = min_cov  
		self.min_cov_pos = min_cov_pos 
		self.sd 	     = sd
		self.proteinID   = proteinID
		self.ORFlen      = ORFlen
		self.CDS_start   = CDS_start
		self.CDS_end     = CDS_end
		self.coding      = coding
		self.FL          = FL
		self.nIndels     = nIndels
		self.nIndelsJunc = nIndelsJunc
		self.isoExp      = isoExp
		self.geneExp     = geneExp
		self.refLen      = refLen
		self.refExons    = refExons
		self.FSM_class   = FSM_class
		self.bite        = bite
		self.percAdownTTS = percAdownTTS

	def get_total_diff(self):
		total_diff = abs(int(self.tss_diff))+abs(int(self.tts_diff))
		return total_diff
	
	def modify(self, transcript, gene, tss_diff, tts_diff, refLen, refExons):
		self.transcript = transcript
		self.genes.append(gene)
		self.tss_diff = tss_diff
		self.tts_diff = tts_diff	
		self.transcripts.append(self.transcript)
		self.refLen = refLen
		self.refExons = refExons

	def geneName(self):
		geneName = "_".join(set(self.genes))
		return geneName

	def ratioExp(self):
		if self.geneExp == 0 or self.geneExp == "NA":
			return "NA"
		else:
			ratio = float(self.isoExp)/float(self.geneExp)	
		return(ratio)

	def CDSlen(self):
		if self.coding == "coding":
			return(str(int(self.CDS_end) - int(self.CDS_start) + 1))
		else:
			return("NA")

	def __str__(self):
		return "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (self.chrom, self.strand, str(self.length), str(self.exon_number), str(self.str_class), "_".join(set(self.genes)), self.transcript, str(self.refLen), str(self.refExons), str(self.tss_diff), str(self.tts_diff), self.subtype, self.RT_switching, self.canonical, str(self.min_samp_cov),str(self.min_cov), str(self.min_cov_pos), str(self.sd), str(self.FL), str(self.nIndels), str(self.nIndelsJunc), self.bite, str(self.isoExp), str(self.geneExp), str(self.ratioExp()), self.FSM_class, self.coding, str(self.ORFlen), str(self.CDSlen()), str(self.CDS_start), str(self.CDS_end), str(self.percAdownTTS))

class myQueryProteins:

	def __init__(self, cds_start, cds_end, orf_length, proteinID="NA"):
		
		self.orf_length  = orf_length 	
		self.cds_start   = cds_start 	
		self.cds_end     = cds_end 
		self.proteinID   = proteinID	

class myQueryJunctions:

	def __init__(self, transcript, junctionNumber, chrom, strand, junc_class, donor, acceptor, diff_start, diff_end, genomCoordStart, genomCoordEnd, transcriptCoord, spliceSite, canonical, indel="NA", coverage=None):

		self.transcript      = transcript
		self.junctionNumber  = junctionNumber
		self.chrom           = chrom 
		self.strand          = strand
		self.junc_class      = junc_class   
		self.donor           = donor
		self.acceptor        = acceptor
		self.diff_start      = diff_start 
		self.diff_end        = diff_end
		self.genomCoordStart = genomCoordStart
		self.genomCoordEnd   = genomCoordEnd
		self.transcriptCoord = transcriptCoord
		self.spliceSite      = spliceSite 
		self.canonical       = canonical
		self.coverage        = []
		self.indel           = indel

	def totalCov(self):
		if "NA" in self.coverage or len(self.coverage)==0:
			totalCov = "NA"
			repCov = "NA"
		else:
			totalCov = sum([int(i) for i in self.coverage])
			repCov = len([i for i in self.coverage if int(i) > 0])
		return(totalCov, repCov)		

	def coverage_files(self):
		if len(self.coverage)==0:
			return("")
		else:
			return("\t"+("\t").join(self.coverage))

	def biteDef(self):
		#if self.diff_donor != 0 and self.diff_acceptor != 0 and self.diff_donor%3 != 0 and self.diff_donor%3 != 0:
		if self.diff_start not in  [0, "NA"] and self.diff_end  not in  [0, "NA"]:
			return ("TRUE")
		else:
			return ("FALSE")

	def __str__(self):

		return "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s%s" % (self.transcript, self.chrom, self.strand, self.junctionNumber, self.genomCoordStart, self.genomCoordEnd, str(self.transcriptCoord), self.junc_class, self.donor, self.acceptor, self.diff_start, self.diff_end, self.biteDef(), self.spliceSite, self.canonical, self.indel, self.totalCov()[1], self.totalCov()[0], self.coverage_files())

class myQueries:

	def __init__(self, line, coordToSort, assoc=None): #CoordToSort can be the start of the first exon or the start of the first junction
		self.line = line
		self.coordToSort = coordToSort
		self.assoc = assoc

class refTranscripts:

	def __init__(self, transcript, gene, chrom, strand, tss=None, tts=None): 
		self.transcript        = transcript
		self.chrom             = chrom
		self.strand            = strand
		self.junctions         = []           #inside two list: starts of junctions and ends of junctions (sorted by coordinates)
		self.tss               = tss
		self.tts               = tts
		self.gene              = gene
		self.exons             = []           #inside two list: starts of exons and ends of exons (sorted by coordinates)
		self.exons_coordinates = []
	
	def add_information(self, line, strand):

		exon_s = line[8].split(",")[0:-1]
		exon_e = line[9].split(",")[0:-1]
		exon_s = [int(i)+1 for i in exon_s]

		#exons 

		self.exons.append(exon_s)
		self.exons.append(exon_e)

		# tss and tts

		if strand == "+":
			tss = exon_s[0]  
			tts = exon_e[-1]	
		else:
			tts = exon_s[0]
			tss = exon_e[-1]

		self.tss = int(tss)
		self.tts = int(tts)

		#junctions

		if len(exon_e) >= 2: #junctions only if multiple-exon isoforms
			junctions_e = exon_s[1:] 
			junctions_e = [int(i)-1 for i in junctions_e]
			junctions_s = exon_e[:-1]	
			junctions_s = [int(i)+1 for i in junctions_s]
			self.junctions.append(junctions_s)
			self.junctions.append(junctions_e)


	def difference_fragment_transcript(self, coord_start, coord_end, ref_coord): 

		if int(coord_start) >= int(ref_coord[0]):
			lack_start = ref_coord.index(int(coord_start))
		else:
			lack_start = int(coord_start) - int(ref_coord[0]) 

		if int(coord_end)  <= int(ref_coord[-1]):
			lack_end  = len(ref_coord) - ref_coord.index(int(coord_end)) -1
		
		else:
			lack_end = int(ref_coord[-1]) - int(coord_end) 

		return (lack_start, lack_end)	


	def transcriptLength(self):
		length=0
		for i in range(len(self.exons[0])):
			length = length + (int(self.exons[1][i])-int(self.exons[0][i])+1)
		return length

class refGenes:

	def __init__(self, gene, chrom, strand):
		self.chrom              = chrom
		self.junctions          = []   #Here we'll save vectors of two elements (start and end of the junction). #The final junctions vector will have as many elements as different juctions exist in that gene.
		self.strand             = strand			
		self.tss                = []
		self.tts                = []
		self.exons              = []   # inside each exon of the gene is in a list (start-end).  #The final junctions vector will have as many elements as different exons exist in that gene.
		self.last3junctions     = []   #all last3junctions for that gene
		self.last5junctions     = []   #all last5junctions for that gene
		self.last3exons         = []   #all last3exons for that gene
		self.last5exons         = []   #all last3exons for that gene
		self.gene               = gene
		self.transcripts        = []
		self.exonsTranscript    = []   #[[[exon1_t1][exon2_t2]], [[exon1_t2][exon2_t2]],....) # Now you can see all the exons for each transcript of the gene
		self.exonsTranscriptTog = []

	def add_transcript(self, line, strand):

		# exons in genomic coordinates

		exon_s = line[8].split(",")[0:-1]
		exon_e = line[9].split(",")[0:-1]
		exon_s = [int(i)+1 for i in exon_s]
		exon_e = [int(i) for i in exon_e]

		subject = line[0]

		exons_per_transcript = []
		exons_per_transcript_total = []

		for i in range(len(exon_s)):
			self.exons.append([exon_s[i],exon_e[i]])
			exons_per_transcript_total.append([exon_s[i],exon_e[i]]) #el primer y ultimo exon lo tenemos en cuenta. ##### Mirar para que lo necesitamos realmente
			if i != 1 and i != len(exon_s):
				exons_per_transcript.append([exon_s[i],exon_e[i]]) #el primer y ultimo exon no lo tenemos en cuenta. ##### Mirar para que lo necesitamos realmente

		self.exonsTranscriptTog.append(exons_per_transcript_total)		#el primer y ultimo exon se tienen en cuenta.
		self.exonsTranscript.append(exons_per_transcript)	#el primer y ultimo exon no se tienen en cuenta.
		self.transcripts.append(subject) #adding transcripts names	
		
		if len(exon_s) >= 2: #last exons only if multiple-exon isoforms
			if strand == "+":
				self.last5exons.append([exon_s[0], exon_e[0]])
				self.last3exons.append([exon_s[-1], exon_e[-1]])
			if strand == "-":
				self.last5exons.append([exon_s[-1], exon_e[-1]])
				self.last3exons.append([exon_s[0], exon_e[0]])

		# tss and tts

		if strand == "+":
			tss = exon_s[0]  
			tts = exon_e[-1]	
		else:
			tts = exon_s[0]
			tss = exon_e[-1]

		self.tss.append(int(tss))
		self.tts.append(int(tts))


		#junctions

		if len(exon_s) >= 2: #junctions only if multiple-exon isoforms
			junctions_e = exon_s[1:] 
			junctions_e = [int(i)-1 for i in junctions_e]
			junctions_s = exon_e[:-1]	
			junctions_s = [int(i)+1 for i in junctions_s]

			for i in range(len(junctions_s)):
				self.junctions.append([junctions_s[i],junctions_e[i]])

			if strand == "+":
				self.last5junctions.append([junctions_s[0], junctions_e[0]])
				self.last3junctions.append([junctions_s[-1], junctions_e[-1]])

			if strand == "-":
				self.last5junctions.append([junctions_s[-1], junctions_e[-1]])
				self.last3junctions.append([junctions_s[0], junctions_e[0]])


	def __str__(self):
		return "%s\t%s\t%s\t%s" % (self.gene, self.chrom, self.strand, self.exons)

def ReverseComplement(seq):
	basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'} 
	seq = seq.upper()
	RCseq = ""
	for base in seq:
		if base not in basecomplement:
			sys.stderr.write("Error: %s NOT a DNA sequence" %(seq))
			return None
			break
		RCseq = basecomplement[base]+RCseq
	return (RCseq)
def transcriptLength(startExons, endExons):
	length=0
	for i in range(len(startExons)):
		length = length + (int(endExons[i])-int(startExons[i])+1)
	return length
def exons_together(starts, ends): #arguments: list of starts of exons and list of ends of exons
	exons_together = []
	for i in range(len(starts)):
		exons_together = exons_together + range(int(starts[i]), int(ends[i])+1)
	return exons_together
def getTranscriptJunctCoordinates(exons_starts, exons_ends, strand, junctions=None):
	exons_unknown_together = []
	starts_transcript = []
	ends_transcript = []
	for i in range(len(exons_starts)):
		exons_unknown_together = exons_unknown_together + range(int(exons_starts[i]), int(exons_ends[i])+1)
	if junctions==None:
		junctions = exons_ends[:-1]
	if strand=="+":
		exons_unknown_together_sort = sorted(exons_unknown_together)
		for l in junctions:	
			for k in range(len(exons_unknown_together_sort)):
				if exons_unknown_together_sort[k] == int(l):
					starts_transcript.append(k+1) 
					ends_transcript.append(k+1+1) #k+1+1
	if strand=="-":
		exons_unknown_together_sort = sorted(exons_unknown_together, reverse=True)
		for l in junctions:	
			for k in range(len(exons_unknown_together_sort)):
				if exons_unknown_together_sort[k] == int(l):
					starts_transcript.append(str(k+1)) 
					ends_transcript.append(str(k)) #k+1-1

	start = min(starts_transcript, ends_transcript)				
	end = max(starts_transcript, ends_transcript)				

	return [start, end]

def fasta_parser(fastaFile):

	try:
		fasta = open(fastaFile, "r")
	except IOError:
		sys.stderr.write('ERROR: Unable to read %s file\n' % fastaFile)
		raise SystemExit(1)
	try:
		seqDicc = {}
		index = 0
		for line in fasta:
			if line.startswith(">"):
				if index > 0:
					seqDicc[name] = seq
				index+=1
				name = line[1:].rstrip()
				#name = line[1:].split()[0].rstrip()
				seq = ''
			elif line.strip()!="":
				seq += line.rstrip()
		if index>0:
			seqDicc[name] = seq
		fasta.close()
	except (RuntimeError, TypeError, NameError):
		sys.stderr.write('ERROR: File %s without fasta format.\n' % fastaFile)
		raise SystemExit(1)

	return(seqDicc)

def genome_parser(fastaFile):

	try:
		fasta = open(fastaFile, "r")
	except IOError:
		sys.stderr.write('ERROR: Unable to read %s file.\n' % fastaFile)
		raise SystemExit(1)
	try:
		seqDicc = {}
		index = 0
		for line in fasta:
			if line.startswith(">"):
				if index > 0:
					seqDicc[name] = seq
				index+=1
				name = line[1:].rstrip().split()[0]
				#name = line[1:].split()[0].rstrip()
				seq = ''
			elif line.strip()!="":
				seq += line.rstrip()
		seqDicc[name] = seq
		fasta.close()
	except (RuntimeError, TypeError, NameError):
		sys.stderr.write('ERROR: File %s without fasta format.\n' % fastaFile)
		raise SystemExit(1)

	return(seqDicc)

def correctionPlusORFpred(args, chrom_names):

	global corrORF
	global corrGTF
	global corrSAM
	global corrFASTA

	corrPathPrefix = args.dir+"/"+os.path.splitext(os.path.basename(args.isoforms))[0]
	corrGFF_gmap = corrPathPrefix +"_corrected_gmap.gff"
	corrGFF = corrPathPrefix +"_corrected.gff"
	corrGTF = corrPathPrefix +"_corrected.gtf"
	corrSAM = corrPathPrefix +"_corrected.sam"
	corrFASTA = corrPathPrefix +"_corrected.fasta"
	ORF =  corrPathPrefix +"_corrected.faa"
	corrORF = os.path.splitext(ORF)[0]+'_ATG.faa'


	if not args.gtf:

		if (os.path.isdir(os.path.dirname(args.gmap_index))): # check index folder
			index_dir = os.path.dirname(args.gmap_index)
			prefix =  os.path.basename(args.gmap_index)

			# modifiying fasta header to keep just transcript ID.
			fastaDicc = fasta_parser(args.isoforms)	
			tmpFasta = os.path.splitext(args.isoforms)[0]+'.tmp'
			tmpFasta_file = open(tmpFasta, "w")
			for ID in fastaDicc:
				if len(ID.split("|"))>2 and ID[0:2]=="PB": # PacBio fasta header (including chained format)
					ID_mod = ID.split("|")[0].split(" ")[0]
				elif len(ID.split("|"))>4: # Refseq fasta header
					ID_mod = ID.split("|")[3]
				else:
					ID_mod = ID.split()[0] # Ensembl fasta header

				tmpFasta_file.write(">"+ID_mod+"\n"+fastaDicc[ID]+"\n")
			tmpFasta_file.close()
			
			sys.stdout.write("\n****Aligning reads with gmap...\n")


			# aligning sequences
			try:
				with open(corrGFF_gmap, 'w') as corrGFF_out:
					subprocess.call(['gmap','-n', '1', '-t', '4', '--cross-species','--min-intronlength=4', '--gff3-add-separators=0','-f', '2', '-D', index_dir,'-d', prefix, tmpFasta], stdout=corrGFF_out)
				with open(corrSAM, 'w') as corrSAM_out:
					subprocess.call(['gmap','-n', '1', '-t', '4', '--cross-species', '--min-intronlength=4', '--gff3-add-separators=0','-f', 'samse', '-D', index_dir,'-d', prefix, tmpFasta], stdout=corrSAM_out)
				os.remove(tmpFasta)
			except (RuntimeError, TypeError, NameError):
				sys.stderr.write('ERROR: Problem during alignment of sequences.\n')
				raise SystemExit(1)



			#GFF GMAP to accurate GFF
			with open(corrGFF, 'w') as corrGFF_out:
				with open(corrGFF_gmap, 'r') as corrGFF_in:
					for i in corrGFF_in:
						if i[0] != "#" and i.split("\t")[2] != "CDS":
							if (".mrna1" in i) or (".path1" in i):
								j = i.replace(".mrna1","")
								k = j.replace(".path1", ".gene")
								corrGFF_out.write(k)

			# GFF to GTF
			subprocess.call([utilitiesPath+"gffread", corrGFF , '-T', '-o', corrGTF])

			os.remove(corrGFF_gmap)
			os.remove(corrGFF)
			del fastaDicc


		else:
			sys.stderr.write("ERROR: '%s' directory containing gmap indexes doesn't exist.\n" % (os.path.dirname(args.gmap_index)))
			sys.exit()



	if args.gtf:

		sys.stdout.write("\nSkipping aligning of sequences because gtf file was provided.\n")

		ind = 0
		with open(args.isoforms, 'r') as isoforms_gtf:
			for line in isoforms_gtf:
				if line[0] != "#" and len(line.split("\t"))!=9:
					sys.stderr.write("\nERROR: input isoforms file with not gtf format.\n")
					sys.exit()
				elif len(line.split("\t"))==9:
					ind += 1
			if ind==0:
					sys.stderr.write("\nERROR: gtf has not annotation lines.\n")
					sys.exit()


		# GFF to GTF (in case the user provides gff instead of gtf)
		corrGTF_tpm = corrGTF+".tmp"
		try:
			subprocess.call([utilitiesPath+"gffread", args.isoforms , '-T', '-o', corrGTF_tpm])
		except (RuntimeError, TypeError, NameError):
			sys.stderr.write('ERROR: File %s without gtf/gff format.\n' % args.isoforms)
			raise SystemExit(1)


		# check if gtf chromosomes inside genome file
		with open(corrGTF, 'w') as corrGTF_out:
			with open(corrGTF_tpm, 'r') as isoforms_gtf:
				for line in isoforms_gtf:
					if line[0] != "#":
						chrom = line.split("\t")[0]
						type = line.split("\t")[2]
						if chrom not in chrom_names:
							sys.stderr.write("\nERROR: gtf \"%s\" chromosome not found in genome reference file.\n" % (chrom))
							sys.exit()
						elif type=="exon":
							corrGTF_out.write(line)
		os.remove(corrGTF_tpm)

		if not os.path.exists(corrSAM):
			sys.stdout.write("\nIndels will be not calculated since you ran SQANTI without alignment step (SQANTI with gtf format as transcriptome input).\n")


	# Corrected gtf to fasta
	subprocess.call([utilitiesPath+"gffread", corrGTF, '-g', args.genome, '-w', corrFASTA]) 

	# ORF generation
	sys.stdout.write("\n**** Predicting ORF sequences...\n")

	gmst_dir = args.dir+"/GMST/"
	if not os.path.exists(gmst_dir):
	 	os.makedirs(gmst_dir)

	subprocess.call(["perl", utilitiesPath+"gmst/gmst.pl", corrFASTA , '-faa', '--strand','direct', '--fnn', '--output',  os.path.splitext(os.path.basename(args.isoforms))[0]+"_corrected"], cwd=gmst_dir) 
	os.rename(gmst_dir+os.path.basename(ORF), ORF) 
	gc.collect()

	# Modifying ORF sequences by removing sequence before ATG

	ORFdicc = fasta_parser(ORF)	
	corrORF_file = open(corrORF, "w")

	orfLenDicc = {}
	for ID in ORFdicc:
		seq = ORFdicc[ID]
		l = len(ID.split("\t"))
		pos = seq.find('M')
		if pos != -1:
			corr_seq = seq[pos:]
			cds_start = int(ID.split("\t")[l-1].split("|")[4])+pos*3
			cds_end = int(ID.split("\t")[l-1].split("|")[5])  
			strand = ID.split("\t")[l-1].split("|")[3]
			ID_mod = ID.split()[0].strip()+"\t"+"isoform|GeneMark.hmm|"+ str(len(corr_seq)) + "_aa|" + strand + "|" + str(cds_start) +"|" + str(cds_end)
			corrORF_file.write(">"+ID_mod+"\n"+corr_seq+"\n")
			orfLenDicc[ID.split()[0].strip()] = myQueryProteins(cds_start, (cds_start+len(corr_seq)*3-1), len(corr_seq))
	
	corrORF_file.close()

	if len(orfLenDicc)==0:
 		sys.stderr.write("WARNING: All input isoforms were predicted as non-coding")

	
	os.remove(ORF)
	os.rename(corrORF, ORF) 

	del ORFdicc
	gc.collect()
	return(orfLenDicc)


def reference_parser(args, chrom_names):

	global referenceFiles
	referenceFiles = os.path.splitext(args.annotation)[0] +".genePred"

 	sys.stdout.write("\n**** Parsing Reference Transcriptome...\n")


	ind = 0
	with open(args.annotation, 'r') as reference_gtf:
		for line in reference_gtf:
			if line[0] != "#" and len(line.split("\t"))!=9:
				sys.stderr.write("\nERROR: reference annotation file with not gtf/gff format.\n")
				sys.exit()
			elif len(line.split("\t"))==9:
				ind += 1
		if ind==0:
				sys.stderr.write("\nERROR: reference annotation file has not annotation lines.\n")
				sys.exit()


	## gtf to genePred
 	if args.name:
 		subprocess.call([utilitiesPath+"gtfToGenePred", args.annotation, referenceFiles, '-genePredExt', '-allErrors', '-ignoreGroupsWithoutExons', '-geneNameAsName2'])
 	else:
 		subprocess.call([utilitiesPath+"gtfToGenePred", args.annotation, referenceFiles, '-genePredExt', '-allErrors', '-ignoreGroupsWithoutExons'])

 	if not os.path.exists(referenceFiles):
		sys.exit()


	## parsing annotation while filtering miRNA RNAs  

	# dicc of transcripts. We differenciate between single-exon transcripts and multi-exon transcripts.
	transcripts_chrom_1exon = {}
	transcripts_chrom_exons = {}

	# dicc of genes
	transcripts_gene = {}
	missingChroms = []
	referenceFiles_tpm = referenceFiles+"notmiRNAtmp"
	flag=True
	with open(referenceFiles_tpm, 'w') as inputAnnotationNotmiRNAs:
		with open(referenceFiles, 'r') as inputAnnotation:
			for i in inputAnnotation:
				line_split = i.split("\t")
				if len(line_split[8].split(",")) < 3:
					if (int(line_split[9].split(",")[0]) - int(line_split[8].split(",")[0])) > 200:
						inputAnnotationNotmiRNAs.write(i)
						flag=True	
					else:
						flag=False	
				else:
						inputAnnotationNotmiRNAs.write(i)
						flag=True	

				if flag==True:
					transcript_id = line_split[0]
					chrom = line_split[1]
					if (chrom not in chrom_names) and (chrom not in missingChroms):
						missingChroms.append(chrom)
					strand = line_split[2]
					gene = line_split[11]

					if gene == "":
						continue

					#add gene info
					if gene not in transcripts_gene:
							transcripts_gene[gene] = refGenes(gene, chrom, strand)
					transcripts_gene[gene].add_transcript(line_split, strand)
					
					#add transcript info
					if len(line_split[8].split(",")) > 2: # multiple-exon isoform
						if chrom not in transcripts_chrom_exons:
							transcripts_chrom_exons[chrom]=[]

						a = refTranscripts(transcript_id, gene, chrom, strand)
						a.add_information(line_split, strand)
						transcripts_chrom_exons[chrom].append(a)

					else: #single-exon isoform
						if chrom not in transcripts_chrom_1exon:
							transcripts_chrom_1exon[chrom]=[]

						a = refTranscripts(transcript_id, gene, chrom, strand)
						a.add_information(line_split, strand)
						transcripts_chrom_1exon[chrom].append(a)

	os.rename(referenceFiles_tpm, referenceFiles) 

	if len(missingChroms)>0:
		sys.stderr.write("\nWARNING: %s sequences not found in genome but containing annotated genes\n" %(",".join(missingChroms)))

	 
	## sorting transcript information by first junction (if mono-exon transcript we sort by first exon)

	sum=0
	for chrom in transcripts_chrom_exons:
		transcripts_chrom_exons[chrom] = sorted(transcripts_chrom_exons[chrom], key = lambda tup:int(tup.junctions[0][0]))
		sum = len(transcripts_chrom_exons[chrom]) + sum
	
	for chrom in transcripts_chrom_1exon:
		transcripts_chrom_1exon[chrom] = sorted(transcripts_chrom_1exon[chrom], key = lambda tup:int(tup.exons[0][0]))
		sum = len(transcripts_chrom_1exon[chrom]) + sum

	## dicc for chomosome and inside gene information. Got from transcripts_gene dicc.
	transcripts_gene_chrom = {}

	for gene in transcripts_gene.keys():
		chrom = transcripts_gene[gene].chrom
		if chrom not in transcripts_gene_chrom:
			transcripts_gene_chrom[chrom]=[]
		transcripts_gene_chrom[chrom].append(transcripts_gene[gene])

	for chrom in transcripts_gene_chrom:
		transcripts_gene_chrom[chrom] = sorted(transcripts_gene_chrom[chrom], key = lambda tup:min(tup.tss + tup.tts))   ##### COMPROBAR QUE ESTA BIEN!!!!

	return (transcripts_chrom_1exon, transcripts_chrom_exons, transcripts_gene, transcripts_gene_chrom)

def isoforms_parser(args):

	global queryFile
 	queryFile = os.path.splitext(corrGTF)[0] +".genePred"

 	sys.stdout.write("\n**** Parsing Isoforms...\n")

 	# gtf to genePred
 	try:
 		subprocess.call([utilitiesPath+"gtfToGenePred", corrGTF, queryFile , '-genePredExt', '-allErrors', '-ignoreGroupsWithoutExons'])
	except (RuntimeError, TypeError, NameError):
		sys.stderr.write('ERROR: Check isoforms file.\n' )
		raise SystemExit(1)

 	## analysis of query transcripts	
 	queryAnnotation = open(queryFile, 'r')

 	myQueryTranscripts_list = []
	
 	for line in queryAnnotation:

 		exon_s = line.split("\t")[8].split(",")[0:-1]
 		exon_e = line.split("\t")[9].split(",")[0:-1]

 		if len(exon_s) >= 2: #junctions only if multiple-exon isoforms
 			coordToSort = int(exon_e[0]) #first genomic splice site

 		else:
 			coordToSort = int(exon_s[0]) #first exon position

 		a = myQueries(line, coordToSort)
 		myQueryTranscripts_list.append(a)

 	queryAnnotation.close()

 	# sorting the subject transcript by start junction or start exon (depending if one or more exons)
	
 	myQueryTranscripts_list = sorted(myQueryTranscripts_list, key = lambda tup:int(tup.coordToSort))

 	return(myQueryTranscripts_list)

def STARcov_parser(coverageFiles): # just valid with unstrand-specific RNA-seq protocols. 

	#may be one file, files separated by comma or a directory where file are.
	if os.path.isdir(coverageFiles)==True:
		cov_paths = [os.path.join(coverageFiles,fn) for fn in next(os.walk(coverageFiles))[2]]
	else: 
		cov_paths = coverageFiles.split(",")

	cov_list_dicc = {}
	for path in cov_paths:
		dicc = {}
		p = open(path.strip(), "r")
		for line in p:
			j = "_".join(line.split()[0:3])
			dicc[j] = str(int(line.split()[6])+int(line.split()[7]))
		cov_list_dicc[os.path.basename(path)] = dicc
		p.close()

	return(cov_list_dicc)

def expression_parser(expressionFile):

	try:
		p = open(expressionFile, "r")
	except IOError:
		sys.stderr.write('ERROR: Unable to read %s expression file\n' % expressionFile)
		raise SystemExit(1)
	try:
		header = p.readline()
		exp_dicc = {}

		for line in p:
			pbid = line.split()[0]
			mean = sum([float(i) for i in line.rstrip().split()[1:]])/len(line.rstrip().split()[1:])
			exp_dicc[pbid] = mean
		p.close()
	
	except IOError:
		sys.stderr.write('File %s without expression matrix format' % expressionFile)
		raise SystemExit(1)

	return(exp_dicc)	
def transcriptsKnownSpliceSites(transcripts_chrom_1exon, transcripts_chrom_exons, line_split, index_chrom_exons, index_chrom_1exon, transcripts_gene, seq, nPolyA):

	# Transcript information for a single query transcript and comparison with reference.

	query_transcript = line_split[0]
	CHROM = line_split[1]
	STRAND = line_split[2]

	EXON_s = [int(i)+1 for i in line_split[8].split(",")[0:-1]]
	EXON_e = [int(i) for i in line_split[9].split(",")[0:-1]]

	isoLength = transcriptLength(EXON_s, EXON_e)

	JUNCTIONS_e = [int(i)-1 for i in EXON_s[1:] ]
	JUNCTIONS_s = [int(i)+1 for i in EXON_e[:-1]]

	ends_positions = []
	myTranscript_Assoc="" 

	JUNCTIONS = []
	[JUNCTIONS.append([JUNCTIONS_s[m],JUNCTIONS_e[m]]) for m in range(len(JUNCTIONS_s))]			
	EXONS_COORD = exons_together(EXON_s, EXON_e)


	# # Intra-priming

	# if STRAND == "+":
	# 	TTS = int(EXON_e[-1])
	# 	downTTS = str(seq[CHROM].seq[TTS:(int(EXON_e[-1])+nPolyA)])
	# 	downTTS = downTTS.upper()
	# 	if any([base not in ["A","T","C","G","N"] for base in downTTS]):
	# 		sys.stderr.write("Error: %s in 3' TTS genomic downstream region is NOT a DNA sequence in genome file" %(downTTS))
	# 		percA = 0

	# elif STRAND == "-":
	# 	TTS = int(EXON_s[0])
	# 	downTTS = ReverseComplement(str(seq[CHROM].seq[int(EXON_s[0]-nPolyA):TTS]))
	# percA = downTTS.count('A')/nPolyA*100
	
	if STRAND == "+":
		TTS = int(EXON_e[-1])
		downTTS = str(seq[CHROM][TTS:(int(EXON_e[-1])+nPolyA)])
		downTTS = downTTS.upper()
		if any([base not in ["A","T","C","G","N"] for base in downTTS]):
			sys.stderr.write("Error: %s in 3' TTS genomic downstream region is NOT a DNA sequence in genome file" %(downTTS))
			percA = 0

	elif STRAND == "-":
		TTS = int(EXON_s[0])
		downTTS = ReverseComplement(str(seq[CHROM][int(EXON_s[0]-(nPolyA+1)):(TTS-1)]))
	percA = float(downTTS.count('A'))/nPolyA*100


	myTranscript_Assoc = myQueryTranscripts("", "NA","NA",len(EXON_s),isoLength, "", chrom=CHROM, strand=STRAND, subtype="no_subcategory",  percAdownTTS=str(percA))



	##***************************************##
	########### SPLICED TRANSCRIPTS ###########
	##***************************************##

	if len(EXON_s) >= 2: #junctions only if multiple-exon isoforms
		if CHROM in transcripts_chrom_exons:
			if CHROM not in index_chrom_exons:
				index_chrom_exons[CHROM] = 0
			
			for k in range(index_chrom_exons[CHROM], len(transcripts_chrom_exons[CHROM])):
				subtype = "no_category"
				ref_transcript = transcripts_chrom_exons[CHROM][k]
				junctions_s = ref_transcript.junctions[0]
				junctions_e = ref_transcript.junctions[1]
				junctions = []
				[junctions.append([junctions_s[n],junctions_e[n]]) for n in range(len(junctions_s))]
				ends_positions.append(int(junctions_e[-1]))


				########### SAME STRAND ###########

				if ref_transcript.strand == STRAND:

					if all(item  < int(JUNCTIONS_s[0]) for item in ends_positions) :
						index_chrom_exons[CHROM] = index_chrom_exons[CHROM] + 1

					elif int(junctions_s[0]) > int(JUNCTIONS_e[-1]):
						break

					else:

						#######################################
						##### 1. Check full-splice match ######
						#######################################
						
						if (JUNCTIONS_e==junctions_e and JUNCTIONS_s==junctions_s): #splice junctions are identical
							tss = int(ref_transcript.tss)
							tts = int(ref_transcript.tts)

							if STRAND == "+":
								TSS = EXON_s[0]  
								TTS = EXON_e[-1]	
								tss_diff = int(TSS) - tss
								tts_diff = tts - int(TTS)  
							else:
								TTS = EXON_s[0]
								TSS = EXON_e[-1]
								tss_diff = tss - int(TSS)
								tts_diff = int(TTS) - tts

							if myTranscript_Assoc.str_class!="full-splice_match":
								myTranscript_Assoc = myQueryTranscripts(ref_transcript.transcript, tss_diff, tts_diff, len(EXON_s), isoLength, "full-splice_match", subtype="multi-exon", chrom=CHROM, strand=STRAND, refLen = ref_transcript.transcriptLength(), refExons= len(ref_transcript.exons[0]), percAdownTTS=str(percA))		
								myTranscript_Assoc.transcripts = []   
								myTranscript_Assoc.transcripts.append(ref_transcript.transcript)
								myTranscript_Assoc.genes = []   
								myTranscript_Assoc.genes.append(ref_transcript.gene)

							elif  (abs(tss_diff) + abs(tts_diff)) < myTranscript_Assoc.get_total_diff():
								myTranscript_Assoc.genes = []  
								myTranscript_Assoc.modify(ref_transcript.transcript, ref_transcript.gene, tss_diff, tts_diff, ref_transcript.transcriptLength(), len(ref_transcript.exons[0]))
							else:
								myTranscript_Assoc.transcripts.append(ref_transcript.transcript)
	

						####################################################################
						##### 2. Check fragment-splice match if none full-splice match ######
						####################################################################

						elif myTranscript_Assoc.str_class!="full-splice_match" and all(s in junctions for s in JUNCTIONS): #check if it's a fragment of a reference transcript if not perfect match found. Just in case you have not found before another perfect match (because we don't break the loop if found a perfect match since we want the similar one (utr differences)) 

							exons_coord = exons_together(ref_transcript.exons[0], ref_transcript.exons[1])
							differential_nt = list(set(EXONS_COORD) - set(exons_coord)) 
							longerUTRs = all((int(d) < int(ref_transcript.exons[0][0]) or int(d) > int(ref_transcript.exons[1][-1]) ) for d in differential_nt)  # extra-nt in query are for differences in utrs
							
							if len(differential_nt)==0 or longerUTRs:
								if STRAND == "+":
									(tss_diff, tts_diff) = ref_transcript.difference_fragment_transcript(EXON_s[0], EXON_e[-1], exons_coord)
									LAST5JUNCTION = [JUNCTIONS_s[0], JUNCTIONS_e[0]]
									LAST3JUNCTION = [JUNCTIONS_e[-1], JUNCTIONS_e[-1]]
									last5junction = [junctions_s[0], junctions_e[0]]
									last3junction = [junctions_s[-1], junctions_e[-1]]
								else:
									(tts_diff, tss_diff) = ref_transcript.difference_fragment_transcript(EXON_s[0], EXON_e[-1], exons_coord)
									LAST3JUNCTION = [JUNCTIONS_s[0], JUNCTIONS_e[0]]
									LAST5JUNCTION = [JUNCTIONS_e[-1], JUNCTIONS_e[-1]]
									last3junction = [junctions_s[0], junctions_e[0]]
									last5junction = [junctions_s[-1], junctions_e[-1]]

								if (LAST5JUNCTION!=last5junction) and (LAST3JUNCTION != last3junction):
									subtype = "internal_fragment"
			
								elif LAST5JUNCTION != last5junction:
									subtype = "5prime_fragment"

								elif LAST3JUNCTION != last3junction:
									subtype = "3prime_fragment"

								else:
									subtype = "complete"

								if myTranscript_Assoc.str_class!="incomplete-splice_match":
									myTranscript_Assoc = myQueryTranscripts(ref_transcript.transcript, tss_diff, tts_diff, len(EXON_s), isoLength,"incomplete-splice_match", subtype=subtype, chrom=CHROM, strand=STRAND, refLen = ref_transcript.transcriptLength(), refExons= len(ref_transcript.exons[0]), percAdownTTS=str(percA))
									myTranscript_Assoc.transcripts.append(ref_transcript.transcript)
									myTranscript_Assoc.genes = []  
									myTranscript_Assoc.genes.append(ref_transcript.gene)
								elif (abs(tss_diff) + abs(tts_diff)) < myTranscript_Assoc.get_total_diff():
									myTranscript_Assoc.genes = []   
									myTranscript_Assoc.modify(ref_transcript.transcript, ref_transcript.gene, tss_diff, tts_diff, ref_transcript.transcriptLength(), len(ref_transcript.exons[0]))	
								else:
									myTranscript_Assoc.transcripts.append(ref_transcript.transcript)

						########################################################################
						##### 3. Association of transcripts to genes by sharing junctions ######
						########################################################################

						if myTranscript_Assoc.str_class!="full-splice_match" and myTranscript_Assoc.str_class!="incomplete-splice_match" and (len(set(junctions_s).intersection(set(JUNCTIONS_s)))>0 or len(set(junctions_e).intersection(set(JUNCTIONS_e)))>0):  # Si no es ni fragment ni perfect match miramos si comparten junctions (asi definiremos el gen al que pertence o si es fusion)
							gene = ref_transcript.gene
							type = "anyKnownSpliceSite"
							if myTranscript_Assoc.str_class=="":
									myTranscript_Assoc = myQueryTranscripts("novel", "NA", "NA", len(EXON_s), isoLength, type, chrom=CHROM, strand=STRAND, subtype="no_subcategory", percAdownTTS=str(percA)) 
							myTranscript_Assoc.genes.append(ref_transcript.gene)


	##***************************************####
	########### UNSPLICED TRANSCRIPTS ###########
	##***************************************####

	else: #If just one exon
		if CHROM in transcripts_chrom_1exon:
			if CHROM not in index_chrom_1exon:
				index_chrom_1exon[CHROM] = 0
			
			for k in range(index_chrom_1exon[CHROM], len(transcripts_chrom_1exon[CHROM])):
					ref_transcript = transcripts_chrom_1exon[CHROM][k]
					
					if ref_transcript.strand == STRAND:
						exon_s = ref_transcript.exons[0]
						exon_e = ref_transcript.exons[1]
						if int(exon_e[-1]) < int(EXON_s[0]):
							index_chrom_1exon[CHROM] = index_chrom_1exon[CHROM] + 1
						elif int(exon_s[0]) > int(EXON_e[-1]):
							break
					# Checking if junctions are identical.
						else:
							tss = int(ref_transcript.tss)
							tts = int(ref_transcript.tts)
							if STRAND == "+":
								TSS = EXON_s[0]  
								TTS = EXON_e[-1]	
								tss_diff = int(TSS) - tss
								tts_diff = tts - int(TTS)  
								rangeRef = set(range(int(tss), int(tts)))
								rangeSub = set(range(int(TSS), int(TTS)))
							else:
								TTS = EXON_s[0]
								TSS = EXON_e[-1]
								tss_diff = tss - int(TSS)
								tts_diff = int(TTS) - tts
								rangeRef = set(range(int(tts), int(tss)))
								rangeSub = set(range(int(TTS), int(TSS)))
							intersect = rangeRef.intersection(rangeSub)
							
							if len(intersect)>0:
								if myTranscript_Assoc.str_class=="":
									myTranscript_Assoc = myQueryTranscripts(ref_transcript.transcript, tss_diff, tts_diff, "1", isoLength, "full-splice_match", subtype="mono-exon", chrom=CHROM, strand=STRAND, refLen=ref_transcript.transcriptLength(), refExons = len(ref_transcript.exons[0]), percAdownTTS=str(percA))
									myTranscript_Assoc.transcripts.append(ref_transcript.transcript)
									myTranscript_Assoc.genes.append(ref_transcript.gene)
								elif (abs(tss_diff) + abs(tts_diff)) < myTranscript_Assoc.get_total_diff():
									myTranscript_Assoc.genes = []   # change!
									myTranscript_Assoc.modify(ref_transcript.transcript, ref_transcript.gene, tss_diff, tts_diff, ref_transcript.transcriptLength(), len(ref_transcript.exons[0]))
								else:
									myTranscript_Assoc.transcripts.append(ref_transcript.transcript)

					
	myJunctions = "" 
	
	return(myTranscript_Assoc, myJunctions, index_chrom_exons, index_chrom_1exon)
def novelIsoformsKnownGenes(assoc, line_split, transcripts_gene):

	query_transcript = line_split[0]
	CHROM = line_split[1]
	STRAND = line_split[2]

	EXON_s = line_split[8].split(",")[0:-1]
	EXON_e = line_split[9].split(",")[0:-1]
	EXON_s = [int(i)+1 for i in EXON_s]
	EXON_e = [int(i) for i in EXON_e]

	JUNCTIONS_e = [int(i)-1 for i in EXON_s[1:] ]
	JUNCTIONS_s = [int(i)+1 for i in EXON_e[:-1]]

	# Si mas de dos genes puede ser fusion gene o genes overlapping

	if len(set(assoc.genes))>1:
		genes = list(set(assoc.genes))
		#Si algun donor/acceptor esta en varios genes no puede ser fusion:
		sites = sum([list(set(sum(transcripts_gene[gene].junctions, []))) for gene in genes], [])
		if any(t > 1 for t in [sites.count(x) for x in (JUNCTIONS_s+JUNCTIONS_e)]): # at least one SITE in more than one gene, so its not a fusion gene
				number_shared_ac_don={}
				for gene in genes:
					number_shared_ac_don[gene] = len(set(JUNCTIONS_s+JUNCTIONS_e).intersection(set(sum(transcripts_gene[gene].junctions, []))))
				maxN = max(number_shared_ac_don.values())
				g = [name for name,j in number_shared_ac_don.iteritems() if j == maxN]
				if len(g)==1: # El que tenga mas junctions del trascrito
					assoc.genes = g
					assoc.str_class = "moreJunctions"
				else: #El que se salga menos de los limites definidos por ambos genes.
					diff = {}
					for gene in genes:
						diff[gene] = abs(EXON_s[0]- min(min([transcripts_gene[gene].tss, transcripts_gene[gene].tts])))+ abs(EXON_e[-1]- max(max([transcripts_gene[gene].tss, transcripts_gene[gene].tts])))
					g = min(diff, key=diff.get)
					assoc.genes = [g]
					assoc.str_class = "maxdiff"
		else:
			assoc.str_class = "fusion"

	if assoc.str_class!="fusion" and len(set(assoc.genes))==1:
		junctions_g = transcripts_gene[assoc.genes[0]].junctions
		junctions_g_s = []
		junctions_g_e = []
		for i in range(len(junctions_g)):
			junctions_g_s.append(junctions_g[i][0])
			junctions_g_e.append(junctions_g[i][1])

		s = all(x in junctions_g_s for x in JUNCTIONS_s)
		e = all(x in junctions_g_e for x in JUNCTIONS_e)
		if s==True and e==True:
				assoc.str_class="novel_in_catalog" # known donors and acceptors
				if all([JUNCTIONS_s[k],JUNCTIONS_e[k]] in junctions_g for k in range(len(JUNCTIONS_s))) == True:
					assoc.subtype = "combination_of_known_junctions"  	
				else:
					assoc.subtype = "no_combination_of_known_junctions"
		else:
			assoc.str_class="novel_not_in_catalog"
			assoc.subtype = "any annotated donor/acceptor"

	return(assoc)
def associationOverlapping(assoc, line_split, transcripts_gene_chrom, index_genes, over_nt):

	assoc.subtype = "no_subcategory"
	assoc.str_class = "intergenic"
	assoc.transcript = "novel"

	query_transcript = line_split[0]
	CHROM = line_split[1]
	STRAND = line_split[2]

	EXON_s = [int(i)+1 for i in line_split[8].split(",")[0:-1]]
	EXON_e = [int(i) for i in line_split[9].split(",")[0:-1]]

	ends_positions = []

	if CHROM in transcripts_gene_chrom:

		if CHROM not in index_genes:
			index_genes[CHROM] = 0
			
		for k in range(index_genes[CHROM], len(transcripts_gene_chrom[CHROM])):
			
			ref_gene = transcripts_gene_chrom[CHROM][k]
			ends_positions.append(max(ref_gene.tss + ref_gene.tts))

			if all(item  < int(EXON_s[0]) for item in ends_positions) :
				index_genes[CHROM] = index_genes[CHROM] + 1
		
			elif min(ref_gene.tts + ref_gene.tss) > int(EXON_e[-1]):
					break

			#checking if overlapping of any trancript
			else:
				exons = ref_gene.exons
				exon_s = [s[0] for s in exons]
				exon_e = [e[0] for e in exons]

				for l in range(len(EXON_s)):
					E_s = int(EXON_s[l])
					E_e = int(EXON_e[l])

					for p in range(len(exons)):
						e_s = int(exons[p][0])
						e_e = int(exons[p][1])
						
						if (E_s >= e_s and E_e <= e_e) or (E_s <= e_s and E_e >= e_e) or (E_s<=e_s<=E_e) or (E_s<=e_e<=E_e):
							
							int_nt = len(set(exons_together(exon_s, exon_e)).intersection(set(exons_together(EXON_s,EXON_e))))

							if ref_gene.strand == STRAND:

								if assoc.str_class != "genic":
									assoc.str_class = "genic"
									assoc.genes = [ref_gene.gene]
									over_nt[query_transcript] = int_nt

								else: #We compare first genic with new one
									if int_nt > over_nt[query_transcript]:
										assoc.genes = [ref_gene.gene]
										over_nt[query_transcript] = int_nt

								if assoc.exon_number > 1:
									min_coordinate = min(ref_gene.tss + ref_gene.tts)
									max_coordinate = max(ref_gene.tss + ref_gene.tts)
									MIN_coordinate = EXON_s[0]
									MAX_coordinate = EXON_e[-1]
									if (MIN_coordinate >= min_coordinate and MAX_coordinate <= max_coordinate):
										assoc.str_class = "novel_not_in_catalog"
										assoc.subtype = "not any annotated donor/acceptor"		

							elif assoc.str_class != "genic":
								assoc.str_class = "antisense" 
								assoc.genes = [] 
								ASgene = "novelGene_"+ref_gene.gene+"_AS"
								assoc.genes.append(ASgene)
							break 

					if assoc.str_class == "genic":
						break
				
				# If not overlapping with exons of genes, we see if it is found in the genomic region of a gene (overlapping introns) 

				if assoc.str_class == "intergenic":
					min_coordinate = min(ref_gene.tss + ref_gene.tts)
					max_coordinate = max(ref_gene.tss + ref_gene.tts)
					for l in range(len(EXON_s)):
						if (E_s >= min_coordinate and E_e <= max_coordinate) or (E_s <= min_coordinate and E_e >= max_coordinate) or (E_s<=min_coordinate<=E_e) or (E_s<=max_coordinate<=E_e):
							assoc.str_class = "genic_intron" #Not associated gene
							break
	return(assoc)
def single_exonAssociations(assoc, line_split, transcripts_gene):

	subject_transcript = line_split[0] 
	subject_gene = assoc.genes[0]
	assoc.subtype = "mono-exon"

	if subject_gene in transcripts_gene.keys():
		
		ref_gene = transcripts_gene[subject_gene]
		transcripts = ref_gene.transcripts
		exons_per_transcript = ref_gene.exonsTranscriptTog
		exons_per_transcript_no_all = ref_gene.exonsTranscript
		#we store the exon information
		EXON_s = [int(i)+1 for i in line_split[8].split(",")[0:-1]]
		EXON_e = [int(i) for i in line_split[9].split(",")[0:-1]]
	
		#si no solapa completamente con algun exon de sus gen, entonces tiene parte genomica.
		assoc.str_class = "genic"
		for i in range(len(exons_per_transcript)):
			overlap = False
			exons = exons_per_transcript[i]
			strand = ref_gene.strand
			for exon in exons:
				if (int(EXON_s[0]) in range(int(exon[0]),int(exon[1])+1)) and (int(EXON_e[0]) in range(int(exon[0]),int(exon[1])+1)):
					overlap = True
					break
			if overlap == True:
				transcript = transcripts[i]
				starts = []
				ends = []
				for j in exons:
					starts.append(j[0])
					ends.append(j[1])
				transcript_range = sorted(exons_together(starts, ends))

				if strand == "+":
					TSS_diff = transcript_range.index(int(EXON_s[0]))
					TTS_diff = len(transcript_range) - transcript_range.index(int(EXON_e[0]))

				if strand=="-":
					TTS_diff = transcript_range.index(int(EXON_s[0]))
					TSS_diff = len(transcript_range) - transcript_range.index(int(EXON_e[0]))

				if assoc.str_class != "incomplete-splice_match":
					assoc.modify(transcript, subject_gene, TSS_diff, TTS_diff, transcriptLength([int(i[0]) for i in exons],[int(i[1]) for i in exons]), len(exons))

				elif  (abs(TSS_diff) + abs(TTS_diff)) < assoc.get_total_diff():
					assoc.modify(transcript, subject_gene, TSS_diff, TTS_diff, transcriptLength([int(i[0]) for i in exons],[int(i[1]) for i in exons]), len(exons))

				else:
					assoc.transcripts.append(transcript)

				assoc.str_class = "incomplete-splice_match"
				assoc.subtype = "mono-exon"


		if assoc.str_class == "genic": # if not "incomplete-splice match" we check if not novel in catalog (example gene with intron retention in all their junctions)
			min_coordinate = min(ref_gene.tss + ref_gene.tts)
			max_coordinate = max(ref_gene.tss + ref_gene.tts)
			MIN_coordinate = EXON_s[0]
			MAX_coordinate = EXON_e[-1]
			if (MIN_coordinate >= min_coordinate and MAX_coordinate <= max_coordinate):
				assoc.str_class = "novel_in_catalog"		
				assoc.subtype = "mono-exon_by_intron_retention/s"


	return(assoc)
def clusteringNovelGeneIsoforms(dicc_novelGenes):

	dicc_novelGenes_sort = sorted(dicc_novelGenes, key = lambda tup:int(tup.coordToSort))
	novel_seen = {}

	index = 1
	for isoform in dicc_novelGenes_sort:
		delete = []
		found = False
		line_split = isoform.line #check if works in other way
		CHROM = line_split[1]
		STRAND = line_split[2]
		EXON_s = [int(i)+1 for i in line_split[8].split(",")[0:-1]]
		EXON_e = [int(i) for i in line_split[9].split(",")[0:-1]]
	
		for x in novel_seen:
			line_split2 = x.line
			chrom = line_split2[1]
			strand = line_split2[2]
			exon_s = [int(i)+1 for i in line_split2[8].split(",")[0:-1]]
			exon_e = [int(i) for i in line_split2[9].split(",")[0:-1]]
			if (EXON_s[0] > exon_e[-1]):
				delete.append(x)
				continue
			
			elif chrom==CHROM and strand==STRAND: 
				for l in range(len(EXON_s)):
					E_s = int(EXON_s[l])
					E_e = int(EXON_e[l])

					for p in range(len(exon_s)):
						e_s = int(exon_s[p])
						e_e = int(exon_e[p])
					
						if (E_s >= e_s and E_e <= e_e) or (E_s <= e_s and E_e >= e_e) or (E_s<=e_s<=E_e) or (E_s<=e_e<=E_e):

							n_g = novel_seen[x].assoc.genes
							for item in n_g:
								isoform.assoc.genes.append(item)

							found = True
							break
					
					if found == True:
						break

		if found == False:
			isoform.assoc.genes.append("novelGene_"+str(index))
			index +=1
			novel_seen[isoform] = isoform
		
		for k in delete:
			novel_seen.pop(k)

	return(dicc_novelGenes_sort)
def junctionInfo(assoc, line_split, transcripts_gene, junctions_list, args, indelInfo, fastaInf, covInf = None):

	sites = list(args.sites.split(",")) 

	genes = list(set(assoc.genes))
	
	query_transcript = line_split[0]
	strand = line_split[2]
	chrom = line_split[1]

	EXON_s = [int(i)+1 for i in line_split[8].split(",")[0:-1]]
	EXON_e = [int(i) for i in line_split[9].split(",")[0:-1]]

	JUNCTIONS_e = [int(i)-1 for i in EXON_s[1:] ]
	JUNCTIONS_s = [int(i)+1 for i in EXON_e[:-1]]

	JUNCTIONS = []
	[JUNCTIONS.append([JUNCTIONS_s[m],JUNCTIONS_e[m]]) for m in range(len(JUNCTIONS_s))]		

	JUNCTIONS_transcriptCoord = getTranscriptJunctCoordinates(EXON_s, EXON_e, assoc.strand)

	junctions = []
	for gene in genes:
		if gene in transcripts_gene:
			ref_gene = transcripts_gene[gene]
			[junctions.append(junction) for junction in ref_gene.junctions]

	s = [y[0] for y in junctions] 
	e = [y[1] for y in junctions]

	for junction_number in range(len(JUNCTIONS)):
		junct = start = end = "known"
		min_diff_s = 0
		min_diff_e = 0 
		site = "NA"
		canonical = "NA"
		tc = int(JUNCTIONS_transcriptCoord[0][junction_number])
		if len(junctions)==0:
			junct = start = end = "novel"
			min_diff_s = "NA"
			min_diff_e = "NA" 
		elif JUNCTIONS[junction_number] not in junctions:
			junct = "novel"
			a = JUNCTIONS[junction_number][0] not in s
			b = JUNCTIONS[junction_number][1] not in e
			if any([a,b]):
				if a==True:
					diff_s = [abs(JUNCTIONS[junction_number][0] - k) for k in s]
					min_diff_s = JUNCTIONS[junction_number][0] - s[diff_s.index(min(diff_s))]
					start = "novel"
				if b==True:
					diff_e = [abs(JUNCTIONS[junction_number][1] - k) for k in e]
					min_diff_e = e[diff_e.index(min(diff_e))] - JUNCTIONS[junction_number][1]
					end = "novel"

		#junctions_list.append(myQueryJunctions(query_transcript, "junction_"+str(junction_number), assoc.chrom, assoc.strand, junct, start, end, min_diff_s, min_diff_e, JUNCTIONS[junction_number][0], JUNCTIONS[junction_number][1], "u"))

		site = fastaInf[chrom][(JUNCTIONS[junction_number][0]-1):(JUNCTIONS[junction_number][0]+1)] + fastaInf[chrom][(JUNCTIONS[junction_number][1]-2):JUNCTIONS[junction_number][1]]
		if strand =="-":
			site = ReverseComplement(site)
		else:
			site = site.upper()
			if any([base not in ["A","T","C","G","N"] for base in site]):
				sys.stderr.write("Error: %s is NOT a DNA sequence in genome file" %(site))
				return None
				break

		if site in sites:
			canonical = "canonical"
		else:
			canonical= "non_canonical"

		coorString = ("-").join([query_transcript,str(JUNCTIONS[junction_number][0]), str(JUNCTIONS[junction_number][1])])

		indel = "NA"
		if indelInfo != None:
			if coorString in indelInfo:
				indel = "TRUE"
			else:
				indel= "FALSE"

		qj = myQueryJunctions(query_transcript, "junction_"+str(junction_number+1), assoc.chrom, assoc.strand, junct, start, end, min_diff_s, min_diff_e, JUNCTIONS[junction_number][0], JUNCTIONS[junction_number][1], tc, site, canonical, indel)	

		if covInf!=None:
			j = "_".join([chrom,str(JUNCTIONS[junction_number][0]),str(JUNCTIONS[junction_number][1])])
			for covfile in covInf:
				if j in covInf[covfile]:
					qj.coverage.append(covInf[covfile][j])
				else:
					qj.coverage.append("0")		

		junctions_list.append(qj)

	return(junctions_list)
def junctionInfo_index(assoc, line_split, transcripts_gene, junctions_list, args, indelInfo, fastaInf, covInf = None):

	sites = list(args.sites.split(",")) 

	genes = list(set(assoc.genes))
	
	query_transcript = line_split[0]
	strand = line_split[2]
	chrom = line_split[1]

	EXON_s = [int(i)+1 for i in line_split[8].split(",")[0:-1]]
	EXON_e = [int(i) for i in line_split[9].split(",")[0:-1]]

	JUNCTIONS_e = [int(i)-1 for i in EXON_s[1:] ]
	JUNCTIONS_s = [int(i)+1 for i in EXON_e[:-1]]

	JUNCTIONS = []
	[JUNCTIONS.append([JUNCTIONS_s[m],JUNCTIONS_e[m]]) for m in range(len(JUNCTIONS_s))]		

	JUNCTIONS_transcriptCoord = getTranscriptJunctCoordinates(EXON_s, EXON_e, assoc.strand)

	junctions = []
	for gene in genes:
		if gene in transcripts_gene:
			ref_gene = transcripts_gene[gene]
			[junctions.append(junction) for junction in ref_gene.junctions]

	s = [y[0] for y in junctions] 
	e = [y[1] for y in junctions]

	for junction_number in range(len(JUNCTIONS)):
		junct = start = end = "known"
		min_diff_s = 0
		min_diff_e = 0 
		site = "NA"
		canonical = "NA"
		tc = int(JUNCTIONS_transcriptCoord[0][junction_number])
		if len(junctions)==0:
			junct = start = end = "novel"
			min_diff_s = "NA"
			min_diff_e = "NA" 
		elif JUNCTIONS[junction_number] not in junctions:
			junct = "novel"
			a = JUNCTIONS[junction_number][0] not in s
			b = JUNCTIONS[junction_number][1] not in e
			if any([a,b]):
				if a==True:
					diff_s = [abs(JUNCTIONS[junction_number][0] - k) for k in s]
					min_diff_s = JUNCTIONS[junction_number][0] - s[diff_s.index(min(diff_s))]
					start = "novel"
				if b==True:
					diff_e = [abs(JUNCTIONS[junction_number][1] - k) for k in e]
					min_diff_e = e[diff_e.index(min(diff_e))] - JUNCTIONS[junction_number][1]
					end = "novel"

		#junctions_list.append(myQueryJunctions(query_transcript, "junction_"+str(junction_number), assoc.chrom, assoc.strand, junct, start, end, min_diff_s, min_diff_e, JUNCTIONS[junction_number][0], JUNCTIONS[junction_number][1], "u"))

		site = fastaInf[chrom].seq[(JUNCTIONS[junction_number][0]-1):(JUNCTIONS[junction_number][0]+1)] + fastaInf[chrom].seq[(JUNCTIONS[junction_number][1]-2):JUNCTIONS[junction_number][1]]
		if strand =="-":
			site = ReverseComplement(site)
		else:
			site = site.upper()
			if any([base not in ["A","T","C","G","N"] for base in site]):
				sys.stderr.write("Error: %s is NOT a DNA sequence in genome file" %(site))
				return None
				break

		if site in sites:
			canonical = "canonical"
		else:
			canonical= "non_canonical"

		coorString = ("-").join([query_transcript,str(JUNCTIONS[junction_number][0]), str(JUNCTIONS[junction_number][1])])

		indel = "FALSE"
		if coorString in indelInfo:
			indel = "TRUE"

		qj = myQueryJunctions(query_transcript, "junction_"+str(junction_number+1), assoc.chrom, assoc.strand, junct, start, end, min_diff_s, min_diff_e, JUNCTIONS[junction_number][0], JUNCTIONS[junction_number][1], tc, site, canonical, indel)	

		if covInf!=None:
			j = "_".join([chrom,str(JUNCTIONS[junction_number][0]),str(JUNCTIONS[junction_number][1])])
			for covfile in covInf:
				if j in covInf[covfile]:
					qj.coverage.append(covInf[covfile][j])
				else:
					qj.coverage.append("0")		

		junctions_list.append(qj)

	return(junctions_list)

def isoformClassification(args, myQueryTranscripts_list, transcripts_chrom_1exon, transcripts_chrom_exons, transcripts_gene, transcripts_gene_chrom, seq, indelsJunc): 

	## read coverage files if provided
 	coverage_header = ""
 	if args.coverage!=None:
 		cov = STARcov_parser(args.coverage)
 		coverage_header = "\t"+("\t").join(cov.keys())
 		sys.stdout.write("\n**** Reading Splice Junctions coverage files.\n")
 	else:
	 	cov = None
 		sys.stdout.write("\nSplice Junction Coverage files not provided.\n")


 	# running classification
 	sys.stdout.write("\n**** Performing Classification of Isoforms...\n")

 	index_genes = {}
 	index_chrom_exons = {}
 	index_chrom_1exon = {}

 	list_novelGenes = []

 	outputPathPrefix = args.dir+"/"+args.output

 	outputClassPath = outputPathPrefix+"_classification.txt"
 	output = open(outputClassPath+"_tmp", "w")

 	outputJuncPath = outputPathPrefix+"_junctions.txt"
 	output_junc = open(outputJuncPath+"_tmp", "w")
 	output_junc.write("isoform\tchrom\tstrand\tjunction_number\tgenomic_start_coord\tgenomic_end_coord\ttranscript_coord\tjunction_category\tstart_site_category\tend_site_category\tdiff_to_Ref_start_site\tdiff_to_Ref_end_site\tbite_junction\tsplice_site\tcanonical\tindel_near_junct\tsample_with_cov\ttotal_coverage"+coverage_header+"\n")

	isoforms_inf = {}

 	for element in myQueryTranscripts_list:
 		novel = 0
 		line = element.line
 		line_split = line.split("\t")
 		query_transcript = line_split[0]

 		junctions_list = [] 
 		over_nt = {}

 		# Association to known genes 
 		(assoc, junctions, index_chrom_exons, index_chrom_1exon) = transcriptsKnownSpliceSites(transcripts_chrom_1exon, transcripts_chrom_exons, line_split, index_chrom_exons, index_chrom_1exon, transcripts_gene, seq, args.window)
 		if assoc.str_class!="":
 			junctions_list = junctionInfo(assoc, line_split, transcripts_gene, junctions_list, args, indelsJunc, seq, covInf = cov)

 		# If isoform share any splice site with any gene (but not is full/incomplete splice match), we divide them into: Novel in catalog, novel not in catalog or Fusion gene.
 		if assoc.str_class == "anyKnownSpliceSite": 			 			
 			assoc = novelIsoformsKnownGenes(assoc, line_split, transcripts_gene)

 		# Association to genes by overlapping (not known junctions). Multi-exon genes divided into: antisense (if overlapping opposite strand), novel in catalog (if not same junctions but inside the bondaries of a known gene), genic_genomic (rest of them).
 		elif assoc.str_class == "":
 			assoc = associationOverlapping(assoc, line_split, transcripts_gene_chrom, index_genes, over_nt)
 			junctions_list = junctionInfo(assoc, line_split, transcripts_gene, junctions_list, args, indelsJunc, seq, covInf = cov)

 		# Finally we can divide monoexonic genic genes between "fragment_match", "complete_intron_retention" or "genic_genomic"
 		if assoc.str_class=="genic" and assoc.exon_number==1:
			assoc = single_exonAssociations(assoc, line_split, transcripts_gene)
		
		if assoc.str_class in ["intergenic", "genic_intron"]: # those are the ones that have not any associated gene so they can share novel genes. Clustering of those isoforms from unknown genes.
			list_novelGenes.append(myQueries(line_split, int(line.split("\t")[8].split(",")[0]), assoc))

		else:
			isoforms_inf[query_transcript] = assoc
			output.write(query_transcript+"\t")	
			output.write(str(assoc))
			output.write("\n")

		# Writing information for each transcript

		if len(junctions_list) > 0:
			for z in junctions_list:
				output_junc.write(str(z))
				output_junc.write("\n")

	# Clustering of novel isoforms from novel genes.
	novelGenes = clusteringNovelGeneIsoforms(list_novelGenes)

	#del transcripts_chrom_1exon, transcripts_chrom_exons, index_chrom_exons, index_chrom_1exon, transcripts_gene, list_novelGenes, junctions_list
	#gc.collect()

	for x in novelGenes:
		isoforms_inf[x.line[0]] = x.assoc

		output.write(x.line[0]+"\t")
		output.write(str(x.assoc))
		output.write("\n")

	output.close()		
	output_junc.close()

	#del novelGenes
	#gc.collect()

	return(isoforms_inf)
def FLcount_parser (fl_files):

	# may be one file, files separated by comma or a directory where file are.
	if os.path.isdir(fl_files)==True:
		cov_paths = [os.path.join(fl_files,fn) for fn in next(os.walk(fl_files))[2]]
	else: 
		cov_paths = fl_files.split(",")

	fl_dicc = {}
	for path in cov_paths:
		if not os.path.exists(path):
			sys.stdout.write("\nERROR: Unable to read %s FL count files.\n" % path)
		else:
			p = open(path.strip(), "r")
			for line in p:
				if line[0] != "#" and line[0:4] != "pbid":
					pbid = line.split("\t")[0]
					fl = int(line.split("\t")[1])
					if pbid in fl_dicc:
						fl_dicc[pbid] = fl_dicc[pbid] + fl
					else:
						fl_dicc[pbid] = fl
			p.close()

	return(fl_dicc)
def mean(data):
	"""Return the sample arithmetic mean of data."""
	n = len(data)
	if n < 1:
		raise ValueError('mean requires at least one data point')
	return sum(data)/n # in Python 2 use sum(data)/float(n)
def _ss(data):
	"""Return sum of square deviations of sequence data."""
	c = mean(data)
	ss = sum((x-c)**2 for x in data)
	return ss
def pstdev(data):
	"""Calculates the population standard deviation."""
	n = len(data)
	if n < 2:
		raise ValueError('variance requires at least two data points')
	ss = _ss(data)
	pvar = ss/n # the population variance
	return pvar**0.5
def run(args):

	start3 = timeit.default_timer()

 	sys.stdout.write("\n**** Parsing provided files...\n")

 	# parse genome
	# seq = SeqIO.index(args.genome,'fasta')
	# chrList=[]
	# [chrList.append(chr) for chr in seq]
	# print(chrList)
	seq = genome_parser(args.genome)

	## correction of sequences and ORF prediction (if gtf provided instead of fasta file, correction of sequences will be skipped)
	orfLenDicc = correctionPlusORFpred(args, seq.keys())

 	## parse reference transcripts
 	transcripts_chrom_1exon, transcripts_chrom_exons, transcripts_gene, transcripts_gene_chrom = reference_parser(args, seq.keys()) 

	## parse query isoforms  
	myQueryTranscripts_list = isoforms_parser(args)

 	## Run indel computation if sam exists
 	if os.path.exists(corrSAM):
		(indelsJunc, indelsTotal) = indels(corrSAM)
	else:
		indelsJunc = None
		indelsTotal = None

	# isoform classification + intra-priming + transcript and junction characterization
	isoforms_inf = isoformClassification(args, myQueryTranscripts_list, transcripts_chrom_1exon, transcripts_chrom_exons, transcripts_gene, transcripts_gene_chrom, seq, indelsJunc)
 	print("\nNumber of classified isoforms "+str(len(isoforms_inf)))
 	outputPathPrefix = args.dir+"/"+args.output
 	outputClassPath = outputPathPrefix+"_classification.txt"
 	outputJuncPath = outputPathPrefix+"_junctions.txt"

	## RT-switching computation
	sys.stdout.write("\n**** RT-switching computation...\n")
	indexRefGenome = args.genome+".fai"
	if not os.path.isfile(indexRefGenome):
		sys.stderr.write("ERROR: %s not find" %(indexRefGenome))					
		sys.exit()
	RTS_list, RTSjunc_list = rts([outputJuncPath+"_tmp", args.genome, indexRefGenome ,"-a"]) # list of RTS-transcripts
	for x in RTS_list:
		if x in isoforms_inf:
			#if isoforms_inf[x].str_class not in ["full-splice_match", "incomplete-splice_match"]:
			isoforms_inf[x].RT_switching = "TRUE"
	del RTS_list
	gc.collect()

 	# ORF information
 	if len(orfLenDicc)>0:
	 	for pbid in orfLenDicc:
	 		if pbid in isoforms_inf:
	 			isoforms_inf[pbid].coding = "coding"
	 			isoforms_inf[pbid].ORFlen = orfLenDicc[pbid].orf_length
	 			isoforms_inf[pbid].CDS_start = orfLenDicc[pbid].cds_start
	 			isoforms_inf[pbid].CDS_end = orfLenDicc[pbid].cds_end

 	del orfLenDicc
 	gc.collect()


	## FSM classification 
 	geneFSM_dicc = defaultdict(list)
 	for iso in isoforms_inf:
 		gene = isoforms_inf[iso].geneName()
 		geneFSM_dicc[gene].append(isoforms_inf[iso].str_class)

	## FL count files
 	if args.fl_count:
 		sys.stdout.write("\n**** Reading Full-length read abundance files.\n")
 		fl = FLcount_parser(args.fl_count)
 		for i in fl:
 			if i not in isoforms_inf:
 				sys.stderr.write("\nWARNING: %s ID found in FL abundance files but not in transcriptome.\n" %i)
 		if len(fl)>0:
 			for x in isoforms_inf:
	 			if x in fl:
	 				isoforms_inf[x].FL = fl[x]
	 			else:
	 				isoforms_inf[x].FL = 0
 	else:
 		sys.stdout.write("\nFull-length read abundance files not provided.\n")


	## Isoform expression information
	if args.expression:
		sys.stdout.write("\n**** Reading Isoform Expression Information.\n")
		exp_dicc = expression_parser(args.expression)
		geneExp_dicc = {}
		for iso in isoforms_inf:
			if iso not in exp_dicc:
				exp_dicc[iso] = 0
				sys.stdout.write("\nIsoform %s not found in expression matrix. Nule expression associated.\n" %(iso))
			gene = isoforms_inf[iso].geneName()
			if gene not in geneExp_dicc:
				geneExp_dicc[gene] = exp_dicc[iso]
			else:
				geneExp_dicc[gene] = geneExp_dicc[gene]+exp_dicc[iso]
	else:
		exp_dicc = None
		geneExp_dicc = None
		sys.stdout.write("\nIsoforms expression files not provided.\n")


 	## Adding indel, FSM class and expression information
 	for iso in isoforms_inf:
 		gene = isoforms_inf[iso].geneName()
 		if exp_dicc!=None and geneExp_dicc!=None:
 			isoforms_inf[iso].geneExp = geneExp_dicc[gene]
 			isoforms_inf[iso].isoExp = exp_dicc[iso]
 		if len(geneFSM_dicc[gene])==1:
 			isoforms_inf[iso].FSM_class = "A"
 		elif "full-splice_match" in geneFSM_dicc[gene]:
 			isoforms_inf[iso].FSM_class = "C"
 		else:
 			isoforms_inf[iso].FSM_class = "B"
 	
 	if indelsTotal!=None:
 		for iso in isoforms_inf:
 			if iso in indelsTotal:
 				isoforms_inf[iso].nIndels = indelsTotal[iso]
 			else:
 				isoforms_inf[iso].nIndels = 0
 				#sys.stdout.write("\nIsoform %s not found in sam file. Nule number of indels associated.\n" %(iso))


 	del  geneFSM_dicc, myQueryTranscripts_list, geneExp_dicc, exp_dicc
	gc.collect()

	## Read junction files and create attributes per transcript

	with open(outputJuncPath+"_tmp", "r") as inFile:
		header = inFile.readline()
		cov_dicc = defaultdict(list)
		for line in inFile:
			line = line.rstrip()
			pbid = line.split("\t")[0]
			can = line.split("\t")[14]
			indel =  line.split("\t")[15]
			bite =  line.split("\t")[12]
			sample_cov = line.split("\t")[16]
			cov = line.split("\t")[17]
			j_c = line.split("\t")[6]

			if cov!= "NA":
				cov_dicc[pbid].append(int(cov))
				cov = int(cov)
			if sample_cov!="NA":
				sample_cov = int(sample_cov)
			if isoforms_inf[pbid].min_cov == "NA" or isoforms_inf[pbid].min_cov > cov:
				isoforms_inf[pbid].min_cov = cov
				isoforms_inf[pbid].min_cov_pos = j_c
			if isoforms_inf[pbid].min_samp_cov == "NA" or isoforms_inf[pbid].min_samp_cov > sample_cov:
				isoforms_inf[pbid].min_samp_cov = sample_cov
			if can == "non_canonical":
				isoforms_inf[pbid].canonical = can
	 		if indel == "TRUE":
	 			if isoforms_inf[pbid].nIndelsJunc != "NA":
	 				isoforms_inf[pbid].nIndelsJunc = isoforms_inf[pbid].nIndelsJunc + 1	
	 			else:
	 				isoforms_inf[pbid].nIndelsJunc = 1
	 		else:
	 			isoforms_inf[pbid].nIndelsJunc = 0
	 		if bite == "TRUE":
	 			isoforms_inf[pbid].bite = "TRUE"
	for i in cov_dicc:
		if len(cov_dicc[i]) > 1:
			isoforms_inf[i].sd = pstdev(cov_dicc[i])


	#### Printing output file:

	sys.stdout.write("\n**** Writing output files...\n")

	with open(outputClassPath, "w") as newFile:

		newFile.write("isoform\tchrom\tstrand\tlength\texons\tstructural_category\tassociated_gene\tassociated_transcript\tref_length\tref_exons\tdiff_to_TSS\tdiff_to_TTS\tsubcategory\tRTS_stage\tall_canonical\tmin_sample_cov\tmin_cov\tmin_cov_pos\tsd_cov\tFL\tn_indels\tn_indels_junc\tbite\tiso_exp\tgene_exp\tratio_exp\tFSM_class\tcoding\tORF_length\tCDS_length\tCDS_start\tCDS_end\tperc_A_downstream_TTS\n")
		for iso in isoforms_inf:
			newFile.write(iso+"\t")
			newFile.write(str(isoforms_inf[iso]))
			newFile.write("\n")		

	with open(outputJuncPath, "w") as newFile:
 		with open(outputJuncPath+"_tmp", "r") as oldFile:
 				line = oldFile.readline()
 				split_l = (line.split("\t"))
 				l = len(split_l)
				newFile.write("\t".join((split_l[0:15]+["RTS_junction"]+split_l[15:l])))

				for line in oldFile:
	 				split_l = (line.split("\t"))
					rts_status = "FALSE"
					if ":".join([split_l[0], split_l[1],split_l[4],split_l[5]]) in RTSjunc_list: 
						rts_status = "TRUE"
					newFile.write( "\t".join( (split_l[0:15]+[rts_status]+split_l[15:l]) ) )


	del isoforms_inf
	del RTSjunc_list
	gc.collect()

	os.remove(outputClassPath+"_tmp")
	os.remove(outputJuncPath+"_tmp")
	os.remove(queryFile)
	os.remove(referenceFiles)


	## Generating report

	sys.stdout.write("\n**** Generating SQANTI report...\n")
	rscriptPath = distutils.spawn.find_executable('Rscript')
	subprocess.call ([rscriptPath, utilitiesPath+"SQANTI_report.R", outputClassPath, outputJuncPath])
	stop3 = timeit.default_timer()


	print stop3 - start3

def main():

	global utilitiesPath

	#arguments
	parser = argparse.ArgumentParser(description="Structural and Quality Annotation of Novel Transcript Isoforms")
	parser.add_argument('isoforms', help='\tIsoforms (Fasta or gtf format; By default "fasta". GTF if specified -g option)') 
	parser.add_argument('annotation', help='\t\tReference annotation file (GTF format)')
	parser.add_argument('genome', help='\t\tReference genome (Fasta format)')
	parser.add_argument('-g', '--gtf', help='\t\tUse when running Sqanti by using as input a gtf of isoforms', action='store_true')
	parser.add_argument('-e','--expression', help='\t\tExpression matrix', required=False)
	parser.add_argument('-x','--gmap_index', help='\t\tPath and prefix of the reference index created by gmap_build. Mandatory unless -m option is specified.')
	parser.add_argument('-o','--output', help='\t\tPrefix for output files.', required=False)
	parser.add_argument('-d','--dir', help='\t\tDirectory for output files. Default: Directory where the script was run.', required=False)
	parser.add_argument('-c','--coverage', help='\t\tJunction coverage files (comma-separated list of SJ.out.tab files generated by STAR or directory where they are in).', required=False)
	parser.add_argument('-s','--sites', default="ATAC,GCAG,GTAG", help='\t\tSet of splice sites to be considered as canonical (comma-separated list of splice sites). Default: GTAG,GCAG,ATAC.', required=False)
	parser.add_argument('-w','--window', default="20", help='\t\tSize of the window in the genomic DNA screened for Adenine content downstream of TTS', required=False, type=int)
	parser.add_argument('-n','--name', help='\t\tUse gene_name tag from GTF to define genes. Default: gene_id used to define genes', action='store_true')
	parser.add_argument('-fl', '--fl_count', help='\t\tFull-length PacBio abundance files (comma-separated list of PacBio abundance files generated by PacBio or directory where there are in).', required=False)
	parser.add_argument("-v", "--version", help="Display program version number.", action='version', version='SQANTI 1.2')

	args = parser.parse_args()

	# path and prefix for output files
	if args.output==None:
		args.output = os.path.splitext(os.path.basename(args.isoforms))[0]

	if args.dir==None:
		args.dir = os.getcwd() 
	else:
		if (os.path.isdir(args.dir)==False):
			sys.stderr.write("ERROR: '%s' directory doesn't exist\n" %(args.dir))
			sys.exit()
		else:
			args.dir = os.path.abspath(args.dir)


	args.genome = os.path.abspath(args.genome)
	if not os.path.isfile(args.genome):
		sys.stderr.write("ERROR: '%s' genome doesn't exist\n" %(args.genome))
		sys.exit()


	args.isoforms = os.path.abspath(args.isoforms)
	if not os.path.isfile(args.isoforms):
		sys.stderr.write("ERROR: '%s' doesn't exist\n" %(args.isoforms))
		sys.exit()

	args.annotation = os.path.abspath(args.annotation)
	if not os.path.isfile(args.annotation):
		sys.stderr.write("ERROR: '%s' doesn't exist\n" %(args.annotation))
		sys.exit()

	if not args.gtf and not args.gmap_index:
		sys.stderr.write("ERROR: for the alignment of isoform sequences a gmap_index must be specified.\n")
		sys.exit()	


	# Running functionality
	sys.stdout.write("\n**** Running SQANTI...\n")

	run(args)	


if __name__ == "__main__":
	main()
