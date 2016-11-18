#!/usr/bin/env python

from __future__ import division
from __future__ import print_function

import sys, os
import subprocess
import argparse
import time
from collections import namedtuple

# Written by Hector del Risco - hdelrisco@ufl.edu

#### Function to check splice junctions for possible RT switching

def rts(args):

    def get_args(args):
        # parse command line arguments
        parser = argparse.ArgumentParser(description="Check splice junctions for possible RT switching")
        parser.add_argument('sjFilepath', type=str, help='file with splice junction information')
        parser.add_argument('mmfaFilepath', type=str, help='path to reference genome')
        parser.add_argument('mmfaiFilepath', type=str, help='path to reference genome index')
        parser.add_argument("-m", "--min_match", type=int, default=8, choices=range(4, 11), help="Minimum number of bases required to match. Default: 8")
        parser.add_argument("-a", "--allow_mismatch", default=False, action="store_true", help="Specify to allow 1 base mismatch in sequences (indels are not allowed)")
        parser.add_argument("-w", "--wiggle_count", type=int, default=1, choices=range(0, 4), help="Number of bases allowed to wiggle on each side of ideal RTS sequence location. Default: 1")
        parser.add_argument("-t", "--include_type", default='a', choices=['a', 'c', 'n'], help="Type of splice junctions to include (a for all, c for canonical, and n for non-canonical). Default: a")
        parser.add_argument("-c", "--include_category", default='a', choices=['a', 'n', 'k'], help="Category of splice junctions to include (a for all, n for novel, and k for known). Default: a")
        parser.add_argument("-v", "--version", help="Display program version number", action='version', version='%(prog)s 0.1')
        return parser.parse_args(args)

    # Load reference genome chromosome index
    # It returns a hash of all the entries with an array containing [fileOffset, length, baseCnt, charCnt]
    def loadRefGenomeIdx(filename):

        errflg = False
        chrIndex = {}
        lines = []
        tfunc = time.time()
        processed = 0
        
        # open input file and read all data into memory
        try:
            tstart = time.time()
            #print("Loading reference genome index from: {}".format(filename))
            f = open(filename, 'r')
        except IOError as e:
            sys.stderr.write("ERROR - Unable to open file, '{0}': {1}\n".format(filename, e))
            errflg = True
        else:
            try:
                lines = f.read().splitlines()
                f.close()
                lenlines = len(lines)
                #print("Loaded {} lines in {} seconds.".format(lenlines, time.time() - tstart))
                idx = 0
                for line in lines:
                    if idx > 0 and line[0] != "#":
                        fields = line.split("\t")
                        if len(fields) == 5:
                            # [offset, length, baseCnt, charCnt]
                            chrIndex[fields[0]] = [int(fields[2]), int(fields[1]), int(fields[3]), int(fields[4])]
                            processed += 1
                        else:
                            sys.stderr.write("Invalid line ({0}): {1}\n".format(idx, line))
                            errflg = True
                            break
                    idx += 1
            except IOError as e:
                sys.stderr.write("ERROR - Unable to read file, '{0}': {1}\n".format(filename, e))
                errflg = True

        if errflg:
            raise SystemExit(1)
        #print("Processed {} records.\nFunction elapsed time: {} seconds.".format(processed, time.time() - tfunc))
        return chrIndex

    #
    # Load splice junctions
    # Expected file format:
    #
    # transcript    junctionNumber  chrom   strand  genomicStartCoord   genomicEndCoord transcriptCoord ...
    #
    # It returns a hash of all the transcripts with a corresponding array of SpliceJunctions
    # (namedtuple is less efficient but it makes the code easier to read - can change to plain array for performance)
    # it also returns SJCounts
    #
    SpliceJunctions = namedtuple("SpliceJunctions", "trans, sjn, chromo, strand, strpos, endpos, transpos, category, startCat, endCat, type")
    SJCounts = namedtuple("SJCounts", "trans, sjTotal, sj, knownCanonical, knownNonCanonical, novelCanonical, novelNonCanonical")

    def loadSpliceJunctions(filepath):
        errflg = False
        index = {}
        sjUnique = {}
        lines = []
        sjTotalCnt = 0
        kcCnt = 0
        kncCnt = 0
        ncCnt = 0
        nncCnt = 0
        tfunc = time.time()

        try:
            # open input file
            tstart = time.time()
            #print("Loading splice junctions from: {}".format(filepath))
            f = open(filepath, 'r')
        except IOError as e:
            sys.stderr.write("ERROR - Unable to open file, '{0}': {1}\n".format(filepath, e))
            errflg = True
        else:
            try:
                # read all file data
                lines = f.read().splitlines()
                f.close()
                lenlines = len(lines)
                #print("Loaded {} lines in {} seconds.".format(lenlines, time.time() - tstart))

                # process all entries, line by line
                idx = 0
                for line in lines:
                    # skip past header and comment lines
                    if idx > 0 and line[0] != "#":
                        fields = line.split("\t")
                        if len(fields) >= 15:
                            # add transcript to hash even if it ends up having no junctions (all dups)
                            trans = fields[0]
                            if trans not in index:
                                index[trans] = []
                                
                            # check for unique splice junctions, does not make sense to use duplicates
                            sjUName = fields[2] + fields[3] + fields[4] + fields[5]
                            if sjUName not in sjUnique:
                                sjUnique[sjUName] = 1
                                index[trans].append(SpliceJunctions(trans, fields[1], fields[2], fields[3], int(fields[4]), int(fields[5]), int(fields[6]), fields[7], fields[8], fields[9], fields[14]))
                                if fields[7] == "known":
                                    if fields[14] == "canonical":
                                        kcCnt += 1
                                    else:
                                        kncCnt += 1
                                else:
                                    if fields[14] == "canonical":
                                        ncCnt += 1
                                    else:
                                        nncCnt += 1
                            else:
                                sjUnique[sjUName] = sjUnique[sjUName] + 1
                            sjTotalCnt += 1
                        else:
                            sys.stderr.write("Invalid splice junction line ({0}): {1}\n".format(idx, line))
                            errflg = True
                            break
                    idx += 1
            except IOError as e:
                sys.stderr.write("ERROR - Unable to read file, '{0}': {1}\n".format(filepath, e))
                errflg = True

        sjCounts = SJCounts(len(index), sjTotalCnt, len(sjUnique), kcCnt, kncCnt, ncCnt, nncCnt)
        if errflg:
            raise SystemExit(1)
        sjucnt = len(sjUnique)
        # print("Found {} transcripts with a total of {} splice junctions. Function elapsed time: {} seconds.".format(sjCounts.trans, sjCounts.sj, time.time() - tfunc))
        # print("Abbrev:\tTrans - Transcripts, SJ - splice junction\n\t\tKC - known canonical, KNC - known non-canonical, NC - novel canonical, NNC - novel non-canonical")
        # print("Fields:\tTrans\tTotalSJ\tSJ\tKC\tKNC\tNC\tNNC")
        # print("Totals:\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(sjCounts.trans, sjCounts.sjTotal, sjCounts.sj, sjCounts.knownCanonical, sjCounts.knownNonCanonical, sjCounts.novelCanonical, sjCounts.novelNonCanonical))
        # print("%:\t\t\t\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}".format(pct(sjCounts.knownCanonical, sjCounts.sj), pct(sjCounts.knownNonCanonical, sjCounts.sj), pct(sjCounts.novelCanonical, sjCounts.sj), pct(sjCounts.novelNonCanonical, sjCounts.sj)))
        return index, sjCounts

    #
    # ***********************************************************
    # All code above was strictly for loading existing data
    # Below is the code to do the actual RTS detection
    # ***********************************************************
    #
    # Visualization of where the pattern repeat is searched for using a wiggle (w) of 1
    #
    #       -------------------------                     ---------------------
    #                       w.....cba|w ------- w.....cba|w
    #       -------------------------                     ----------------------
    #                5' exon                intron              3' exon
    #
    # A perfect match is from the end (3') of the first exon, to the end (3') of the first intron
    # This applies to both the + and - strand (- strand is modified by the read function for 5' to '3 and complement)
    # In the case above a perfect match is "abc....." matching "abc....."
    #

    #
    # Check splice junctions for RTS and write results to results file
    #
    def checkSJforRTS(sjIdx, sjCounts, chrIdx, args, filename):
        rtsTrancList = []
        transcnt = 0
        transrts = 0
        sjCnt = 0
        kcCnt = 0
        kncCnt = 0
        ncCnt = 0
        nncCnt = 0
        exonSeq = ""
        intronSeq = ""
        wiggle = args.wiggle_count
        mmfaFilepath = args.mmfaFilepath
        cnt = PATSEQLEN + (2 * wiggle)

        tfunc = time.time()
        with open(filename, 'w') as f:
            f.write("trans\tsjNumber\tchromo\tstrand\tstrpos\tendpos\tcategory\ttype\texonSeq\tintronSeq\tmatchLen\tmatchPat\texonWiggle\tintronWiggle\tmismatch\n")
            for trans in sjIdx:
                transcnt += 1
                transflg = False
                # process all splice junctions
                for sj in sjIdx[trans]:
                    if args.include_type != 'a':
                        if args.include_type == 'c' and sj.type != "canonical":
                            continue
                        elif args.include_type == 'n' and sj.type != "non_canonical":
                            continue
                    if args.include_category != 'a':
                        if args.include_category == 'n' and sj.category != "novel":
                            continue
                        elif args.include_category == 'k' and sj.category != "known":
                            continue

                    # get sequences for pattern and search area
                    # the SJ start and end position are positioned at the start/end of the intron
                    # ths SJ start is always a lower position than the end regardless of strand
                    if sj.strand == "+":
                        # we always subtract to get to starting position
                        # the count includes 2 wiggles, both ends, so we adjust by 1 wiggle when positioning
                        # sequence data on disk: lowpos ----> hipos
                        # 5' -----exonSeq(SJstrpos)--------intronSeq(SJendpos) 3'
                        exonSeq = getRefSequence(chrIdx, sj.chromo, sj.strand, (sj.strpos - cnt + wiggle), cnt, mmfaFilepath)
                        intronSeq = getRefSequence(chrIdx, sj.chromo, sj.strand, (sj.endpos - cnt + wiggle + 1), cnt, mmfaFilepath)
                    else:
                        # we are almost on the starting position so just a minor adjustment
                        # sequence data on disk: lowpos ----> hipos
                        # 3' -----(SJstrpos)intronSeq--------(SJendpos)exonSeq 5'
                        intronSeq = getRefSequence(chrIdx, sj.chromo, sj.strand, (sj.strpos - wiggle), cnt, mmfaFilepath)
                        exonSeq = getRefSequence(chrIdx, sj.chromo, sj.strand, (sj.endpos - wiggle + 1), cnt, mmfaFilepath)
                        #if trans == "PB.308.3" and sj.sjn == "junction_22":
                        #    print("Sequences[{}-{}]: {}, {}]".format(sj.strpos, sj.endpos, intronSeq, exonSeq))

                    # check for RTS repeats and save results to file
                    if len(exonSeq) > 0 and len(intronSeq) > 0:
                        rtsflg, matchLen, matchPat, exwiggle, inwiggle, mismatch = checkForRepeatPat(trans, sj.sjn, exonSeq, intronSeq, args)
                        if rtsflg:
                            if not transflg:
                                transflg = True
                                transrts += 1
                            sjCnt += 1
                            if sj.category == "known":
                                if sj.type == "canonical":
                                    kcCnt += 1
                                else:
                                    kncCnt += 1
                            else:
                                if sj.type == "canonical":
                                    ncCnt += 1
                                else:
                                    nncCnt += 1
                            rtsTrancList.append(trans)
                            f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(trans, sj.sjn, sj.chromo, sj.strand, sj.strpos, sj.endpos, sj.category, sj.type, exonSeq, intronSeq, matchLen, matchPat, exwiggle, inwiggle, mismatch))

            sjRTSCounts = SJCounts(transrts, sjCnt, sjCnt, kcCnt, kncCnt, ncCnt, nncCnt)
            # trans values are not valid since only those with unique SJs are included
            #print("Found {} transcripts with a total of {} splice junctions. Function elapsed time: {} seconds.".format(sjRTSCounts.trans, sjRTSCounts.sj, time.time() - tfunc))
            #print("Fields:\tTrans\tSJ\tKC\tKNC\tNC\tNNC")
            #print("RTS:\t{}\t{}\t{}\t{}\t{}\t{}".format(sjRTSCounts.trans, sjRTSCounts.sj, sjRTSCounts.knownCanonical, sjRTSCounts.knownNonCanonical, sjRTSCounts.novelCanonical, sjRTSCounts.novelNonCanonical))
            #print("%:\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}".format(pct(sjRTSCounts.trans, sjCounts.trans), pct(sjRTSCounts.sj, sjCounts.sj), pct(sjRTSCounts.knownCanonical, sjCounts.knownCanonical), pct(sjRTSCounts.knownNonCanonical, sjCounts.knownNonCanonical), pct(sjRTSCounts.novelCanonical, sjCounts.novelCanonical), pct(sjRTSCounts.novelNonCanonical, sjCounts.novelNonCanonical)))
            # print("Found a total of {} RTS splice junctions. Function elapsed time: {} seconds.".format(sjRTSCounts.sj, time.time() - tfunc))
            # print("Fields:\tSJ\tKC\tKNC\tNC\tNNC")
            # print("RTS:\t{}\t{}\t{}\t{}\t{}".format(sjRTSCounts.sj, sjRTSCounts.knownCanonical, sjRTSCounts.knownNonCanonical, sjRTSCounts.novelCanonical, sjRTSCounts.novelNonCanonical))
            # print("%:\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}".format(pct(sjRTSCounts.sj, sjCounts.sj), pct(sjRTSCounts.knownCanonical, sjCounts.knownCanonical), pct(sjRTSCounts.knownNonCanonical, sjCounts.knownNonCanonical), pct(sjRTSCounts.novelCanonical, sjCounts.novelCanonical), pct(sjRTSCounts.novelNonCanonical, sjCounts.novelNonCanonical)))
        return rtsTrancList,sjRTSCounts
        
    def pct(value, total):
        return (value / total * 100)

    #
    # Check for possible RTS
    #
    # Look for the minimum pattern length match required using the following priorities:
    #
    #   1. try to find it without any mismatches and without any wiggles
    #   2. try to find it with wiggle
    #   3. try to find it with mismatches but without wiggle
    #   4. try to find it with mismatches and wiggle
    #
    # The reason for doing this is that you can end up having sequences matches
    # in different locations and of different lengths when you allow for wiggle and mismatches
    # and we do not want to set the wiggle count or mismatch flags unless needed for matching
    #
    def checkForRepeatPat(trans, sjnum, exonSeq, intronSeq, args):
        # look for pattern matches
        matchlen = 0
        matchpat = ""
        rtsflg = False
        exwiggle = 0
        inwiggle = 0
        mismatch = False
        maxloop = 1
        if args.allow_mismatch:
            maxloop = 2
        exseqlen = len(exonSeq)
        exseqend = exseqlen - args.wiggle_count - 1
        inseqlen = len(intronSeq)
        inseqend = inseqlen - args.wiggle_count - 1
        for miscnt in range(0, maxloop):
            if miscnt > 0:
                mismatch = True
            for ewiggle in range(0, args.wiggle_count + 1):
                for iwiggle in range(0, args.wiggle_count + 1):
                    # use 0 or positive ewiggle
                    expat = exonSeq[exseqend - args.min_match + ewiggle + 1:exseqend + ewiggle + 1]
                    # check vs the 0 or positive iwiggle
                    inpat = intronSeq[inseqend - args.min_match + iwiggle + 1:inseqend + iwiggle + 1]
                    rtsflg = isMatch(expat, inpat, mismatch)
                    if not rtsflg and iwiggle != 0:
                        # check vs the negative iwiggle
                        inpat = intronSeq[inseqend - args.min_match - iwiggle + 1:inseqend - iwiggle + 1]
                        rtsflg = isMatch(expat, inpat, mismatch)
                        if rtsflg:
                            iwiggle = -iwiggle

                    if not rtsflg and ewiggle != 0:
                        # use negative ewiggle
                        expat = exonSeq[exseqend - args.min_match - ewiggle + 1:exseqend - ewiggle + 1]
                        # check vs the 0 or positive iwiggle
                        inpat = intronSeq[inseqend - args.min_match + iwiggle + 1:inseqend + iwiggle + 1]
                        rtsflg = isMatch(expat, inpat, mismatch)
                        if rtsflg:
                            ewiggle = -ewiggle
                        elif iwiggle != 0:
                            # check vs the negative iwiggle
                            inpat = intronSeq[inseqend - args.min_match - iwiggle + 1:inseqend - iwiggle + 1]
                            rtsflg = isMatch(expat, inpat, mismatch)
                            if rtsflg:
                                ewiggle = -ewiggle
                                iwiggle = -iwiggle

                    if rtsflg:
                        matchlen = args.min_match
                        matchpat = expat
                        exwiggle = ewiggle
                        inwiggle = iwiggle
                        #if trans == "PB.5159.1" and sjnum == "junction_1":
                        #    print("Match found[{}, {}-{}; {}, {}-{}] miscnt({}): '{}'".format(ewiggle, (exseqend - args.min_match + ewiggle + 1), (exseqend + ewiggle + 1), iwiggle, (inseqend - args.min_match + iwiggle + 1), (inseqend + iwiggle + 1), miscnt, matchpat))
                        break
                if rtsflg:
                    break
            if rtsflg:
                break
        return rtsflg, matchlen, matchpat, exwiggle, inwiggle, mismatch

    #
    # Check if sequences match - sequences must have the same length
    #
    # Note: If mismatch flag is set, will allow 1 mismatch in comparison
    #       Regardless of mismacth flag value, indels are not allowed
    #
    def isMatch(exseq, inseq, allowMismatch):
        rtsflg = False
        
        if exseq == inseq:
            rtsflg = True
        elif allowMismatch:
            exseqlen = len(exseq)
            inseqlen = len(inseq)
            if exseqlen == inseqlen:
                cntmatch = 0
                mismatch = False
                for idx in range(0, exseqlen):
                    if exseq[idx:idx+1] == inseq[idx:idx+1]:
                        cntmatch += 1
                    else:
                        if mismatch:
                            break
                        else:
                            mismatch = True
                if cntmatch >= (exseqlen - 1):
                    rtsflg = True
        return rtsflg

    #
    # Get reference sequence for given genomic position
    #
    # Note: Negative strand sequences will be processed as follows:
    #       * the complement values will be returned, e.g TAC to ATG
    #       * the values will be reversed for 5' to 3' direction, e.g. ATG to GTA
    #         (this is needed to ensure the pattern match is checked from the 3' end for both strands)
    #
    def getRefSequence(chrIdx, chrName, strand, strpos, cnt, mmfaFilepath):
        seq = ""

        # check if we have index for given chromosome
        if chrName in chrIdx:
            # [offset, length, ntCount, charCnt]
            data = chrIdx[chrName]
            offset = data[0]
            basecnt = data[2]
            chrcnt = data[3]
            
            # read reference sequence data
            seq = readRefSeqData(strand, strpos, cnt, offset, basecnt, chrcnt, mmfaFilepath)
            if strand == "-":
                # reverse direction to show 5' to 3', seq already has the complements
                seq = seq[::-1]
        return seq

    #
    # Read reference sequence data
    #
    # Note: Negative strand sequences will be processed as follows:
    #       * the complement values will be returned, e.g TAC to ATG
    #
    def readRefSeqData(strand, strpos, cnt, offset, basecnt, chrcnt, mmfaFilepath):
        seq = ""
        
        with open(mmfaFilepath, 'r') as f:
            strline = int((strpos - 1) / basecnt)
            seekpos = int(offset + strline * chrcnt + ((strpos - 1) % basecnt))
            size = cnt + int(cnt/basecnt) + 4
            #print("refRead({}): str: {}, cnt: {}, off: {}, basecnt: {}, chrcnt: {}".format(strand, strpos, cnt, offset, basecnt, chrcnt))
            #print("strline: {}, seekpos: {}, size: {}".format(strline, seekpos, size))
            f.seek(seekpos)
            strdata = f.read(size)
            strdata = strdata.replace('\n', '')
            strdata = strdata.upper()
            seq = list(strdata[0:cnt])
            if strand == "-":
                # we need to get the complements
                idx = 0
                for idx in range(0, len(seq)):
                    if seq[idx] == 'A':
                        seq[idx] = 'T'
                    elif seq[idx] == 'C':
                        seq[idx] = 'G'
                    elif seq[idx] == 'T':
                        seq[idx] = 'A'
                    elif seq[idx] == 'G':
                        seq[idx] = 'C'
        return ''.join(seq)



    # Process command line request
    # get command line arguments
    args = get_args(args)
    # Output file
    rts_dir = "./RTS"
    if not os.path.exists(rts_dir):
        os.makedirs(rts_dir)

    rtsResultsFilepath = rts_dir+"/sj.rts.results.tsv"

    # Definitions
    PATSEQLEN = 10    # max size of sequence area to search for match (does not include wiggle)

    # load required data

    chrIdx = loadRefGenomeIdx(args.mmfaiFilepath)
    sjIdx, sjCounts = loadSpliceJunctions(args.sjFilepath)

    # # perform RTS analysis
    rtsTrancList, sjRTSCounts = checkSJforRTS(sjIdx, sjCounts, chrIdx, args, rtsResultsFilepath)

    return(set(rtsTrancList))


if __name__ == "__main__":
    import sys
    rts(sys.argv[1:])
