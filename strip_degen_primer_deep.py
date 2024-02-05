#!/usr/bin/python

__author__ = "Patrick Denis Browne"
__email__ = "pdbr@plen.ku.dk"
__date__ = "02/2024"
__license__ = "GPLv3"

degen_table = {'A':['A'], 'T':['T'], 'G':['G'], 'C':['C']}
degen_table['R'] = ['A', 'G']
degen_table['Y'] = ['C', 'T']
degen_table['M'] = ['A', 'C']
degen_table['K'] = ['G', 'T']
degen_table['S'] = ['C', 'G']
degen_table['W'] = ['A', 'T']
degen_table['H'] = ['A', 'C', 'T']
degen_table['B'] = ['C', 'G', 'T']
degen_table['V'] = ['A', 'C', 'G']
degen_table['D'] = ['A', 'G', 'T']
degen_table['N'] = ['A', 'C', 'G', 'T']

comp_nucls = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
comp_nucls['R'] = 'Y'
comp_nucls['Y'] = 'R'
comp_nucls['M'] = 'M'
comp_nucls['K'] = 'K'
comp_nucls['S'] = 'S'
comp_nucls['W'] = 'W'
comp_nucls['H'] = 'D'
comp_nucls['B'] = 'V'
comp_nucls['V'] = 'B'
comp_nucls['D'] = 'H'
comp_nucls['N'] = 'N'


def str_to_list(nuclseq):
    return [nucl.upper() for nucl in nuclseq]

def check_score(query, subject, maxmismatch=3):
    qnucls = str_to_list(query)
    snucls = str_to_list(subject)
    match = 0
    mismatch = 0
    for i,qnucl in enumerate(qnucls):
        snucl = snucls[i]
        if snucl in degen_table[qnucl]:
            match += 1
        else:
            mismatch += 1
        if mismatch > maxmismatch:
            return None, None
    return match, mismatch

def reverse_complement(seqstr):
    nucls = [comp_nucls[nucl] for nucl in seqstr]
    nucls.reverse()
    revc = ''.join(nucls)
    return revc

def search_forward(query, subject, depth=0, maxmismatch=3):
    best_score = 0
    best_index = None
    for i in range(depth+1):
        subject_subseq = subject[i : i+len(query)]
        match, mismatch = check_score(query, subject_subseq, maxmismatch=maxmismatch)
        if isinstance(match, int):
            if match > best_score:
                best_score = match
                best_index = i
                if best_score == len(query):
                    break
    return best_index

def search_reverse(query, subject, depth=0, maxmismatch=3):
    revcsubject = reverse_complement(subject)
    best_index = search_forward(query, revcsubject, depth=depth, maxmismatch=maxmismatch)
    return best_index

def find_barcode(query, barcodes, maxmismatch=2):
    bc = None
    for barcode in barcodes:
        i, j = check_score(query, barcode)
        if isinstance(i, int):
            bc = barcode
            break
    return bc

if __name__ == '__main__':
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    import sys
    import gzip
    import itertools
    import argparse

    parser = argparse.ArgumentParser(description='A tool to remove artefactual sequences from reads')
    parser.add_argument('-i', '--infile', default=sys.stdin,
                        help="Input file. Fastq format only. Optionally gzip compressed. Defaults to reading from STDIN")
    parser.add_argument('-o', '--outfile', default=sys.stdout,
                        help='Output file name. Will be written as uncompressed fastq. Default is STDOUT')
    parser.add_argument('-f', '--forward', required=True,
                        help='Forward primer/adapter seqeunce. IUPAC degeneracies supported.')
    parser.add_argument('-r', '--reverse', required=True,
                        help='Reverse primer/adapter sequence. Reverse complement orientation to "forward". IUPAC supported')
    parser.add_argument('-d', '--depth', default=100, type=int,
                        help='Positive integer. Searches <INT> bases from the ends for primers/adapters. \
                            Higher values may drastically reduce runtimes. Use a low value if the primer/adapter \
                                is generally found near the extremes of the read.')
    parser.add_argument('-m', '--max_mismatches', default=0, type=int,
                        help='Maximum number of mismatches to either primer/adapter. Default: 0')
    parser.add_argument('--fasta', action="store_true",
                        help='Use fastA format. Default is to use fastQ format.')
    
    args = parser.parse_args()
    infile = args.infile
    outfile = args.outfile
    fwd_construct = args.forward.upper()
    rev_construct = args.reverse.upper()
    DEPTH = args.depth
    maxmismatches = args.max_mismatches
    infmt = 'fastq'
    outfmt = 'fastq'
    if args.fasta:
        infmt = 'fasta'
        outfmt = 'fasta'

    try:
        if infile.endswith('.gz'):
            fh_in = gzip.open(infile, 'rt')
        else:
            fh_in = open(infile, 'rt')
    except AttributeError:
        fh_in = infile
    
    outrecs = []
    numrecs = 0
    
    for rec in SeqIO.parse(fh_in, infmt):
        numrecs += 1
        seqstr = str(rec.seq).upper()
        fwd_i = search_forward(fwd_construct, seqstr, depth=DEPTH, maxmismatch=maxmismatches)
        rev_i = search_reverse(rev_construct, seqstr, depth=DEPTH, maxmismatch=maxmismatches)
        revc_fwd_i = search_reverse(fwd_construct, seqstr, depth=DEPTH, maxmismatch=maxmismatches)
        revc_rev_i = search_forward(rev_construct, seqstr, depth=DEPTH, maxmismatch=maxmismatches)
        if isinstance(fwd_i, int) and isinstance(rev_i, int):
            idx_start = fwd_i + len(fwd_construct)
            idx_end = len(rec) - ( rev_i + len(rev_construct) )
            newrec = rec[idx_start : idx_end]
            outrecs.append(newrec)
        elif isinstance(revc_fwd_i, int) and isinstance(revc_rev_i, int):
            idx_start = revc_rev_i + len(rev_construct)
            idx_end = len(rec) - ( revc_fwd_i + len(fwd_construct) )
            newrec = rec[idx_start : idx_end]
            newrec.seq = newrec.seq.reverse_complement()
            for k in newrec.letter_annotations.keys():
                newrec.letter_annotations[k].reverse()
            outrecs.append(newrec)
    
    numfilt = SeqIO.write(outrecs, outfile, outfmt)
    
    print(f"Processed {numrecs} records. {numfilt} records written to output.", file=sys.stderr)


