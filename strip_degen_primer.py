#! /usr/bin/python

from Bio import SeqIO
import gzip
import sys

def expand_degenerate(primer):
    """
    Takes a degenerate primer sequence and returns all possible 'A','T',
    'G' and 'C' only permutations of that degerate primer as a list.
    """
    degen_table = {}
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
    
    bc = primer.upper()
    bcs = []
    if bc[0] in degen_table:
        for base in degen_table[bc[0]]:
            bcs.append(base)
    else:
        bcs.append(bc[0])
    for i in range(1, len(bc)):
        if not bc[i] in degen_table:
            for j in range(len(bcs)):
                bcs[j] = bcs[j] + bc[i]
        if bc[i] in degen_table:
            appends = degen_table[bc[i]]
            original_length = len(bcs)
            bcs = bcs * len(appends)
            for j in range(len(bcs)):
                degen_index = int(j // original_length)
                append = appends[degen_index]
                bcs[j] = bcs[j] + append
    return(bcs)

def primer_match(fprimers, rprimers, R1seq_obj, R2seq_obj, mismatches): # primers is a list (containing one item if a primer is not degenerate)
    '''
    This basically accomplishes inexact string matching, while looking
    for the best match. It will search for an inexact forward primer in 
    the first i positions of the R1 read and for an inexact reverse
    primer in the first j positions of the R2 read, with i and j being
    the lengths of the (degenerage) forward and reverse primers,
    respectively.
    
    Returs a tuple with two biopython seq objects
    (R1_striped, R2_striped) which will have the primers stripted.
    '''
    
    R1_len = len(fprimers[0])
    R2_len = len(rprimers[0])
    R1_fragment = str(R1seq_obj.seq)[:R1_len].upper()
    R2_fragment = str(R2seq_obj.seq)[:R2_len].upper()
    R1_critical = R1_len - mismatches
    R2_critical = R2_len - mismatches
    R1_found = False
    R2_found = False
    for fprimer in fprimers:
        if R1_found == True:
            break
        count = 0
        for i in range(R1_len):
            if fprimer[i] == R1_fragment[i]:
                count += 1
                if count >= R1_critical:
                    R1_found = True
                    break
    for rprimer in rprimers:
        if R2_found == True:
            break
        count = 0
        for i in range(R2_len):
            if rprimer[i] == R2_fragment[i]:
                count += 1
                if count >= R2_critical:
                    R2_found = True
                    break
    if R1_found == True and R2_found == True:
        R1_return = R1seq_obj[R1_len:]
        R2_return = R2seq_obj[R2_len:]
        return(R1_return,R2_return)
    else:
        return(False,False)

usage = """
strip_degen_primer.py -R1 <R1infile> -R2 <R2infile> -p1 <R1primer> -p2 \
<R2primer> -R1o <R1outfile> -R2o <R2outfile> -f <format> \
-max_mismatches <integer>

<infiles>:  FASTQ or FASTA* format paired reads. May optionally be gzip
            compressed if the file name ends with '.gz'. Primers will be
            striped from the reads in these files. If no primer is found
            in either of the pairs, the read pair will be discarded.

<primers>:  Degenerate primer sequence (IUPAC DNA degenerate code). It
            is fine if the primer is not degereate once it then contains
            only (A, T, G, or C). The primer sequences are case
            insensitive.

<format>:   FASTQ is assumed (default). Fasta can be specified instead.
            Only fastq and fasta are supported here. Input and output
            will both be in this format.

<outfiles>: The output files for the R1 and R2 striped reads. Gzip
            compression is not supported for these.

max_mismatches: Up to <int> maximum mismatches are allowed to the primer.
                Default = 3
"""

if '--help' in sys.argv or '-h' in sys.argv:
    print(usage)
    quit()
opts = {}
for arg in ['-R1','-R2','-p1','-p2','-R1o','-R2o']:
    if not arg in sys.argv:
        print "Error: you must specify an option with '%s'" % arg
        print "Use '--help' or '-h' for more details"
        quit()
    opts[arg.lstrip('-')] = sys.argv[sys.argv.index(arg)+1]
if '-f' in sys.argv:
    seq_fmt = sys.argv[sys.argv.index('-f')+1].lower()
    if not seq_fmt in ['fastq','fasta','fas']:
        print "Invalid format specified. It must be either fastq or fasta"
        quit()
    if seq_fmt == 'fas':
        seq_fmt = 'fasta'
else:
    seq_fmt = 'fastq'
if '-max_mismatches' in sys.argv:
    max_mismatches = int(sys.argv[sys.argv.index('-max_mismatches')+1])
else:
    max_mismatches = 3

R1outs = []
R2outs = []
if opts['R1'].endswith('.gz'):
    R1_reads = SeqIO.parse(gzip.open(opts['R1']),seq_fmt)
else:
    R1_reads = SeqIO.parse(opts['R1'],seq_fmt)
if opts['R2'].endswith('.gz'):
    R2_reads = SeqIO.parse(gzip.open(opts['R2']),seq_fmt)
else:
    R2_reads = SeqIO.parse(opts['R2'],seq_fmt)
fprimers = expand_degenerate(opts['p1'])
rprimers = expand_degenerate(opts['p2'])

total = 0
matched = 0
not_matched = 0
for R1rec in R1_reads:
    R2rec = R2_reads.next()
    R1out,R2out = primer_match(fprimers, rprimers, R1rec, R2rec, max_mismatches)
    total += 1
    if R1out and R2out:
        matched += 1
        R1outs.append(R1out)
        R2outs.append(R2out)
    else:
        not_matched += 1
print "Processed %i read pairs" % total
print "Striped primers from %i pairs" % matched
print "Discarded %i pairs not matching to primers" % not_matched
print "Now writing output"
_ = SeqIO.write(R1outs,opts['R1o'],seq_fmt)
_ = SeqIO.write(R2outs,opts['R2o'],seq_fmt)
