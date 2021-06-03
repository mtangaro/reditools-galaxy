#!/usr/bin/env python
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
import sys
import gzip
import argparse
import os
from multiprocessing import Process, Queue
from queue import Empty
import time

list_var = ['AC', 'AG', 'AT', 'CA', 'CG', 'CT', 'GA', 'GC', 'GT', 'TA', 'TC', 'TG']

try:
    import pysam
except:
    sys.exit('Pysam module not found.')
# from multiprocessing import Process, Queue

pysamVersion = pysam.__version__

sys.stderr.write('Pysam version used: %s\n' % pysamVersion)
version = '1.2'

parser = argparse.ArgumentParser(description="REDItoolDnaRna")

parser.add_argument("-i", "--input-rna", help="RNA-Seq BAM file")
parser.add_argument("-j", "--input-dna", nargs='*', help="DNA-Seq BAM file(s separated by space) or folder")
parser.add_argument("-I", "--input-sort-rna", action='store_true', help="Sort input RNA-Seq BAM file")
parser.add_argument("-J", "--input-sort-dna", action='store_true', help="Sort input DNA-Seq BAM file")
parser.add_argument("-f", "--reference", help="Reference in fasta file")
parser.add_argument("-C", "--base-interval", type=int, default=100000, help="Base interval to explore (100000)")
parser.add_argument("-k", "--list-chrm", nargs='*', type=str,
                    help="List of chromosomes to skip (separated by space ex: chr1 chr2)")
parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads")
parser.add_argument("-Y", "--region", type=str,
                    help="Work Only On Region: chrxx:start-end (positions are distributed by the number of threads)")
parser.add_argument("-o", "--output-folder", help="Output folder")
parser.add_argument("-F", "--internal-folder", type=str, default='reditools' ,help="Internal folder name")
parser.add_argument("-M", "--save-list-qual", action="store_true", help="Save a list of columns with quality scores")
parser.add_argument("-c", "--min-read-coverage", nargs=2, default=[10, 10], type=int,
                    help="Min. read coverage (dna rna) [10 10]")
parser.add_argument("-q", "--min-qual-score", nargs=2, default=[30, 30], type=int, help="Min. quality score (dna rna) [30 30]")
parser.add_argument("-m", "--min-map-qual-score", nargs=2, default=[30, 30], type=int,
                    help="Min. mapping quality score (dna rna) [30 30] | Bowtie use 255 | Bowtie2 use 40 | BWA use 30 | RNA-STAR use 255 | HiSAT2 use 60 | Tophat1 use 255 | Tophat2 use 50 | GSNAP use 30")
parser.add_argument("-O", "--min-homo-length", nargs=2, default=[5, 5], type=int,
                    help="Min. homoplymeric length (dna rna) [5 5]")
parser.add_argument("-s", "--infer-strand", type=int, default=1, help="Infer strand (for strand oriented reads) [1]")
parser.add_argument("-g", "--strand-type", type=int, default=1, help="Strand inference (1:maxValue) (2:useConfidence) [1]")
parser.add_argument("-x", "--strand-conf", type=float, default=0.70, help="Strand confidence [0.70]")
parser.add_argument("-S", "--strand-corr", action="store_true", help="Strand correction")
parser.add_argument("-G", "--infer-gff",
                    help="Infer strand by GFF annotation (must be GFF and sorted, otherwise use -X)")
parser.add_argument("-K", "--exclude-gff",
                    help="GFF File with positions to exclude (must be GFF and sorted, otherwise use -X)")
parser.add_argument("-T", "--include-gff",
                    help="Work only on given GFF positions (must be GFF and sorted, otherwise use -X)")
parser.add_argument("-X", "--sort-ann-file", action="store_true", help="Sort annotation files")
parser.add_argument("-e", "--ex-mh-rna", action="store_true", help="Exclude multi hits in RNA-Seq")
parser.add_argument("-E", "--ex-mh-dna", action="store_true", help="Exclude multi hits in DNA-Seq")
parser.add_argument("-d", "--ex-dup-rna", action="store_true", help="Exclude duplicates in RNA-Seq")
parser.add_argument("-D", "--ex-dup-dna", action="store_true", help="Exclude duplicates in DNA-Seq")
parser.add_argument("-p", "--use-pair-rna", action="store_true", help="Use paired concardant reads only in RNA-Seq")
parser.add_argument("-P", "--use-pair-dna", action="store_true", help="Use paired concardant reads only in DNA-Seq")
parser.add_argument("-u", "--cons-map-qual-rna", action="store_true", help="Consider mapping quality in RNA-Seq")
parser.add_argument("-U", "--cons-map-qual-dna", action="store_true", help="Consider mapping quality in DNA-Seq")
parser.add_argument("-a", "--trim-rna", nargs=2, default=[0, 0], type=int,
                    help="Trim x bases up and y bases down per read [0 0] in RNA-Seq")
parser.add_argument("-A", "--trim-dna", nargs=2, default=[0, 0], type=int,
                    help="Trim x bases up and y bases down per read [0 0] in DNA-Seq")
parser.add_argument("-b", "--blat-corr-rna", action="store_true", help="Blat folder for correction in RNA-Seq")
parser.add_argument("-B", "--blat-corr-dna", action="store_true", help="Blat folder for correction in DNA-Seq")
parser.add_argument("-l", "--rm-subs-homo-rna", action="store_true",
                    help="Remove substitutions in homopolymeric regions in RNA-Seq")
parser.add_argument("-L", "--rm-subs-homo-dna", action="store_true",
                    help="Remove substitutions in homopolymeric regions in DNA-Seq")
parser.add_argument("-v", "--min-reads-sup-var", type=int, default=3,
                    help="Min. num. of reads supporting the variation for RNA-Seq [3]")
parser.add_argument("-n", "--min-ed-freq-rna", type=float, default=0.10, help="Min. editing frequency for RNA-Seq [0.10]")
parser.add_argument("-N", "--min-ed-freq-dna", type=float, default=0.10, help="Min. editing frequency for DNA-Seq [0.10]")
parser.add_argument("-z", "--ex-pos-mult-ch-rna", action="store_true",
                    help="Exclude positions with multiple changes in RNA-Seq")
parser.add_argument("-Z", "--ex-pos-mult-ch-dna", action="store_true",
                    help="Exclude positions with multiple changes in DNA-Seq")
parser.add_argument("-W", "--list-var", nargs='*', type=str, default=list_var,
                    help="Select RNA-Seq positions with defined changes (separated by space ex: AG TC) [default all]")
parser.add_argument("-R", "--ex-invar-pos", action="store_true", help="Exclude invariant RNA-Seq positions")
parser.add_argument("-V", "--ex-sites-not-dna", action="store_true", help="Exclude sites not supported by DNA-Seq")
parser.add_argument("-w", "--file-splice", help="File containing splice sites annotations")
parser.add_argument("-r", "--n-bases-near-splice", type=int, default=4,
                    help="Num. of bases near splice sites to explore [4]")
parser.add_argument("-H", "--no-head", action="store_true", help="No Table Header")
parser.add_argument("-GZ", "--gzip", action="store_true", help="Gzip output files")
parser.add_argument("-RR", "--reads", action="store_true", help="Get reads containing nuc. changes")
parser.add_argument("-FQ", "--fastq", nargs='*',
                    help="Fastq to get reads [requires -RR or --reads], separated by space [if paired]")
parser.add_argument("-AP", "--addP", action="store_true", help="Add positions for reads")
parser.add_argument("-RMO", "--rmOver", action="store_true", help="Remove overlapping reads")
parser.add_argument("-RMI", "--rmIndels", action="store_true", help="Remove positions with surrounding Indels")

args = parser.parse_known_args()[0]

# -i --input-rna
bamfile = args.input_rna  # input rna
# -j --input-dna
CTRLbamfile = args.input_dna  # input dna or folder
if CTRLbamfile:
    if os.path.isdir(CTRLbamfile[0]):
        gbamfile = [(os.path.join(CTRLbamfile, x), 0) for x in os.listdir(CTRLbamfile) if x[:-4] == '.bam']
    else:
        gbamfile = [(x, 0) for x in CTRLbamfile]
    dgbamfile = dict(gbamfile)
else:
    gbamfile = []
    dgbamfile = {}
# -I --input-sort-rna
sortbam = args.input_sort_rna
# -J --input-sort-dna
sortgbam = args.input_sort_dna
# -f --reference
fastafile = args.reference
# -c --base-interval
chunckval = args.base_interval
# -k --list-chrm
CTRLnochrs = args.list_chrm
if CTRLnochrs:
    nochrs = CTRLnochrs
else:
    nochrs = []
# -t --threads
NCPU = args.threads
# -Y --region
fworkR = False
CTRLworkR = args.region
if CTRLworkR:
    workR = ('', [0, 0])
    try:
        wline = CTRLworkR.split(':')
        wchr = wline[0]
        wcoord = [int(x.replace(',', '')) for x in wline[1].split('-')]
        workR = (wchr, wcoord)
    except:
        sys.exit('Working region not correct. Use the format chrxx:start-end')
    fworkR = True
# -o --output-folder
outfolder_ = args.output_folder
# -F --internal-folder
infolder = args.internal_folder
# -M --save-list-qual
slist = args.save_list_qual
# -c --min-read-coverage
CTRLclist = args.min_read_coverage
MINCOV = CTRLclist[1]
gMINCOV = CTRLclist[0]
# -q --min-qual-score
CTRLqlist = args.min_qual_score
MQUAL = CTRLqlist[1]
gMQUAL = CTRLqlist[0]
# -m --min-map-qual-score
CTRLmlist = args.min_map_qual_score
MAPQ = CTRLmlist[1]
gMAPQ = CTRLmlist[0]
# -O --min-homo-length
CTRLhomo = args.min_homo_length
homo = CTRLhomo[1]
ghomo = CTRLhomo[0]
# -s --infer-strand
CTRLstrand = args.infer_strand
if CTRLstrand == 1:
    getstrand, unchange1, unchange2 = 1, 1, 0
elif CTRLstrand == 0:
    getstrand, unchange1, unchange2 = 1, 0, 0
elif CTRLstrand == 2:
    getstrand, unchange1, unchange2 = 1, 0, 1
elif CTRLstrand == 12:
    getstrand, unchange1, unchange2 = 1, 1, 1
# -g --strand-type
CTRLuseconf = args.strand_type
if CTRLstrand == 1:
    useconf = 0
elif CTRLstrand == 2:
    useconf = 1
# -x --strand-conf
strconf = args.strand_conf
# -S --strand-corr
corrstr = args.strand_corr
# -G --infer-gff
uann = False
CTRLannfile = args.infer_gff
if CTRLannfile:
    annfile = CTRLannfile
    uann = True
# -K --exclude-gff
expos = False
CTRLexfile = args.exclude_gff
if CTRLexfile:
    exfile = CTRLexfile
    if os.path.exists(exfile):
        expos = True
# -T --include-gff
uwf = False
CTRLwfile = args.include_gff
if CTRLwfile:
    wfile = CTRLwfile
    if os.path.exists(wfile):
        uwf = True
# -X --sort-ann-file
sortann = args.sort_ann_file
# -e --ex-mh-rna
exh = args.ex_mh_rna
# -E --ex-mh-dna
gexh = args.ex_mh_dna
# -d --ex-dup-rna
exd = args.ex_dup_rna
# -D --ex-dup-dna
gexd = args.ex_dup_dna
# -p --use-pair-rna
conc = args.use_pair_rna
# -P --use-pair-dna
gconc = args.use_pair_dna
# -u --cons-map-qual-rna
mq = args.cons_map_qual_rna
# -U --cons-map-qual-dna
gmq = args.cons_map_qual_dna
# -a --trim-rna
rmnuc = False
rmp = args.trim_rna
if rmp != [0, 0]:
    rmnuc = True
# -A --trim-dna
grmnuc = False
grmp = args.trim_dna
if grmp != [0, 0]:
    grmnuc = True
# -b --blat-corr-rna
blatr = False
CTRLblatfolder = args.blat_corr_rna
if CTRLblatfolder:
    blatfolder = CTRLblatfolder
    if os.path.exists(blatfolder):
        blatr = True
# -B --blat-corr-dna
gblatr = False
CTRLgblatfolder = args.blat_corr_dna
if CTRLgblatfolder:
    gblatfolder = CTRLgblatfolder
    if os.path.exists(gblatfolder):
        gblatr = True
# -l --rm-subs-homo-rna
rmsh = args.rm_subs_homo_rna
# -L --rm-subs-homo-dna
grmsh = args.rm_subs_homo_dna
# -v --min-reads-sup-var
vnuc = args.min_reads_sup_var
# -n --min-ed-freq-rna
mmf = args.min_ed_freq_rna
# -N --min-ed-freq-dna
gmmf = args.min_ed_freq_dna
# -z --ex-pos-mult-ch-rna
exms = args.ex_pos_mult_ch_rna
# -Z --ex-pos-mult-ch-dna
exnonh = args.ex_pos_mult_ch_dna
# -W --list-var
usubs = args.list_var
# -R --ex-invar-pos
exinv = args.ex_invar_pos
# -V --ex-sites-not-dna
exnosupp = args.ex_sites_not_dna
# -w --file-splice
exss = False
CTRLsplicefile = args.file_splice
if CTRLsplicefile:
    splicefile = CTRLsplicefile
    if os.path.exists(splicefile):
        exss = True
# -r --n-bases-near-splice
nss = args.n_bases_near_splice
# -H --no-head
noheader = args.no_head
# --gzip -GZ
gziptab = args.gzip
# --reads -RR
greads = args.reads
# --fastq -FQ
fastq = args.fastq
# --addP -AP
addP = args.addP
# --rmOver -RMO
rmOver = args.rmOver
# --rmIndels -RMI
rmIndel = args.rmIndels

MAX_DEPTH = 100000

script_time = time.strftime("%d/%m/%Y %H:%M:%S", time.localtime(time.time()))

param = args.__dict__
params = []
for x in param:
    s='%s : %s\n' %(x, param[x])
    params.append(s)

def BaseCount(seq, ref, mfr, VNUC):
    b = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    subs = []
    subv = []
    for i in seq.upper():
        if b.has_key(i): b[i] += 1
    for i in b:
        if not b.has_key(ref): continue
        if b[i] != 0 and i != ref:
            vv = float(b[i]) / (b[i] + b[ref])
            subv.append((b[i], vv, ref + i))
    subv.sort()
    subv.reverse()
    for i in subv:
        if i[0] >= VNUC and i[1] >= mfr: subs.append(i[2])
    freq = 0.0
    if len(subs) == 0:
        subs.append('-')
    else:
        freq = subv[0][1]
    return sum(b.values()), [b['A'], b['C'], b['G'], b['T']], ' '.join(subs), '%.2f' % (freq)

def meanq(v, n):
    try:
        m = float(v) / n
    except:
        m = 0.0
    return '%.2f' % (m)

def rmHomo(sequp, seqdw, gh, ref):
    if len(sequp) == 0 and len(seqdw) == 0: return 0
    up, dw = 0, 0
    for i in seqdw:
        if i == ref:
            dw += 1
        else:
            break
    for i in sequp[::-1]:
        if i == ref:
            up += 1
        else:
            break
    hlen = up + dw + 1
    if hlen >= gh: return 1
    return 0


def prop(tot, va):
    try:
        av = float(va) / tot
    except:
        av = 0.0
    return av


def vstand(strand):
    vv = [(strand.count('+'), '+'), (strand.count('-'), '-'), (strand.count('*'), '*')]
    if vv[0][0] == 0 and vv[1][0] == 0: return '*'
    if useconf:
        totvv = sum([x[0] for x in vv[:2]])
        if prop(totvv, vv[0][0]) >= strconf: return '+'
        if prop(totvv, vv[1][0]) >= strconf: return '-'
        return '*'
    else:
        if vv[0][0] == vv[1][0] and vv[2][0] == 0: return '+'
        return max(vv)[1]


def comp(s):
    a = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    ss = ''
    for i in s.upper():
        if a.has_key(i):
            ss += a[i]
        elif i == ' ':
            ss += ' '
        elif i == '-':
            ss += '-'
        else:
            ss += 'N'
    return ss


def comp2(s):
    ss = {}
    a = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    for i, j in enumerate('ACGT'): ss[a[j]] = s[i]
    return str([ss['A'], ss['C'], ss['G'], ss['T']])


def whereis(program):
    for path in os.environ.get('PATH', '').split(':'):
        if os.path.exists(os.path.join(path, program)) and not os.path.isdir(os.path.join(path, program)): return 1
    return 0

def vstrand(lista):
    if len(lista) == 0: return '2'
    p = lista.count('+')
    m = lista.count('-')
    if p == len(lista):
        return '1'
    elif m == len(lista):
        return '0'
    else:
        return '2'

def checkSubs(s):
    if s == '-': return 1
    for i in s.split():
        if i in usubs: return 1
    return 0

def makeCluster(allcoord):
    cluster = []
    remaining = []
    c1 = allcoord[0][0]
    c2 = allcoord[0][1]
    for i in range(len(allcoord)):
        if allcoord[i] != (c1, c2):
            if c1 <= allcoord[i][0] <= c2:
                cluster.append(allcoord[i])
                if allcoord[i][1] > c2:
                    c2 = allcoord[i][1]
            else:
                remaining.append(allcoord[i])
        else:
            cluster.append((c1, c2))
    return (c1, c2), remaining


def newCoords(interval, start, end):
    coords = []
    interval.sort()
    while len(interval) != 0:
        coord, interval = makeCluster(interval)
        coords.append(coord)
    c1, c2 = coords[0][0], coords[-1][1]
    if c1 < start: c1 = start
    if c2 > end: c2 = end
    if c1 == c2: c1 = start - 1  # MODIFICATO
    return coords, c1, c2


def checkPos(coords, pos):
    for i in coords:
        if i[0] <= pos <= i[1]: return 1
    return 0


def parseFeat(line):
    l = line.split('\t')
    cc = (int(l[3]) - 1, int(l[4]) - 1)
    return cc


def normByStrand(seq_, strand_, squal_, mystrand_):
    st = '+'
    if mystrand_ == '0': st = '-'
    seq, qual, squal = '', 0, []
    for i in range(len(seq_)):
        if strand_[i] == st:
            seq += seq_[i]
            qual += squal_[i] - 33
            squal.append(squal_[i])
    return seq, qual, squal


def normByBlat(seq_, strand_, squal_, blatc_, qqval):
    seq, qual, squal, strand = '', 0, [], ''
    for i in range(len(seq_)):
        if blatc_[i] == '1':
            seq += seq_[i]
            qual += squal_[i]
            squal.append(squal_[i])
            strand += strand_[i]
    return seq, qual, squal, strand


def normByOverlap(seq_, strand_, squal_, blatc_, qqval, over_):
    seq, qual, squal, strand, blatc = '', 0, [], '', []
    for i in range(len(seq_)):
        if over_[i] == 0:
            seq += seq_[i]
            qual += squal_[i]
            squal.append(squal_[i])
            strand += strand_[i]
            blatc.append(blatc_[i])
    return seq, qual, squal, strand, blatc


def testBlat(blc):
    if blc.count('1') > blc.count('0'): return 1
    return 0


def countIndels(lista):
    for i in lista:
        if i.count(None) > 0: return 1
    return 0


def getOverlap(lista):
    r = [0 for x in range(len(lista))]
    l = [x[0] for x in lista]
    us = {}
    x = 0
    for i in lista:
        if l.count(i[0]) == 2:
            s = '='
            if i[1] != i[2]: s = '!'
            if us.has_key(i[0]):
                us[i[0]].append((x, s))
            else:
                us[i[0]] = [(x, s)]
        x += 1
    for i in us:
        v = us[i]
        if v[0][1] == v[1][1]:
            r[v[0][0]] = 1
        else:
            if v[0][1] == '!':
                r[v[0][0]] = 1
            elif v[1][1] == '!':
                r[v[1][0]] = 1
    return r

###########################################################
###########################################################
script_time = time.strftime("%d/%m/%Y %H:%M:%S", time.localtime(time.time()))
sys.stderr.write("Script time --> START: %s\n" % (script_time))
sys.stderr.write("Analysis name: %s\n" % (outfolder_))
###########################################################


if not os.path.exists(bamfile):
    sys.exit('RNA-Seq BAM file %s not found.' % (bamfile))
if sortbam:
    sys.stderr.write('Sorting RNA-Seq BAM file.\n')
    pysam.sort(bamfile, '%s_sorted' % (outfolder_))
    os.rename(bamfile, bamfile + '_old')
    os.rename('%s_sorted.bam' % (outfolder_), bamfile)
    sys.stderr.write('Indexing RNA-Seq BAM file.\n')
    pysam.index(bamfile)
if not os.path.exists(bamfile + '.bai') and not sortbam:
    sys.stderr.write('Indexing RNA-Seq BAM file.\n')
    pysam.index(bamfile)
###########################################################

dgdic = {}  # dizionario chr:bam file
for i in gbamfile:
    if not os.path.exists(i[0]):
        sys.stderr.write('DNA-Seq BAM file %s not found.\n' % (i[0]))
        sys.stderr.write('Working without DNA-Seq BAM file %s.\n' % (i[0]))
        del dgbamfile[i[0]]
    else:
        if sortgbam:
            sys.stderr.write('Sorting DNA-Seq BAM file %s.\n' % (i[0]))
            pysam.sort(i[0], '%s_sorted' % (outfolder_))
            os.rename(i[0], i[0] + '_old')
            os.rename('%s_sorted.bam' % (outfolder_), i[0])
            sys.stderr.write('Indexing DNA-Seq BAM file %s.\n' % (i[0]))
            pysam.index(i[0])
        if not os.path.exists(i[0] + '.bai') and not sortgbam:
            sys.stderr.write('Indexing DNA-Seq BAM file %s.\n' % (i[0]))
            pysam.index(i[0])
if len(gbamfile) == 0:
    sys.stderr.write('Working without DNA-Seq BAM file(s).\n')
    nogbam = 1
else:
    for i in dgbamfile:
        idxinfo = pysam.idxstats(i)
        for j in idxinfo.split('\n'):  # MOD
            l = (j.strip()).split('\t')
            if l[0] in ['*', '']: continue  # MOD
            if int(l[2]) + int(l[3]) > 0: dgdic[l[0]] = i

###########################################################
if not os.path.exists(fastafile):
    sys.exit('Fasta file %s not found.' % (fastafile))
if not os.path.exists(fastafile + '.fai'):
    sys.stderr.write('Indexing Fasta file.\n')
    pysam.faidx(fastafile)
###########################################################
# Check reference for name consistency
grefs = dgdic.keys()
rrefs = {}
ridxinfo = pysam.idxstats(bamfile)
for j in ridxinfo.split('\n'):  # MOD
    l = (j.strip()).split('\t')
    if l[0] in ['*', '']: continue  # MOD
    if int(l[2]) + int(l[3]) > 0: rrefs[l[0]] = int(l[1])
frefs = []
fidxinfo = open(fastafile + '.fai')
for j in fidxinfo:
    l = (j.strip()).split('\t')
    if l[0] == '': continue
    frefs.append(l[0])
fidxinfo.close()
# in rna-seq
rnof = []
for i in rrefs.keys():
    if i not in frefs: sys.stderr.write('WARNING: Region %s in RNA-Seq not found in reference file.\n' % (i))
if len(gbamfile) != 0:
    for i in grefs:
        if i not in frefs: sys.stderr.write('WARNING: Region %s in DNA-Seq not found in reference file.\n' % (i))
###########################################################

###########################################################
# Annotation file for working regions
if uwf:
    if not os.path.exists(wfile):
        sys.exit('GFF file %s not found.' % (wfile))
    if sortann:
        if not whereis('grep'): sys.exit('grep command not found.')
        if not whereis('sort'): sys.exit('sort command not found.')
        sys.stderr.write('Sorting GFF file.\n')
        scmd = 'grep ^"#" %s; grep -v ^"#" %s | sort -k1,1 -k4,4n > %s' % (wfile, wfile, '%s_workf' % (outfolder_))
        os.system(scmd)
        os.rename(wfile, wfile + '_old')
        os.rename('%s_workf' % (outfolder_), wfile)
    if not os.path.exists(wfile + '.tbi'):
        sys.stderr.write('Indexing GFF file.\n')
        wfile = pysam.tabix_index(wfile, preset='gff')
###########################################################
# Annotation file for strand detection
if uann:
    getstrand = 0
    if not os.path.exists(annfile):
        sys.exit('Annotation file %s not found.' % (annfile))
    if sortann:
        if not whereis('grep'): sys.exit('grep command not found.')
        if not whereis('sort'): sys.exit('sort command not found.')
        sys.stderr.write('Sorting annotation file.\n')
        scmd = 'grep ^"#" %s; grep -v ^"#" %s | sort -k1,1 -k4,4n > %s' % (annfile, annfile, '%s_annotation' % (outfolder_))
        os.system(scmd)
        os.rename(annfile, annfile + '_old')
        os.rename('%s_annotation' % (outfolder_), annfile)
    if not os.path.exists(annfile + '.tbi'):
        sys.stderr.write('Indexing annotation file.\n')
        annfile = pysam.tabix_index(annfile, preset='gff')
###########################################################
# Annotation file to exclude positions
if expos:
    if not os.path.exists(exfile):
        sys.exit('File %s not found.' % (exfile))
    if sortann:
        if not whereis('grep'): sys.exit('grep command not found.')
        if not whereis('sort'): sys.exit('sort command not found.')
        sys.stderr.write('Sorting file.\n')
        scmd = 'grep ^"#" %s; grep -v ^"#" %s | sort -k1,1 -k4,4n > %s' % (exfile, exfile, '%s_exfile' % (outfolder_))
        os.system(scmd)
        os.rename(exfile, exfile + '_old')
        os.rename('%s_exfile' % (outfolder_), exfile)
    if not os.path.exists(exfile + '.tbi'):
        sys.stderr.write('Indexing %s file.\n' % (exfile))
        exfile = pysam.tabix_index(exfile, preset='gff')
###########################################################
dicregions = dict(rrefs.items())
chrs = [x for x in dicregions.keys() if x not in nochrs]
if fworkR:
    sys.stderr.write('Analysis on region %s:%i-%i.\n' % (workR[0], workR[1][0], workR[1][1]))
else:
    sys.stderr.write('Analysis on %i regions.\n' % (len(chrs)))
###########################################################
if infolder != '':
    outfolder = os.path.join(outfolder_, 'DnaRna_%s_%s' % (infolder, outfolder_))
else:
    outfolder = os.path.join(outfolder_, 'DnaRna_%s' % (outfolder_))
if not os.path.exists(outfolder):
    splitfolder = os.path.split(outfolder)
    if not os.path.exists(splitfolder[0]): os.mkdir(splitfolder[0])
    os.mkdir(outfolder)
outtable = os.path.join(outfolder, '%s_outTable' % (outfolder_))
if slist:
    slistfile = os.path.join(outfolder, '%s_outPileupRNA' % (outfolder_))
    if len(gbamfile) != 0: gslistfile = os.path.join(outfolder, '%s_outPileupDNA' % (outfolder_))
# write command line and input parameters
f = open(os.path.join(outfolder, 'parameters.txt'), 'w')
f.writelines(params)
f.close()

###########################################################
d = {}
if blatr:
    badblat = blatfolder  # os.path.join(blatfolder,'blatseqs_%s.bad'%(chr))
    if os.path.exists(badblat):
        sys.stderr.write('Using Blat mapping...\n')
        f = open(badblat)
        for i in f:
            l = (i.strip()).split()
            d[l[0] + '_' + l[1]] = int(l[1])
        f.close()
        sys.stderr.write('Found %i reads.\n' % (len(d)))

def exploreBAM(myinput):
    isgbam = 1
    inputs = myinput.split('$')
    chr, bamfile, start_region, lenregion, suff_ = inputs[0], inputs[1], int(inputs[2]), int(inputs[3]), inputs[4]
    if not dgdic.has_key(chr): isgbam = 0
    outfile = os.path.join(outfolder, 'table_%s' % (suff_))
    if slist:
        if gziptab:
            outrna = gzip.open(os.path.join(outfolder, 'pileupRNA_%s.gz' % (suff)), 'wb')
        else:
            outrna = open(os.path.join(outfolder, 'pileupRNA_%s' % (suff)), 'w')
        if not nogbam and isgbam:
            if gziptab:
                outdna = gzip.open(os.path.join(outfolder, 'pileupDNA_%s.gz' % (suff)), 'wb')
            else:
                outdna = open(os.path.join(outfolder, 'pileupDNA_%s' % (suff)), 'w')
    # d,di,gd={},{},{}
    di, gd = {}, {}
    bam = pysam.Samfile(bamfile, "rb")
    if not nogbam and isgbam:
        gbam = pysam.Samfile(dgdic[chr], "rb")
    fasta = pysam.Fastafile(fastafile)
    if uann: tabix = pysam.Tabixfile(annfile)
    if expos: extabix = pysam.Tabixfile(exfile)
    if uwf: wtabix = pysam.Tabixfile(wfile)
    if gziptab:
        outfile = outfile + '.gz'
        out = gzip.open(outfile, 'wb')
    else:
        out = open(outfile, 'w')
    sys.stderr.write('Started analysis on region: %s:%i-%i\n' % (chr, start_region + 1, lenregion))
    # sys.stderr.write('OUTFILE: %s\n' %(outfile))
    #	if blatr:
    #		badblat=os.path.join(blatfolder,'blatseqs_%s.bad'%(chr))
    #		if os.path.exists(badblat):
    #			sys.stderr.write('Using Blat mapping for region %s\n'%(chr))
    #			f=open(badblat)
    #			for i in f:
    #				l=(i.strip()).split()
    #				d[l[0]+'_'+l[1]]=int(l[1])
    #			f.close()
    #			sys.stderr.write('Found %i reads for region %s\n'%(len(d),chr))
    if gblatr:
        gbadblat = os.path.join(gblatfolder, 'blatseqs_%s.bad' % (chr))
        if os.path.exists(gbadblat):
            sys.stderr.write('Using Blat mapping for DNA region %s\n' % (chr))
            f = open(gbadblat)
            for i in f:
                l = (i.strip()).split()
                gd[l[0] + '_' + l[1]] = int(l[1])
            f.close()
            sys.stderr.write('Found %i reads for region %s\n' % (len(gd), chr))
    if exss:
        if os.path.exists(splicefile):
            sys.stderr.write('Loading known splice sites for region %s\n' % (chr))
            f = open(splicefile)
            for i in f:
                l = (i.strip()).split()
                if l[0] != chr: continue
                st, tp, cc = l[4], l[3], int(l[1])
                if st == '+' and tp == 'D':
                    for j in range(nss): di[cc + (j + 1)] = 0
                if st == '+' and tp == 'A':
                    for j in range(nss): di[cc - (j + 1)] = 0
                if st == '-' and tp == 'D':
                    for j in range(nss): di[cc - (j + 1)] = 0
                if st == '-' and tp == 'A':
                    for j in range(nss): di[cc + (j + 1)] = 0
            f.close()
            sys.stderr.write('Loaded %i positions for %s\n' % (len(di), chr))
    if greads:
        outreads = open(os.path.join(outfolder, 'readsRNA_%s' % (suff_)), 'w')
        grdb = {}
    if addP:
        outAddP = open(os.path.join(outfolder, 'readsPosRNA_%s' % (suff_)), 'w')
        grdb2 = {}
    for kpos in range(start_region, lenregion, chunckval):
        startk, endk = kpos, (kpos + chunckval) - 1
        if endk > lenregion: endk = lenregion - 1
        # sys.stderr.write('%i %i\n'%(startk,endk))
        # check features in the give region if GFF provided
        if uwf:
            if chr in wtabix.contigs:
                wfeat = [parseFeat(feat) for feat in wtabix.fetch(reference=chr, start=startk, end=endk)]
                # print wfeat
                if len(wfeat) == 0: continue
                wcoords, startk, endk = newCoords(wfeat, startk, endk)
            else:
                continue
        # get FASTA sequence
        # print startk,endk
        # refgenome=fasta.fetch(chr,startk,endk+1).upper()
        # print refgenome
        # explore dna-seq bam
        #####################
        gdic = {}
        if not nogbam and isgbam:
            for pileupcolumn in gbam.pileup(chr, startk, endk, stepper='nofilter', max_depth=MAX_DEPTH):
                if uwf and not checkPos(wcoords, pileupcolumn.reference_pos): continue
                if not startk <= pileupcolumn.reference_pos <= endk: continue
                gref = fasta.fetch(chr, pileupcolumn.reference_pos, pileupcolumn.reference_pos + 1).upper()
                gseq, gqual, gstrand, gblatc, gsqual = '', 0, '', '', []
                if grmsh:
                    if ((pileupcolumn.reference_pos + 1) - ghomo) - 1 < 0:
                        sequp = ''
                    else:
                        sequp = (fasta.fetch(chr, ((pileupcolumn.reference_pos + 1) - ghomo) - 1,
                                             (pileupcolumn.reference_pos + 1) - 1)).upper()
                    seqdw = (fasta.fetch(chr, pileupcolumn.reference_pos + 1,
                                         (pileupcolumn.reference_pos + 1) + ghomo)).upper()
                for pileupread in pileupcolumn.pileups:  # per ogni base dell'allineamento multiplo
                    if pileupread.is_del: continue
                    if pileupread.alignment.is_qcfail: continue
                    # gs,gq,gt,gqq=pileupread.alignment.seq[pileupread.qpos].upper(),ord(pileupread.alignment.qual[pileupread.qpos])-gQVAL,'*',pileupread.alignment.qual[pileupread.qpos]
                    # multiple hits
                    if gexh:
                        if pileupread.alignment.is_secondary: continue
                        if pileupread.alignment.has_tag('NH'):
                            if pileupread.alignment.get_tag('NH') > 1: continue
                    # duplicates
                    if gexd and pileupread.alignment.is_duplicate: continue
                    # se paired end
                    if gconc:  # se devi usare solo le paired
                        # se non sono paired
                        if not pileupread.alignment.is_paired: continue
                        # se non sono concordanti
                        if not pileupread.alignment.is_proper_pair: continue
                        # se concordanti ma nello stesso orientamento
                        flag = pileupread.alignment.flag
                        if pileupread.alignment.is_duplicate: flag = flag - 1024
                        if pileupread.alignment.is_secondary: flag = flag - 256
                        if flag in [67, 131, 115, 179]: continue
                    # mapping quality
                    if gmq and pileupread.alignment.mapping_quality < gMAPQ: continue
                    # se la qualita' >= della qualita' minima
                    gs, gq, gt, gqq = pileupread.alignment.query_sequence[pileupread.query_position].upper(), \
                                      pileupread.alignment.query_qualities[pileupread.query_position], '*', \
                                      pileupread.alignment.query_qualities[pileupread.query_position]
                    if gq >= gMQUAL and pileupcolumn.reference_pos in pileupread.alignment.get_reference_positions():
                        if grmnuc:
                            grlen = pileupread.alignment.query_length  # pileupread.alignment.qlen #lunghezza della specifica read
                            gqp = pileupread.query_position
                            if grmp[0] > 0:  # rimuovi posizioni al 5'
                                if pileupread.alignment.is_reverse:
                                    if (pileupread.alignment.query_alignment_end - grmp[
                                        1]) + 1 <= gqp <= pileupread.alignment.query_alignment_end: continue
                                else:
                                    if pileupread.alignment.query_alignment_start <= gqp <= (
                                            pileupread.alignment.query_alignment_start + grmp[0]) - 1: continue
                            if grmp[1] > 0:  # rimuovi posizioni al 3'
                                if pileupread.alignment.is_reverse:
                                    if pileupread.alignment.query_alignment_start <= gqp <= (
                                            pileupread.alignment.query_alignment_start + grmp[0]) - 1: continue
                                else:
                                    if (pileupread.alignment.query_alignment_end - grmp[
                                        1]) + 1 <= gqp <= pileupread.alignment.query_alignment_end: continue
                        # se la read di appartenenza non mappa in modo univoco con Blat
                        if gblatr:
                            rt = 0
                            if pileupread.alignment.is_read1:
                                rt = 1
                            elif pileupread.alignment.is_read2:
                                rt = 2
                            rname = pileupread.alignment.query_name + '_%i' % (rt)
                            if gd.has_key(rname):
                                gblatc += '0'  # continue
                            else:
                                gblatc += '1'
                        # se la base e' diversa dal reference
                        # se in regione omopolimerica scarta
                        if grmsh and rmHomo(sequp, seqdw, ghomo, gref): continue
                        gseq += gs
                        gqual += gq
                        gstrand += gt
                        gsqual.append(gqq)
                if gseq.strip() != '':
                    if gblatr:
                        if testBlat(gblatc):
                            gseq, gqual, gsqual, gstrand = normByBlat(gseq, gstrand, gsqual, gblatc, 33)
                        else:
                            continue
                    gcov, gbcomp, gsubs, gfreq = BaseCount(gseq, gref, gmmf, 0)
                    if gcov < gMINCOV: continue
                    gmqua = meanq(gqual, len(gseq))
                    ghinfo = 0  # non omozigote
                    if gsubs == '-': ghinfo = 1  # omozigote
                    gdic[pileupcolumn.reference_pos] = ([str(gcov), gmqua, str(gbcomp), gsubs, gfreq], ghinfo)
                    if slist:
                        if not nogbam and isgbam: outdna.write(
                            '\t'.join([chr, str(pileupcolumn.reference_pos + 1), gref, gseq, gsqual]) + '\n')
        #####################
        # explore rna-seq bam
        # print startk,endk
        for pileupcolumn in bam.pileup(chr, startk, endk, stepper='nofilter', max_depth=MAX_DEPTH):
            # print chr,startk,endk
            # print dir(pileupcolumn)
            if uwf and not checkPos(wcoords, pileupcolumn.reference_pos): continue
            if not startk <= pileupcolumn.reference_pos <= endk: continue
            # print
            # print chr,pileupcolumn.reference_pos+1
            ref = fasta.fetch(chr, pileupcolumn.reference_pos, pileupcolumn.reference_pos + 1).upper()
            # seq,qual,strand,squal,blatc='',0,'','',''
            seq, qual, strand, blatc, squal = '', 0, '', '', []
            if rmOver: rall = []
            if rmIndel: indels = []
            if rmsh:
                if ((pileupcolumn.reference_pos + 1) - homo) - 1 < 0:
                    sequp = ''
                else:
                    sequp = (fasta.fetch(chr, ((pileupcolumn.reference_pos + 1) - homo) - 1,
                                         (pileupcolumn.reference_pos + 1) - 1)).upper()
                seqdw = (
                    fasta.fetch(chr, pileupcolumn.reference_pos + 1, (pileupcolumn.reference_pos + 1) + homo)).upper()
            for pileupread in pileupcolumn.pileups:  # per ogni base dell'allineamento multiplo
                # if pileupread.alignment.is_supplementary and not pileupread.is_del: print pileupread.alignment
                if pileupread.is_del: continue
                if pileupread.alignment.is_qcfail: continue
                if pileupread.alignment.is_supplementary: continue
                if pileupread.alignment.has_tag('SA'): continue
                # print pileupread
                # print dir(pileupread)
                # print
                # print dir(pileupread.alignment)
                # print pileupread.alignment.get_tag('NM')
                # s,q,t,qq=pileupread.alignment.query_sequence[pileupread.query_position].upper(),pileupread.alignment.query_qualities[pileupread.query_position],'*',pileupread.alignment.qual[pileupread.query_position]
                # s,q,t,qq=pileupread.alignment.seq[pileupread.qpos].upper(),ord(pileupread.alignment.qual[pileupread.qpos])-QVAL,'*',pileupread.alignment.qual[pileupread.qpos]
                # escludi posizioni introniche nei pressi di splice sites
                if exss and di.has_key(pileupcolumn.reference_pos + 1): continue
                # multiple hit
                if exh:
                    if pileupread.alignment.is_secondary: continue
                    if pileupread.alignment.has_tag('NH'):
                        if pileupread.alignment.get_tag('NH') > 1: continue
                # duplicates
                if exd and pileupread.alignment.is_duplicate: continue
                # se paired end
                if conc:  # se devi usare solo le paired
                    # se non sono paired
                    if not pileupread.alignment.is_paired: continue
                    # se non sono concordanti
                    if not pileupread.alignment.is_proper_pair: continue
                    # se concordanti ma nello stesso orientamento
                    flag = pileupread.alignment.flag
                    if pileupread.alignment.is_duplicate: flag = flag - 1024
                    if pileupread.alignment.is_secondary: flag = flag - 256
                    if flag in [67, 131, 115, 179]: continue
                # print pileupread.alignment.qual
                # print pileupread.alignment.flag
                # print pileupread.alignment.is_paired
                # mapping quality
                if mq and pileupread.alignment.mapping_quality < MAPQ: continue
                # print pileupread.alignment.query_sequence
                if not pileupread.alignment.query_qualities: pileupread.alignment.query_qualities = [30 for vn in range(
                    len(pileupread.alignment.query_sequence))]
                # print pileupread.query_position
                # s,q,t,qq=pileupread.alignment.query_sequence[pileupread.query_position].upper(),pileupread.alignment.query_qualities[pileupread.query_position],'*',pileupread.alignment.qual[pileupread.query_position]
                s, q, t, qq = pileupread.alignment.query_sequence[pileupread.query_position].upper(), \
                              pileupread.alignment.query_qualities[pileupread.query_position], '*', \
                              pileupread.alignment.query_qualities[pileupread.query_position]
                if rmIndel:
                    indelALN = pileupread.alignment.get_aligned_pairs(matches_only=False, with_seq=False)
                    indelreg = indelALN[pileupread.query_position - 5:pileupread.query_position] + indelALN[
                                                                                                   pileupread.query_position + 1:pileupread.query_position + 6]
                    indel = countIndels(indelreg)
                # se la qualita' >= alla qualita' minima
                if q >= MQUAL and pileupcolumn.reference_pos in pileupread.alignment.get_reference_positions():
                    # tags=dict(pileupread.alignment.tags)
                    # deduci la strand per ogni posizione
                    if getstrand:
                        # usa le info del mapping se strand oriented
                        if pileupread.alignment.is_read1:
                            if unchange1:
                                if pileupread.alignment.is_reverse:
                                    t = '-'
                                else:
                                    t = '+'
                            else:
                                if pileupread.alignment.is_reverse:
                                    t = '+'
                                else:
                                    t = '-'
                        elif pileupread.alignment.is_read2:
                            if unchange2:
                                if pileupread.alignment.is_reverse:
                                    t = '-'
                                else:
                                    t = '+'
                            else:
                                if pileupread.alignment.is_reverse:
                                    t = '+'
                                else:
                                    t = '-'
                        else:  # for single ends
                            if unchange1:
                                if pileupread.alignment.is_reverse:
                                    t = '-'
                                else:
                                    t = '+'
                            else:
                                if pileupread.alignment.is_reverse:
                                    t = '+'
                                else:
                                    t = '-'
                    if rmnuc:
                        rlen = pileupread.alignment.query_length  # pileupread.alignment.qlen #lunghezza della specifica read
                        qp = pileupread.query_position
                        # print pileupcolumn.reference_pos+1, qp
                        # print pileupread.alignment
                        # print pileupread.alignment.query_name,pileupread.alignment.get_aligned_pairs(matches_only=False, with_seq=False)
                        # print pileupread.alignment, qp , pileupread.alignment.is_reverse,pileupread.alignment.query_sequence[pileupread.query_position].upper()
                        # print (pileupread.alignment.query_alignment_end-rmp[1]),pileupread.alignment.query_alignment_end-1
                        # print pileupread.alignment.query_alignment_start, (pileupread.alignment.query_alignment_start+rmp[0])-1
                        if rmp[0] > 0:  # rimuovi posizioni al 5'
                            if pileupread.alignment.is_reverse:
                                if (pileupread.alignment.query_alignment_end - rmp[
                                    1]) <= qp <= pileupread.alignment.query_alignment_end - 1: continue
                            else:
                                if pileupread.alignment.query_alignment_start <= qp <= (
                                        pileupread.alignment.query_alignment_start + rmp[0]) - 1: continue
                        if rmp[1] > 0:  # rimuovi posizioni al 3'
                            if pileupread.alignment.is_reverse:
                                if pileupread.alignment.query_alignment_start <= qp <= (
                                        pileupread.alignment.query_alignment_start + rmp[0]) - 1: continue
                            else:
                                if (pileupread.alignment.query_alignment_end - rmp[
                                    1]) <= qp <= pileupread.alignment.query_alignment_end - 1: continue
                    # print qp, rmp
                    # if pileupread.alignment.is_reverse:
                    #	if qp>(rlen-rmp[0])-1: continue
                    #	if qp<rmp[1]:continue
                    # else:
                    #	if qp<rmp[0]:continue
                    #	if qp>(rlen-rmp[1])-1: continue
                    # se la read di appartenenza non mappa in modo univoco con Blat
                    if blatr:
                        rt = 0
                        if pileupread.alignment.is_read1:
                            rt = 1
                        elif pileupread.alignment.is_read2:
                            rt = 2
                        rname = pileupread.alignment.query_name + '_%i' % (rt)
                        if d.has_key(rname):
                            blatc += '0'  # continue
                        else:
                            blatc += '1'
                    # se la base e' diversa dal reference
                    # se in regione omopolimerica scarta
                    if rmsh and rmHomo(sequp, seqdw, homo, ref): continue
                    seq += s
                    qual += q
                    strand += t
                    squal.append(qq)
                    if rmIndel: indels.append(indel)
                    if rmOver: rall.append((pileupread.alignment.query_name, ref, s))
                    if greads:  # --reads option
                        if ref != s:
                            rt = 0
                            if pileupread.alignment.is_read1:
                                rt = 1
                            elif pileupread.alignment.is_read2:
                                rt = 2
                            rqname = pileupread.alignment.query_name + '_%i' % (rt)
                            rname = pileupread.alignment.query_name
                            rseqname = pileupread.alignment.query_sequence
                            # print rqname, rseqname, chr, pileupread.alignment.reference_start,pileupread.alignment.reference_end
                            # rqname_comp=rqname+'$'+chr+'$'+str(pileupread.alignment.reference_start)+'$'+str(pileupread.alignment.reference_end)
                            # print pileupread.alignment.query_name,chr,pileupread.alignment.reference_start,pileupread.alignment.reference_end
                            # print chr,pileupread.alignment.reference_name
                            if not pileupread.alignment.is_unmapped:
                                mate = bam.mate(pileupread.alignment)
                                addpos = (
                                pileupread.alignment.query_name, mate.query_name, pileupread.alignment.reference_name,
                                mate.reference_name, pileupread.alignment.reference_start,
                                pileupread.alignment.reference_end, mate.reference_start, mate.reference_end)
                            # print mate.query_name, mate.reference_start , mate.reference_end
                            else:
                                addpos = (
                                pileupread.alignment.query_name, '-', pileupread.alignment.reference_name, '-',
                                pileupread.alignment.reference_start, pileupread.alignment.reference_end, 0, 0)
                            # print 'MATE',bam.mate(pileupread.alignment)
                            rqname_comp = rqname + '$' + pileupread.alignment.reference_name + '$' + str(
                                pileupcolumn.reference_pos + 1)
                            # addpos=(chr+'_'+str(pileupcolumn.reference_pos+1),pileupcolumn.reference_pos+1)
                            if not grdb.has_key(rqname):
                                # print rqname reference_start
                                outreads.write('>' + rqname_comp + '\n' + rseqname + '\n')
                            # grdb[rqname]=[addpos]
                            # else:
                            #	if addpos not in grdb[rqname]:
                            #		grdb[rqname].append(addpos)
                            if not grdb2.has_key(rname): grdb2[rname] = addpos
            if seq.strip() != '':
                if rmIndel:
                    # print 'Indels:',indels
                    if indels.count(1) > 0:
                        # print 'REMOVED by indels'
                        continue
                if rmOver:
                    over_ = getOverlap(rall)
                    seq, qual, squal, strand, blatc = normByOverlap(seq, strand, squal, blatc, 33, over_)
                # if over_.count(1)>0:
                #	print chr,pileupcolumn.reference_pos+1
                #	print seq
                #	print 'Over:',over_
                #	print blatc
                if blatr:
                    if testBlat(blatc):
                        seq, qual, squal, strand = normByBlat(seq, strand, squal, blatc, 33)
                    else:
                        continue
                mystrand = '2'
                # print seq,strand,strand.count('+'),strand.count('-')
                if uann and not getstrand:
                    if chr in tabix.contigs:
                        sres = [kk.strand for kk in tabix.fetch(reference=chr, start=(pileupcolumn.reference_pos),
                                                                end=(pileupcolumn.reference_pos + 1),
                                                                parser=pysam.asGTF())]
                        mystrand = vstrand(sres)
                if getstrand and not uann:
                    mystr = vstand(strand)
                    if mystr == '-':
                        mystrand = '0'
                    elif mystr == '+':
                        mystrand = '1'
                    else:
                        mystrand = '2'
                if mystrand == '0':
                    seq = comp(seq)
                    ref = comp(ref)
                if mystrand in ['0', '1'] and corrstr:
                    seq, qual, squal = normByStrand(seq, strand, squal, mystrand)
                cov, bcomp, subs, freq = BaseCount(seq, ref, mmf, vnuc)
                if cov < MINCOV: continue
                if exms and subs.count(' ') > 0: continue
                mqua = meanq(qual, len(seq))
                if expos:
                    if chr in extabix.contigs:
                        exres = [kk for kk in extabix.fetch(reference=chr, start=(pileupcolumn.reference_pos),
                                                            end=(pileupcolumn.reference_pos + 1))]
                        if len(exres) > 0: continue
                # se la sostituzione non e' in usubs
                if exinv and subs == '-': continue
                if not checkSubs(subs): continue
                # print out rna-seq info + dna-seq
                if gdic.has_key(pileupcolumn.reference_pos):  # abbiamo l'informazione genomica
                    if exnonh and not gdic[pileupcolumn.reference_pos][1]: continue
                    if mystrand == '0':
                        gdic[pileupcolumn.reference_pos][0][2] = comp2(eval(gdic[pileupcolumn.reference_pos][0][2]))
                        gdic[pileupcolumn.reference_pos][0][3] = comp(gdic[pileupcolumn.reference_pos][0][3])
                    line = '\t'.join(
                        [chr, str(pileupcolumn.reference_pos + 1), ref, mystrand, str(cov), mqua, str(bcomp), subs,
                         freq] + gdic[pileupcolumn.reference_pos][0]) + '\n'
                    out.write(line)
                else:
                    if exnosupp: continue
                    line = '\t'.join(
                        [chr, str(pileupcolumn.reference_pos + 1), ref, mystrand, str(cov), mqua, str(bcomp), subs,
                         freq] + ['-', '-', '-', '-', '-']) + '\n'
                    out.write(line)
                if slist: outrna.write('\t'.join([chr, str(pileupcolumn.reference_pos + 1), ref, seq, squal]) + '\n')
    bam.close()
    if not nogbam and isgbam: gbam.close()
    fasta.close()
    out.close()
    if uwf: wtabix.close()
    if uann: tabix.close()
    if expos: extabix.close()
    if slist:
        outrna.close()
        if not nogbam and isgbam: outdna.close()
    if os.path.getsize(outfile) == 0: os.remove(outfile)
    if greads: outreads.close()
    if addP:
        for Name in grdb2:
            pn = grdb2[Name]
            if pn[1] == '-':
                pcoo = [pn[4], pn[5]]
                pcoo.sort()
                outAddP.write('%s\t%i\t%i\n' % (pn[2], pcoo[0] - 100, pcoo[-1] + 100))
            else:
                if pn[0] == pn[1] and pn[2] == pn[3]:
                    pcoo = [xy for xy in pn[4:]]
                    pcoo.sort()
                    outAddP.write('%s\t%i\t%i\n' % (pn[2], pcoo[0] - 100, pcoo[-1] + 100))
        # outAddP.write('%s\t%s\n' %(Name,str(grdb2[Name])))
        outAddP.close()
    sys.stderr.write('Job completed for region: %s:%i-%i\n' % (chr, start_region + 1, lenregion))

def do_work(q):
    while True:
        x = q.get(block=True)
        if x == None: break
        exploreBAM(x)


####
wRegions = []
if fworkR:
    if NCPU == 1:
        wRegions.append((workR[0], workR[1][0] - 1, workR[1][1]))
    elif NCPU > 1:
        wlen = workR[1][1] - workR[1][0]
        wint = wlen / NCPU
        wstart = workR[1][0] - 1
        for i in range(NCPU - 1):
            wRegions.append((workR[0], wstart, wstart + wint))
            wstart = (wstart + wint)
        wRegions.append((workR[0], wstart, workR[1][1]))
####
work_queue = Queue()
suffix = []
kkn = 0
if fworkR:
    for i in wRegions:
        suff = '%s_%s_%i' % (i[0], outfolder_, kkn)
        suffix.append(suff)
        strinput = '$'.join([i[0], bamfile, str(i[1]), str(i[2]), suff])  # i+'$'+bamfile
        # print strinput
        work_queue.put(strinput)
        kkn += 1
else:
    for i in chrs:
        suff = '%s_%s_%i' % (i, outfolder_, kkn)
        suffix.append(suff)
        strinput = '$'.join([i, bamfile, '0', str(dicregions[i]), suff])  # i+'$'+bamfile
        # print strinput
        work_queue.put(strinput)
        kkn += 1
processes = []
for i in range(NCPU):
    work_queue.put(None)
    t = Process(target=do_work, args=(work_queue,))
    t.start()
    processes.append(t)
for t in processes:
    t.join()
work_queue.empty()
#
head = 'Region\tPosition\tReference\tStrand\tCoverage-q%i\tMeanQ\tBaseCount[A,C,G,T]\tAllSubs\tFrequency\tgCoverage-q%i\tgMeanQ\tgBaseCount[A,C,G,T]\tgAllSubs\tgFrequency\n' % (
MQUAL, gMQUAL)
sys.stderr.write('Merging Tables.\n')
if gziptab:
    o = gzip.open(outtable + '.gz', 'wb')
else:
    o = open(outtable, 'w')
if noheader == 0: o.write(head)
if slist:
    if gziptab:
        o2 = gzip.open(slistfile + '.gz', 'wb')
    else:
        o2 = open(slistfile, 'w')
    if len(gbamfile) != 0:
        if gziptab:
            o3 = gzip.open(gslistfile + '.gz', 'wb')
        else:
            o3 = open(gslistfile, 'w')
if greads:
    outReadsFile = os.path.join(outfolder, '%s_outReads' % (outfolder_))
    o4 = open(outReadsFile, 'w')
if addP:
    outPosFile = os.path.join(outfolder, '%s_outPosReads' % (outfolder_))
    o5 = open(outPosFile, 'w')
for i in suffix:
    if gziptab:
        tabfile = os.path.join(outfolder, 'table_%s.gz' % (i))
    else:
        tabfile = os.path.join(outfolder, 'table_%s' % (i))
    if os.path.exists(tabfile):
        if gziptab:
            f = gzip.open(tabfile, 'rb')
        else:
            f = open(tabfile)
        for j in f: o.write(j)
        f.close()
        os.remove(tabfile)
    if slist:
        if len(gbamfile) != 0:
            if gziptab:
                dnafile = os.path.join(outfolder, 'pileupDNA_%s.gz' % (i))
            else:
                dnafile = os.path.join(outfolder, 'pileupDNA_%s' % (i))
            if os.path.exists(dnafile):
                if gziptab:
                    f = gzip.open(dnafile, 'rb')
                else:
                    f = open(dnafile)
                for j in f: o3.write(j)
                f.close()
                os.remove(dnafile)
        if gziptab:
            rnafile = os.path.join(outfolder, 'pileupRNA_%s.gz' % (i))
        else:
            rnafile = os.path.join(outfolder, 'pileupRNA_%s' % (i))
        if os.path.exists(rnafile):
            if gziptab:
                f = gzip.open(rnafile, 'rb')
            else:
                f = open(rnafile)
            for j in f: o2.write(j)
            f.close()
            os.remove(rnafile)
    if greads:
        readsfile = os.path.join(outfolder, 'readsRNA_%s' % (i))
        if os.path.exists(readsfile):
            f = open(readsfile)
            for j in f: o4.write(j)
            f.close()
            os.remove(readsfile)
    if addP:
        addPfile = os.path.join(outfolder, 'readsPosRNA_%s' % (i))
        if os.path.exists(addPfile):
            f = open(addPfile)
            for j in f: o5.write(j)
            f.close()
            os.remove(addPfile)

o.close()
if slist:
    o2.close()
    if len(gbamfile) != 0: o3.close()
sys.stderr.write('Results saved on %s\n' % (outtable))
if slist:
    if len(gbamfile) != 0: sys.stderr.write('Pileup for DNA saved on %s\n' % (gslistfile))
    sys.stderr.write('Pileup for RNA saved on %s\n' % (slistfile))
if greads:
    o4.close()
    sys.stderr.write('RNA reads saved on %s\n' % (outReadsFile))
    fastqFiles = {'r1': os.path.join(outfolder, '%s_R1.fq' % (outfolder_)), 'r2': os.path.join(outfolder, '%s_R2.fq' % (outfolder_))}
    for i in range(1, len(fastq) + 1):
        sys.stderr.write('Getting reads R%i\n' % (i))
        cmd = '/opt/exp_soft/biomed/epicardi/bin/seqtk subseq %s %s > %s' % (
        fastq[i - 1], outReadsFile, fastqFiles['r' + str(i)])
        os.system(cmd)
    for i in range(1, len(fastq) + 1):
        sys.stderr.write('RNA reads in FASTQ saved on %s\n' % (fastqFiles['r' + str(i)]))
if addP: o5.close()
script_time = time.strftime("%d/%m/%Y %H:%M:%S", time.localtime(time.time()))
sys.stderr.write("Script time --> END: %s\n" % (script_time))
