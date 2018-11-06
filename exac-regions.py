from __future__ import print_function

# from UCSC. see data/get-chain.py, pipe output to sort -k1,1 -k2,2n | uniq | bgzip -c > data/self-chains.id90.bed.gz

SELF_CHAINS = "data/self-chains.id90.bed.gz"

# from UCSC.  data/segmental.bed.gz

SEGDUPS = "data/segmental.bed.gz"

import sys
import itertools as it
import operator

import numpy as np
import utils as u
from cyvcf2 import VCF
from pyfaidx import Fasta

import argparse

parser=argparse.ArgumentParser()
parser.add_argument("-n", "--nosingletons", help="if you do NOT want singletons", action="store_true", default=False)
parser.add_argument("-w", "--varflag", help="if you want separation by variant flags", action="store_true", default=False)
parser.add_argument("-x", "--variants", help="ExAC or some other such variant file (VCF.gz)") # again, -v prints a stupid doctest message
# ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3/ExAC.r0.3.sites.vep.vcf.gz "toyexac.vcf.gz"
parser.set_defaults(variants = 'data/gnomad-vep-vt.vcf.gz')
parser.add_argument("-e", "--exons", help="File of exons, or genome space in which you are interested (GTF.gz)")
# ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
parser.set_defaults(exons = 'data/Homo_sapiens.GRCh37.75.gtf.gz')
parser.add_argument("-c", "--coverage", help="Location of coverage files with {chrom} in name, or genome space in which you are interested (txt.gz)")
parser.add_argument("-d", "--depth", help="Coverage depth", default=10, type=int)
parser.add_argument("-l", "--limit", help="Coverage cutoff/limit", default=0.5, type=float) 
parser.add_argument("-f", "--fasta", help="Fasta file, or genome space in which you are interested (GTF.gz)")
# fasta combined from: ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/
#if [ ! -s hg19.fa ]; then
#    sed 's/^>chr/^>/g' $DATA/hg19.fasta > data/hg19.fa # needed for fetalcoords.py; really in GRCh37 format, but from hg19 origin
#fi
parser.set_defaults(fasta = 'data/grch37.fa')
# arguments in exclude must be in bed format and tabixed
parser.add_argument("-s", "--exclude", default=[], help="Exclude files, tabixed and in BED format; e.g. self-chains and segdups", nargs='*') 
# ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3/coverage
parser.set_defaults(coverage = 'data/exacv2.chr{chrom}.cov.txt.gz')
args=parser.parse_args()
nosingletons=args.nosingletons
varflag=args.varflag
VCF_PATH = args.variants
GTF_PATH = args.exons
COVERAGE_PATH = args.coverage
depth = args.depth
cutoff = args.limit
exclude = args.exclude
FASTA_PATH = args.fasta

zip = it.izip


exac = VCF(VCF_PATH)
kcsq = exac["CSQ"]["Description"].split(":")[1].strip(' "').split("|")

#exac = exac("2:112538945-112551053")


header = "chrom\tstart\tend\taf\tfunctional\tgene\ttranscript\texon\timpact\tvstart\tvend\tn_bases\tcg_content\tranges\tcoverage\tposns\tvarflag"
print("#" + header)
keys = header.split("\t")
global mranges, splitter


def merge_rows(rows):
    """
    >>> merge_rows([dict(gene='ABC', vstart=1, vend=3), dict(gene='ABC', vstart=2, vend=5), dict(gene='ABC', vstart=4, vend=6)])
    [{'vstart': 1, 'vend': 6, 'gene': 'ABC'}]
    >>> merge_rows([dict(gene='ABC', vstart=1, vend=6), dict(gene='ABC', vstart=2, vend=6), dict(gene='ABC', vstart=4, vend=6)])
    [{'vstart': 1, 'vend': 6, 'gene': 'ABC'}]
    >>> merge_rows([dict(gene='ABC', vstart=1, vend=2), dict(gene='ABC', vstart=1, vend=4), dict(gene='ABC', vstart=1, vend=6)])
    [{'vstart': 1, 'vend': 6, 'gene': 'ABC'}]
    >>> merge_rows([dict(gene='ABC', vstart=1, vend=2), dict(gene='Easy as 123', vstart=1, vend=4), dict(gene='Easy as 123', vstart=1, vend=6)])
    [{'vstart': 1, 'vend': 2, 'gene': 'ABC'}, {'vstart': 1, 'vend': 6, 'gene': 'Easy as 123'}]
    >>> merge_rows([dict(gene='ABC', vstart=1, vend=2), dict(gene='ABC', vstart=2, vend=4), dict(gene='ABC', vstart=4, vend=6)])
    [{'vstart': 1, 'vend': 2, 'gene': 'ABC'}, {'vstart': 2, 'vend': 4, 'gene': 'ABC'}, {'vstart': 4, 'vend': 6, 'gene': 'ABC'}]
    >>> merge_rows([dict(gene='ABC', vstart=1, vend=2), dict(gene='Easy as', vstart=1, vend=4), dict(gene='123', vstart=1, vend=6)])
    [{'vstart': 1, 'vend': 2, 'gene': 'ABC'}, {'vstart': 1, 'vend': 4, 'gene': 'Easy as'}, {'vstart': 1, 'vend': 6, 'gene': '123'}]
    >>> merge_rows([dict(gene='ABC', vstart=1, vend=2), dict(gene='ABC', vstart=1, vend=6), dict(gene='ABC', vstart=1, vend=4)])
    [{'vstart': 1, 'vend': 6, 'gene': 'ABC'}]
    """
    new_rows = [rows[0]]
    for row in rows[1:]:
        if new_rows[-1]['gene'] == row['gene'] and row['vstart'] < new_rows[-1]['vend'] and row['vend'] > new_rows[-1]['vend']:
            new_rows[-1]['vend'] = row['vend']
        elif new_rows[-1]['gene'] != row['gene'] or new_rows[-1]['vend'] <= row['vstart']:
            new_rows.append(row)
    return new_rows

def separate_ranges(ranges, varflags): # for putting each VARTRUE range in its own range list, so that coverage and "cg_content" are only for true regions
    """
    >>> separate_ranges([[(26782, 26890)], [(45349, 45487), (58320, 58416), (58687, 58777), (60611, 60618)]], [['VARFALSE'], ['VARFALSE', 'VARFALSE', 'VARFALSE', 'VARTRUE']])
    ([[(26782, 26890)], [(45349, 45487), (58320, 58416), (58687, 58777)], [(60611, 60618)]], [['VARFALSE'], ['VARFALSE', 'VARFALSE', 'VARFALSE'], ['VARTRUE']])
    >>> separate_ranges([[(26782, 26890)], [(45349, 45487), (58320, 58416), (58687, 58777), (60611, 60618), (60422, 60443)]], [['VARFALSE'], ['VARFALSE', 'VARFALSE', 'VARFALSE', 'VARTRUE', 'VARTRUE']])
    ([[(26782, 26890)], [(45349, 45487), (58320, 58416), (58687, 58777)], [(60611, 60618), (60422, 60443)]], [['VARFALSE'], ['VARFALSE', 'VARFALSE', 'VARFALSE'], ['VARTRUE', 'VARTRUE']])
    """
    newranges=[]; newvarflags=[]
    for rangelist, flaglist in zip(ranges, varflags):
        fr, fv, tr, tv = [], [], [], []
        for r, v in zip(rangelist, flaglist):
            if v=='VARFALSE':
                fr.append(r); fv.append(v) 
            if v=='VARTRUE':
                tr.append(r); tv.append(v)
        if fr and fv:
            newranges.append(fr)
            newvarflags.append(fv)
        if tr and tv:
            newranges.append(tr)
            newvarflags.append(tv)
  
    return newranges, newvarflags

import doctest
res = doctest.testmod(verbose=False)
if res.failed != 0:
    sys.exit(1)

chroms = [str(x) for x in range(1, 23)] + ["X", "Y"]

#for chrom, viter in it.groupby(exac, operator.attrgetter("CHROM")):

def checkac(info, idx):
    if isinstance(info['AC'], tuple):
        return info['AC_AFR'][idx] == 0 and info['AC_AMR'][idx] == 0 and info['AC_EAS'][idx] == 0 and info['AC_NFE'][idx] == 0 and info['AC_OTH'][idx] == 0 and info['AC_SAS'][idx] == 0
    else:
        return info['AC_AFR'] == 0 and info['AC_AMR'] == 0 and info['AC_EAS'] == 0 and info['AC_NFE'] == 0 and info['AC_OTH'] == 0 and info['AC_SAS'] == 0

def perchrom(vcf_chrom):
    vcf, chrom = vcf_chrom

    viter = VCF(VCF_PATH)(chrom)
    chrom=str(chrom)
    rows = []
    print("reading chrom " + chrom, file=sys.stderr)

    fasta = Fasta(FASTA_PATH, as_raw=True)
    fa = str(fasta[chrom])
    coverage_array = u.read_coverage(chrom, length=len(fa), cov=depth,
                        path=COVERAGE_PATH)

    gene_exon_starts, gene_exon_ends, splitters = u.read_exons("|tabix {gtf} {chrom}"
                                                            .format(chrom=chrom,
                                                                gtf=GTF_PATH),
                                                            chrom, cutoff,
                                                            coverage_array,
                                                            exclude)
                                            #"|tabix {bed} {chrom}".format(chrom=chrom, bed=SELF_CHAINS),"|tabix {bed} {chrom}".format(chrom=chrom, bed=SEGDUPS))

    prevpos=-1; idx=0
    for v in viter:
        if not (v.FILTER is None or v.FILTER in ["PASS", "SEGDUP", "LCR"]):
            continue
        info = v.INFO
        try:
            as_filter=info['AS_FilterStatus'].split(",")[0]
            if as_filter not in ["PASS", "SEGDUP", "LCR"] :
                continue
        except KeyError:
            pass
        try:
            csqs = [dict(zip(kcsq, c.split("|"))) for c in info['CSQ'].split(",")]
        except KeyError:
            continue
        # NOTE: using max here for alternates to be conservative
        try: # gnomad doesn't have adj like exacv1
            ac = info['AC_Adj']
        except KeyError:
            ac = info['AC']
        if not isinstance(ac, (int, long)):
            ac = max(ac)
        try:
            af = ac / float(info['AN_Adj'] or 1)
        except KeyError:
            af = ac / float(info['AN'] or 1)
        if ac == 1: #self-explanatory, but filters out singletons
            if nosingletons: continue
        # if checkac(info,idx):
            # if isinstance(info['AC_ASJ'], tuple):
                # if info['AC_ASJ'][idx] !=0:
                    # continue
            # else:
                # if info['AC_ASJ'] !=0:
                    # continue
            # if isinstance(info['AC_ASJ'], tuple):
                # if info['AC_FIN'][idx] !=0:
                    # continue
            # else:
                # if info['AC_FIN'] != 0:
                    # continue
            
        # NOTE: not requiring canonical or requiring the csq to match the
        # particular alt that we chose.
        for csq in (c for c in csqs if c['BIOTYPE'] == 'protein_coding'): # getting duplicate rows because of this, wastes memory and potentially compute time, could remove and replace with just if isfunctional, add to rows then move on?
            # skipping intronic
            if csq['Feature'] == '' or csq['EXON'] == '': continue #or csq['cDNA_position'] == '': continue
            if not u.isfunctional(csq): continue
            try:
                if csq['cDNA_position']:
                    cdna_start, cdna_end = u.get_cdna_start_end(csq['cDNA_position'], v)
            except KeyError:
                cdna_start, cdna_end = 'na','na' # apparently, sometimes ENSEMBL doesn't annotate splice_donor_variant&coding_sequence_variant combinations with cdna coords

            rows.append(dict(chrom=v.CHROM, vstart=v.start, vend=v.end, af=af,
                functional=int(u.isfunctional(csq)),
                gene=csq['SYMBOL'], transcript=csq['Feature'], exon=csq['EXON'],
                impact=csq['Consequence'],
                cdna_start=cdna_start,   cdna_end=cdna_end))

        

    # now we need to sort and then group by gene so we know the gaps.
    rows.sort(key=operator.itemgetter('gene', 'vstart', 'vend'))

#TODO:we're only using vstart. maybe we should be using end, normalizing and decomposing?
# we may be misrepresenting deletions as in the case of vstart, vend, pos: 145838622 145838695 145838623; shouldn't this whole area not be covered?
    
    rows = merge_rows(rows)

    out = []
    for chrom_gene, trows in it.groupby(rows, lambda row: (row['chrom'], row['gene'])):
        if not chrom_gene in gene_exon_starts: 
            sys.stderr.write("skipping:" + str(chrom_gene) + " as it was not pulled from GTF\n")
            continue
        exon_starts = gene_exon_starts[chrom_gene]
        exon_ends = gene_exon_ends[chrom_gene]
        last = exon_starts[0]
        splitter = splitters.get(chrom_gene, None)

        for i, row in enumerate(trows, start=1):
            # istart and iend determine if we need to span exons.

            assert row['vstart'] <= exon_ends[-1], (row, exon_ends) # maybe use POS instead of vstart, so we can normalize and decompose?; should i check if end is less?
            row['vstart']=row['vstart']+1 # vstart is bed format variant coordinate, still true maybe use POS instead of vstart?
            last2 = last
            if varflag:
                mranges, last, varflags = u.get_ranges(last, row['vstart'], row['vend'], exon_starts, exon_ends, row['chrom'])#TODO: fix get_ranges to do what split_ranges does, and land behind vend because it ends at vstart
            else:
                 mranges, last, varflags = u.get_ranges_w_variant(last, row['vstart'], row['vend'], exon_starts, exon_ends, row['chrom'])
            #print (last, row['vstart'], row['vend'], mranges, splitter, varflags, exon_starts, exon_ends)
            mranges2, varflags2 = u.split_ranges(mranges, splitter, varflags)
            if varflag:
                mranges2, varflags2 = separate_ranges(mranges2, varflags2)
            for ranges, vf in zip(mranges2, varflags2):
                row['coverage'] = ",".join(",".join(u.floatfmt(g) for g in coverage_array[s:e]) for s, e in ranges)
                row['posns'] = list(it.chain.from_iterable([range(s+1, e+1) for s, e in ranges])) # since range is not inclusive at the end add +1, need to add +1 to start
                row['ranges'] = ["%d-%d" % (s, e) for s, e in ranges]
                row['varflag'] = ",".join(vf)
                #print (row)
                #print (last,row['vstart'], "last n' vstart")
                #print (ranges, row['ranges'])
                #print (exon_starts, exon_ends, "starts n' ends")
                #print (splitter, "splitta!")
                seqs = [fa[s:e] for s, e in ranges]
                # this can happen for UTR variants since we can't really get
                # anything upstream of them.
                if row['posns'] == []:  # UTR:
                    p = row['vstart']
                    row['coverage'] = ",".join(u.floatfmt(g) for g in coverage_array[p:p+1])
                    row['posns'] = [p]

                # post-hoc sanity check
                exon_bases = set(it.chain.from_iterable(range(s, e) for s, e in zip(exon_starts, exon_ends)))
                ranges = set(it.chain.from_iterable(range(int(x[0]), int(x[1])) for x in (z.split("-") for z in row['ranges'])))
                m = len(ranges - exon_bases)
                if m > len(row['ranges']):
                    #print (last2, row['vstart'], row['vend'], exon_starts, exon_ends)
                    #print (row['ranges'], ranges - exon_bases)
                    print(last, row['chrom'], row['vstart'], row['vend'], row['ranges'], len(ranges -
                        exon_bases), mranges2, ranges, 
                        file=sys.stderr)

                # start or end? if we use end then can have - diff.
                row['ranges'] = ",".join(row['ranges'])
                row['n_bases'] = len(row['posns'])
                row['start'] = str(min(row['posns'])-1) #I put -1 because I am not including the position of the start coordinate, as it is 0-based.  however, I want to make it the start coordinate
                row['end'] = str(max(row['posns'])) # base on ranges
                row['posns'] = ",".join(map(str, row['posns']))
                row['cg_content'] = u.floatfmt(np.mean([u.cg_content(s) for s in seqs]))
                if row['cg_content'] == 'nan':
                    row['cg_content'] = '0'
                # we are re-using the dict for each loop so force a copy.
                try:
                    if row['ranges']:
                        endrange=int(row['ranges'].split('-')[-1])
                        if last < endrange:
                            last = endrange  #so we can start at where the last range ended
                        out.append(dict(row))
                except KeyError:
                    pass
            #else:
            #    row['ranges']=row['start']+"-"+row['end']
            #    row['end']=str(int(row['end'])+1)
            #    out.append(dict(row))
            # last = row['vstart']

        for exon_ending in exon_ends: # when variant is the last variant in a gene, it fails to make regions for the end of the gene, thus this block of code
            try:
                if mranges[-1][-1] > exon_ending: # do if mranges[-1][-1] < exon_ending or mranges is empty; if it is equal, we want it because it could be a 1 bp region
                    continue
            except IndexError:
                pass 
            if varflag:
                mranges, last, varflags = u.get_ranges(last, exon_ends[-1]+1, exon_ends[-1]+1, exon_starts, exon_ends, row['chrom']) #TODO: fix vend?
            else:
                 mranges, last, varflags = u.get_ranges_w_variant(last, exon_ends[-1]+1, exon_ends[-1]+1, exon_starts, exon_ends, row['chrom'])
            #print (last, row['vstart'], row['vend'], mranges, splitter, varflags, exon_starts, exon_ends)
            mranges2, varflags2 = u.split_ranges(mranges, splitter, varflags)
            if varflag:
                mranges2, varflags2 = separate_ranges(mranges2, varflags2)
            for ranges, vf in zip(mranges2, varflags2):
                row['coverage'] = ",".join(",".join(u.floatfmt(g) for g in coverage_array[s:e]) for s, e in ranges)
                row['posns'] = list(it.chain.from_iterable([range(s+1, e+1) for s, e in ranges])) #range is not inclusive at the end, need to add +1 to s
                row['ranges'] = ["%d-%d" % (s, e) for s, e in ranges]
                row['varflag'] = ",".join(vf)
                #print (row)
                #print (last,row['vstart'], "last n' vstart")
                #print (ranges, row['ranges'])
                #print (exon_starts, exon_ends, "starts n' ends")
                #print (splitter, "splitta!")
                seqs = [fa[s:e] for s, e in ranges]
                # this can happen for UTR variants since we can't really get
                # anything upstream of them.
                if row['posns'] == []:  # UTR:
                    p = row['vstart']
                    row['coverage'] = ",".join(u.floatfmt(g) for g in coverage_array[p:p+1])
                    row['posns'] = [p]
                # post-hoc sanity check
                exon_bases = set(it.chain.from_iterable(range(s, e) for s, e in zip(exon_starts, exon_ends)))
                ranges = set(it.chain.from_iterable(range(int(x[0]), int(x[1])) for x in (z.split("-") for z in row['ranges'])))
                m = len(ranges - exon_bases)
                if m > len(row['ranges']):
                    print(last, row['chrom'], row['vstart'], row['vend'], row['ranges'], len(ranges -
                        exon_bases), mranges2, ranges, 
                        file=sys.stderr)

                # start or end? if we use end then can have - diff.
                row['ranges'] = ",".join(row['ranges'])
                row['n_bases'] = len(row['posns'])
                row['start'] = str(min(row['posns'])-1) #I put -1 because I am not including the position of the start coordinate, as it is 0-based.  however, I want to make it he start coordinate 
                row['end'] = str(max(row['posns']))
                row['posns'] = ",".join(map(str, row['posns']))
                row['cg_content'] = u.floatfmt(np.mean([u.cg_content(s) for s in seqs]))
                if row['cg_content'] == 'nan':
                    row['cg_content'] = '0'
                # we are re-using the dict for each loop so force a copy.
                try:
                    if row['ranges']:
                        last = int(row['ranges'].split('-')[-1]) #so we can start at where the last range ended
                        out.append(dict(row))
                except KeyError:
                    pass

    # still print in sorted order
    out.sort(key=operator.itemgetter('start'))
    last = (None, None)
    outs = []
    for d in out:
        key = d['start'], d['end'], d['gene'] # added d['gene'] so it doesn't throw away longer regions without variants in a different gene overlapping the same genome space
        if key == last:
            continue
        last = key
        outs.append(d)
    return outs

from itertools import imap
import multiprocessing as mp
p = mp.Pool(22)

for outs in p.imap_unordered(perchrom, ((VCF, str(chrom)) for chrom in chroms)):
#for outs in imap(perchrom, ((VCF, str(chrom)) for chrom in chroms if chrom == '18')):
    for d in outs:
        print("\t".join(map(str, (d[k] for k in keys))))
