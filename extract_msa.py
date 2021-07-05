# split trex bed files, reformat correctly in order to re-aligned 
# build a uniq identifier for each motif 
import sys 
tr={}
prefix = sys.argv[2] # where to put all fasta files 

# reading bed file 
with open (sys.argv[1], 'r' ) as bedfile:
    for line in bedfile:
        if line.startswith('#'):
            continue 
        fld = line.strip().split()
        gid, aln = (fld[0], fld[5])
        if gid not in tr: 
            tr[gid] = [aln]
        else: 
            tr[gid].append(aln)
# dump each tr in a separate fasta file 
for gid in tr.keys():
    trs=tr[gid]
    for i,t in enumerate(trs):
        # for each TR in gene gid 
        outfn = 'TR%d_%s.fasta' % (i,gid)
        fd = open (prefix+"/" +outfn,'w')
        motifs = t.split(',')
        # for each motif in TR 
        for pos, m in enumerate(motifs):
            m = m.translate(None, '-')
            head = '>motif_%d_TR%d_%s\n' % (pos,i, gid)
            fd.write(head)
            fd.write(m+'\n')
        fd.close()
