# serialization of dict containing genbank annotations 
# for a gene -> list on exons -> sequences  
# script args : folder where gbk are located, prefix 

from Bio import SeqIO, SeqFeature, Seq
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import sys
import glob
import re
import gzip



# -- func
def mk_uniq_id(feat):
    '''
    (Hopefully) return an uniq id resulting by the concatanation of gene and transcript_id feilds -- pretty ugly isn't it ?
    '''
    gg = feat.qualifiers['gene'][0].split('.')[0]
    tt = feat.qualifiers['note'][0].split("=")[1].split('.')[0] # why transcript_id is hidden in note ?
    return "%s\t%s" % (gg,tt)

def format_exons(l):
    ex = [ "{}:{}".format(x[0],x[1]) for x in l]
    return ','.join(ex)

# iter on all genbanks file
gbanks = glob.glob('%s/*.dat.gz' % sys.argv[1])

for f in gbanks :
   print >> sys.stderr, 'Alloc : %s ' % f
   handle = gzip.open(f, "r")
   gbank=SeqIO.parse(handle,'genbank')
   cid = re.search(r".*\.chromosome\.(.+)\..*", f, re.IGNORECASE)

   for k in gbank:
      for feat in k.features:
         if(feat.type =="CDS"):
            gid = mk_uniq_id(feat)
            start=feat.location.start.position
            end=feat.location.end.position
            if len(feat.sub_features) == 0: # if gene have no intron mk a big exon
               exons = [(start, end)]
               exons = [(x[0]-start,x[1]-start) for x in exons ]
               # BUG : print twice 
               #print '%s\t%s\t%s' % ( gid, format_exons(exons), k.seq[start:end] )
            else:
               exons = [ (sf.location.start.position,sf.location.end.position) for sf in feat.sub_features ]
               # change location according the start of the CDS
               exons = [ (x[0]-start,x[1]-start) for x in exons ]
            if feat.strand == 1:
               #print '%s\t%s\t%s' % (gid, exons, k.seq[start:end])
               print '%s\t%s\t%s' % ( gid, format_exons(exons), k.seq[start:end] )
            else:
               l = feat.location.end.position - feat.location.start.position
               exons = [(l-x[1],l-x[0]) for x in exons ] # flip the locations
               #print '%s\t%s\t%s' % (gid, exons, k.seq[start:end].reverse_complement())
               print '%s\t%s\t%s' % ( gid, format_exons(exons), k.seq[start:end].reverse_complement() )


