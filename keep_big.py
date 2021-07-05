# Load an ensembl multifasta peptides file, filer alt transcribe according to length    
# -> keep only the biggest one by gene id. 

from Bio import SeqIO
import sys

records = [ x for x in SeqIO.parse(sys.argv[1],"fasta") ]

seqs = {}
for r in records:
   gid_raw = [x for x in r.description.split() if x.startswith("gene:")]
   tid_raw = [x for x in r.description.split() if x.startswith("transcript:")]
   gid = gid_raw[0].split(':')[1] #.split('.')[0]
   tid = tid_raw[0].split(':')[1]

   print gid, tid 
   if gid in seqs:
      if len(r.seq) > seqs[gid][1] :
         seqs[gid] = (tid, len(r.seq))
   else: 
      seqs[gid] = (tid, len(r.seq))

for i in seqs.keys():
   print "%s_%s" % (i, seqs[i][0])
   