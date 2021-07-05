# -*- coding: utf-8 -*-

# parse hmmsearch results
# reify each modules, & dump module sequence in a db ready format 
## usage : cat table | python tbl2fasta.py <TAG>

import sys
import pickle
import os.path
import sqlite3
import slicer

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

from itertools import chain

# --------------------------------------------------------  CLASSES
class Module(object):
   def __init__(self, lst):
      self.target_name = lst[0]
      self.target_accession = lst[1]
      self.tlen = int(lst[2])
      self.query_name = lst[3]
      self.query_accession = lst[4]
      self.qlen = int(lst[5])
      self.E_value = float(lst[6])
      self.score  = float(lst[7])
      self.bias = float(lst[8])
      self.order = int(lst[9])
      self.total =  int(lst[10])
      self.c_Evalue  = float(lst[11])
      self.i_Evalue  = float(lst[12])
      self.dom_score =  float(lst[13])
      self.dom_bias  =  float(lst[14])
      self.hmm_from  =  int(lst[15])
      self.hmm_to  = int(lst[16])
      self.ali_from  = int(lst[17])
      self.ali_to  = int(lst[18])
      self.env_from  = int(lst[19])
      self.env_to  = int(lst[20])
      self.target_description = lst[21]

   def __str__(self):
      return '%s:[%d-%d]' %(self.target_name, self.env_from, self.env_to)

   # return the len of the detetcted module
   def cov(self):
      return (self.env_to - self.env_from)
      
# --------------------------------------------------------  GENE MODEL
class gene_model:
   
   def __init__(self, gid, tid, exons, seq):
      self.gid = gid
      self.tid = tid
      self.exons = exons
      self.seq = seq
      
   def __str__(self):
      return '%s\t%s\t%s\t%s...' % (self.gid, self.gid, self.exons, self.seq [0:10])
   
   def protein_len (self):
      ll = str2list(self.exons)
      return sum([y[1]-y[0] for y in ll])
      
   def __cmp__(self, o):
      if self.protein_len() < o.protein_len():
         return -1
      elif self.protein_len() > o.protein_len():
         return 1
      else:return 0
   
   def to_string(self):
      return [self.gid, self.tid, self.exons, self.seq]

def str2list(str_exons):
   splited_loc = [ x.split(':') for x in str_exons.split(',') ]
   return  [ ( int(x[0]), int(x[1]) ) for x in splited_loc ]

# --------------------------------------------------------  TBL parsing 
def group_by_id (lines):
   # extract ids from lines
   ids = set(map(lambda x:x[3], lines))
   # group by ids
   groups = [[y for y in lines if y[3]==x] for x in ids]
   return groups

def groups2dict(groups):
   # dict storing all modules for a key (=a prot id)
   fmods= {}
   # alloc objets
   for g in groups:
      m=[Module(x) for x in g]
      k = g[0][3]
      fmods[k] = m
   return fmods

# --------------------------------------------------------  DB RELATED
##MDB = "/Users/sdescorp/Desktop/FAMS/ens_data.db" 
##MDB="/Users/sdescorp/Desktop/BE_TR/FAMS.0/ens_data.db"
##MDB="/Users/sdescorp/Projets/gfr/BE_TR/FAMS.0/ens_data.db"
MDB="/Volumes/PROJETS/EUTR/trwww/db/ens_data.db"

def get_model_from_id(gid):
   ret = None 
   
   # connection to db
   conn = sqlite3.connect(MDB)
   
   try: 
      c = conn.cursor()
      c.execute("SELECT * FROM models WHERE gid = '%s';" % gid)
      ret =c.fetchall()
   except sqlite3.Error as er:
       print 'er:', er.message
   finally:
       conn.close()
   models = [ gene_model(*x) for x in ret ]
   models.sort()
   ##print [x.protein_len() for x in models]
   return models[-1]

# --------------------------------------------------------  STARTS HERE
if __name__ == '__main__':

      # read lines form sdtin skip lines starting with \#
      lines  = [x.strip().split() for x in sys.stdin.readlines() if not x.startswith('#')]
      groups = group_by_id(lines)
      mods   = groups2dict(groups)
      
      for k in mods:
          all=mods[k]
          print k,len(all)

      # sys.exit(1)
      #outf = open('%s_modules.fasta' % sys.argv[1], 'w')
      
      # loop over all genes 
      for k in mods: 
         #print >> sys.stderr, '+ searching  %s' % k
         gid,tid,sp,exons = (None,None,None,None)
         # loop over all modules 
         for i, m in enumerate(mods[k]):
            if (m.i_Evalue > 0.00001):
                #print "val is :"  + str(m.i_Evalue)
                continue

            gid,pid,sp = ( m.query_name.split('_')[0], m.query_name.split('_')[1].split('|')[0], m.query_name.split('|')[1] )
            
            # FIXME
            if sp.startswith('<SEED>'):
               sp = m.query_name.split('|')[2]
               
            # get sequence from db 
            try: 
               # get sequence 
               gm = get_model_from_id(gid)
               _seq=[]
               for e in gm.exons.strip().split(','):
                  start,end = ( int(x) for x in e.split(':') )
                  _seq.append( gm.seq[start:end] )
               cds_seq = "".join(_seq)
               ocds = Seq(cds_seq, generic_dna)
               # AA sequence 
               gene_seq = str( ocds.translate(to_stop=True) )
              
               # dump seq 
               if len(gene_seq[m.env_from:m.env_to])>5:
                  cds_start = slicer.p2c(m.env_from)[0]
                  cds_end = slicer.p2c(m.env_to)[1]
                  #print cds_end 
                  cds_frag_seq = ocds[cds_start:cds_end]
                  tag = sys.argv[1]
                  print "%s\t%s\t%s\t%d\t%s\t%s" % (sp, tag, gid, i, str(cds_frag_seq) ,gene_seq[m.env_from:m.env_to] )
               else: 
                  print >> sys.stderr, "!! big trouble sequence length of db and ortho are not consistant"
            except IndexError as e: 
               print >> sys.stderr, 'No model for : %s' % gid 
               continue
