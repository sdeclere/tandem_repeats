# -*- coding: utf-8 -*-

# parse hmmsearch results
# reify each modules
# filter hits according distribution
# load fasta and check if
## usage : cat table | python parse_ortho_fam.py <fasta ortho fam>
## for f in *.tbl; do echo "cat $f | python parse_ortho_fam.py ${f%.*}.fasta"; done


import sys
import slicer
import pickle
import os.path
import sqlite3

from Bio import SeqIO
from itertools import chain

## usage :
## find ../dataset/  -name "*tbl"|while read f ; do sbatch --qos=fast --partition=dedicated,common --wrap="cat $f \
## | python mk_ortho_fam_table.py ${f%.*}.fasta"; done


# --------------------------------------------------------  DB RELATED
def get_exons_from_id(gid):
   ret = None 
   
   # connection to db
   conn = sqlite3.connect(MDB)
   try: 
      c = conn.cursor()
      c.execute("SELECT exons FROM models WHERE gid = '%s';" % gid)
      ret =c.fetchall()
   except sqlite3.Error as er:
       print 'er:', er.message
   finally:
       conn.close()
   return ret

def str2list(str_exons):
   splited_loc = [ x.split(':') for x in str_exons.split(',') ]
   return  [ ( int(x[0]), int(x[1]) ) for x in splited_loc ]


# --------------------------------------------------------  MISC MATH
# from https://stackoverflow.com/questions/15389768/standard-deviation-of-a-list
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
    pvar = ss/(n-1) # sample variance (for pop n)
    return pvar**0.5

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

# --------------------------------------------------------  FUN
def group_by_id (lines):
   # extract ids from lines
   ids = set(map(lambda x:x[0], lines))
   # group by ids
   groups = [[y for y in lines if y[0]==x] for x in ids]
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

# --------------------------------------------------------  STARTS HERE
##MDB = "/Users/sdescorp/Desktop/FAMS/ens_data.db"
MDB="/Users/sdescorp/Desktop/BE_TR/FAMS.0/ens_data.db"

if __name__ == '__main__':

      # read lines form sdtin skip lines starting with \#
      lines  = [x.strip().split() for x in sys.stdin.readlines() if not x.startswith('#')]
      groups = group_by_id(lines)
      mods   = groups2dict(groups)

      # create tabular module outputs
      #fout = sys.argv[1].split('.')[-2].split('/')[-1] # hacky
      #out = file(fout+".tab", 'w')
      #print >> sys.stderr, '_Starting analysis for : %s ' % fout
      
      for k in  mods.keys():
          print >> sys.stderr, '+ searching  %s' % k
          gid,tid,sp,exons = (None,None,None,None)
          
          for i, m in enumerate(mods[k]):
             gid,pid,sp = ( m.query_name.split('_')[0], m.query_name.split('_')[1].split('|')[0], m.query_name.split('|')[1] )

             if (m.i_Evalue > 0.00001):
                #print "val is :"  + str(m.i_Evalue)
                continue
               
             # FIXME
             if sp.startswith('<SEED>'):
                sp = m.query_name.split('|')[2]

             try:
                 # get gene model from db
                 str_lst = get_exons_from_id( gid ) # query returns a list of tuples
                 str_lst = [x[0] for x in str_lst]
                 all_exons = [str2list(x) for x in str_lst]             
                 all_exons.sort(key=lambda x: sum([y[1]-y[0] for y in x])) # keeps only the largest model
                 exons = all_exons[-1] # keeps only the largest model
             except IndexError:
                print >> sys.stderr, '!! No gene model avaible for %s' % gid
                continue

             # convert prot pos -> genomic pos
             rna = slicer.p2c(m.env_from)[1], slicer.p2c(m.env_to)[0]
             gstart, gend = (slicer.c2g(x, exons) for x in rna)

             # overlaping with exons
             ex_slice = slicer.find_overlaps(exons, (gstart,gend))
             str_ex_slice = [ '%d:%d' % (x[0],x[1]) for x in ex_slice[0]] if len(ex_slice[0])>0  else ['NA']
             str_ex_index = ','.join( [str(x) for x in ex_slice[1]] ) if len(ex_slice[1])>0 else 'NA'
             
             #print "%s\t%s\t%d\t%d\t%d\t%s"  % (gid, sp, i, m.env_from, m.env_to, ','.join(str_ex_slice) + '\t' + str_ex_index)
             mid_mot = ((m.env_to - m.env_from)/float(2)) + m.env_from
             r_pos = float(mid_mot)/m.qlen 
             #print gid, m.env_from, m.env_to, mid_mot, r_pos
             #print  (gid, sp, r_pos, m.env_from, m.env_to, gstart, gend, str_ex_index)
             print "%s\t%s\t%.3f\t%d\t%d\t%d\t%d\t%s"  % (gid, sp, r_pos, m.env_from, m.env_to, gstart, gend, str_ex_index)
             
