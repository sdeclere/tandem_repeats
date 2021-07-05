#!/usr/bin/env python
# -*- coding: ascii -*-

import sys 

"""
extract E-value in HRR file. 
"""
## FUNC 
def do_parse(fd):
   count=0
   a,b = ["",""]
   desc = ""
   lines = open(fd,"r").readlines() 
   for l in lines : 
      if l.startswith('Match_columns'):
         nb_col = l.strip().split()[1]
      if l.startswith("Query"):
         a= l.strip().split()[1]
      if l.startswith(">"):
         count = count +1
         b= l.strip().split('>')[1]
      if l.startswith("Probab"):
         if count >= 2: 
            desc = l.strip().split()
            r = [x.split("=")[1] for x in desc]
            print a + "\t" + b + "\t" + "\t".join(r) + "\t" + nb_col
            b=""
            desc =""
         # elif count > 2:
         #    sys.exit(1)
         
if __name__=='__main__':
    do_parse(sys.argv[1])
    
