import sys
import pickle
import pdb

def p2c(ppos):
    ''' Return a tuple which contains (spliced)RNA projection of a protein location. '''
    lpos = ppos*3
    return (lpos,lpos+2)

def c2g(cpos, exons):
    ''' Return the a RNA location projected  on a gene'''
    d = 0
    for s,e in exons:
        l = e - s # len of the exon
        if cpos <= d+l: # test if the current exon hold the postion
            # cpos-d is the location of cpos into the current exon
            return s + cpos - d
        d += l
    return None

def over_ex(gstart,gend, exons):
    '''
    Return list of intervals overlapping gstart->gend.
    This fn suppose that inputed exons list is sorted.
    '''
    ret = []
    for s,e in exons:
        # if exons before first pos
        if e < gstart:
            continue
        # gstart and gend included in a single interval
        elif gstart >= s and gend < e:
            return [(gstart, gend)]
        # if gstart inculde in exon
        elif s <= gstart <= e:
            ret.append ( (gstart,e) )
        elif s <= gend <= e:
            ret.append( (s, gend) )
            return ret
        else:
            ret.append( (s,e) )

def find_overlaps(input_tuple_list, search_interval):
   """
   Collect all interval overlapping a bigger one.
   """
   # make sure that exons are sorted
   input_tuple_list_s = sorted(input_tuple_list, key=lambda tup: tup[0])
   over_exons = []
   index_exons = []
   for index, tup in enumerate(input_tuple_list_s):
       if overlaps(tup, search_interval):
           over_exons.append(tup)
           index_exons.append(index)
   return (over_exons, [x+1 for x in index_exons])

def overlaps (a, b):
    i1, i2  = sorted([a,b], key=lambda tup: tup[0])
    return ( i1[0] <= i2[0] <= i1[1] ) or (  i1[0] <= i2[1] <= i1[1] )
