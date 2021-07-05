import networkx as nx 
import itertools
import sys 

FAMS={}

input = "fams_compo.csv"
with open(sys.argv[1]) as fd:

   for line in fd:
       l = line.strip().split(',')
       famid = l[0]
       gids= [ x.strip('"') for x in l[1:] ]

       # skip good fams
       #if len(gids)<=2 : continue
       FAMS[famid] = set(gids)

PROTS={}
for k in FAMS:
    for p in FAMS[k]:
        if p in PROTS:
            PROTS[p].append(k)
        else:
            PROTS[p]=[k]

# build GRAPH 
G = nx.Graph()
G.add_nodes_from(FAMS.keys()) # add nodes 

# add edges 
for p in PROTS.keys():
    for pair in itertools.product(PROTS[p], repeat=2):
        #print (pair)
        G.add_edge(*pair)

# emunarates connected components
for e in nx.connected_components(G):
    fname = [x+".fasta" for x in e]
    fname.sort()
    f0 = fname[0]

    #print(e)
    # discard concatanated seeds 
    print ( 'mv %s ORIG/' % ' '.join(fname) )

    # concat all seeds 
    print ( 'cat %s >| %s' % (' '.join(['ORIG/'+x for x in fname]), f0) )
    
    #sys.exit(1)