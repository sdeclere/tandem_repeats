# To be used with : 
# for f in FAM_*fasta; do grep -H ">" $f ;done | tr ":" "\t " | sed -e 's/\.fasta//g' |tr '_' '\t' | awk '{print $1"_"$2"\t"$6}' > famid_prot.txt

import sys 

FAMS={}
with open (sys.argv[1]) as fd: 
    for line in fd:
        f,p = line.strip().split('\t')
        if f in FAMS:
            FAMS[f].update([p])
        else: 
            FAMS[f]=set([p])
for f in FAMS:
    content = FAMS[f] 
    print ("%s,%s" % (f, ",".join(content)))
