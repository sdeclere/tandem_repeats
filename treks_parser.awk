# parse (poorly formated) T-Reks output
# produce a bed-like file 
BEGIN {
    OFS="\t";
    cutoff=5;
    print "#gene_id\tfrom\tto\tlength\texp";
}

# global swtich use to store alignement
/^\*.*$/ { inside_aln = 0; print atr, substr(aln,2); aln="" ; atr=""  }

inside_aln {
    aln = aln "," $0;
}

/^>.*$/ {
    atr="";
	split($4,a,":");
	split($5,b,":");
	id=a[2] "_" b[2];
}

/^Length.*$/ {
    atr =  id "\t" $8 "\t" $10 "\t" $2 "\t" $6;
    inside_aln=1; 
}



END { }
