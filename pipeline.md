# Tandem repeats discovery pipeline 
This file schematically summarises the main stages of the analysis, explaining the command lines used. 

Version 9 as freezed on June the 5th. 

## 00 TREKS
java -Djava.io.tmpdir=${params.root}/tmp -jar $trex_jar -infile=$f -overlapfilter -nosplit > ${params.sp}_trex.run 
awk -f ${scripts}/treks_parser.awk ${trex_out} | awk '\$4 > 35' > ${params.sp}_trex.bed

## Make genes model's DB  
python scripts/dump_annotation.py ${gen_path} ${params.sp} > ens_models.txt 
python extract_msa.py ${deb_file} > ${tr_fasta} 

## run a script selecting only longest transcribes in ensembl annotation
find . -name "*pep.all.fa" -exec python scripts/keep_big.py {} \; > white_list_of_pep.txt  
while read line; do if [ -f TR0_${line}.fasta ]; then cp *$line.fasta NO_ALTER/; fi; done < white_list_of_pep.txt 

## Alignement of all motifs
for f in *.fasta; do t_coffee $f ; done 
mv *.aln T_COFFEE/
for f in *.aln; do t_coffee -other_pg seq_reformat -in $f -action +rm_gap 50 > ${f%.*}_50.aln ; done # col with more 50% of "-" 
for f in *_50.aln; do t_coffee -other_pg seq_reformat -in $f -output fasta_aln > ${f%.*}.fasta ; done # reformat to fasta 
for f in *_50.fasta ; do a=`echo $f |sed -e 's/_50//g'`; mv $f $a ; done 

## Generation of HMM
export HHLIB=/Users/sdescorp/projets/hhsuite-2.0.16/lib/hh/
for f in *.fasta; do hhsuite-2.0.16/bin/hhmake -M 50 -i $f -name $f; done 
mkdir HHM && mv *.hhm HHM/
cat *.hhm > db.modules

# Comparaison of HMMs
for f in *.hhm; do /Users/sdescorp/projets/hhsuite-2.0.16/bin/hhsearch -i $f -d modules.db ; done 

## Network & clustering 
for f in *.hhr; do python hhr2abc.py $f; done > all_vs_all.2
cat all_vs_all.2| awk '{if ($1 != $2) {print $0}}' | awk '{if ( ($4<1e-5) && ($6/$10 > 0.7) ) {print $0}}' > filter_105_07.txt
tr -d '%' < filter_105_07.txt | awk '$7 > 70 ' > filter_105_07_07.txt
cut -f1,2,4 filter_105_07_07.txt > filter_105_07_07.abc
mcxload -abc filter_105_07_07.abc  --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o filter_105_07_07.mci -write-tab filter_105_07_07.tab 
mcl filter_105_09_07.mci -I 1.7
mcxdump -icl out.filter_105_07_07.mci.I17 -o dump.filter_105_07_07.mci.I17 -tabr filter_105_07_07.tab 

## Families Construction 
var=0
while read p ; do var=$((var+1)); echo $p | awk -v idx=$var '{for (i=1; i<=NF; i++) {o = o $i " "}} END {print "cat " o " | tr -d \"-\" > FAM_" idx}'  ; done < dump.mcl

## FAM -> ORTHO_FAM 
for f in *.fasta ; do cat $f | python2.7  rest_one2one.py > ORTHOS/ORTHO_${f}; done 

## merge ortho_fams containing shared members 
for f in FAM_*fasta; do grep -H ">" $f ;done | tr ":" "\t " | awk -F '\t' -v OFS='\t' '{split($2,a,"_");print $1, a[1]}' | 
tr -d '>' | sed -e 's/\.fasta//g' > prot_map.txt 
python famid_prot2.py prot_map.txt > prot_map.csv
python connected_fams.py prot_map.csv >disinter.sh
 
## Remove duplicated sequences in merged fams 
for f in *.fasta ; do seqkit rmdup $f -o UNIQ/$f; done 

## final detection of tandem repeats with HMM 
for f in FAM_*.fasta; do 
    echo "hmmscan --domtblout ORTHOS/UNIQ/ORTHO_${f%.*}.tbl AFA/${f%.*}-90.hmm ORTHOS/UNIQ/ORTHO_${f}" 
; done

## build motifs db table
for f in *.tbl; do ff=$(echo $f |sed -e "s/\.tbl//") ; cat $f | python2.7 ../tbl2db.py $ff; done > motifs.table 




