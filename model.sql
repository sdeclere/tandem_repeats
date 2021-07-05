----------------------------------
-- TANDEM REP IN FAMS (DB) 
----------------------------------

--  stef 11/04/2018 (v1)
--

-- input file separator
.separator "\t"

-- TABLES
DROP TABLE IF EXISTS motifs;
CREATE TABLE motifs (
       sp     varchar(255) ,
       famid  varchar(255) ,
       gid    varchar(255) ,
       pos    varchar(255) ,
       nuseq  text,
	aaseq  text, 
       efid   varchar(255), 
	blamed BOOLEAN CHECK (blamed IN (0,1)), 
       annot  varchar(255));
.import test4 motifs

DROP TABLE IF EXISTS sp;
CREATE TABLE sp (
       spid      varchar(255) NOT NULL,
       vulgar    varchar(255) NOT NULL,
       clade     varchar(255) NOT NULL ) ;
.import sp.tab sp

DROP TABLE IF EXISTS pfam;
CREATE TABLE pfam (
       famid   varchar(255) NOT NULL,
       pfam    text,   
       annot   text,   
	blamed  BOOLEAN CHECK (blamed IN (0,1)) ) ;
.import pfam_ok_ok.txt pfam

-- DROP TABLE IF EXISTS fam_annots;
-- CREATE TABLE fam_annots (
--              famid      varchar(255) NOT NULL,
--              exonidx    float NOT NULL,
--              tag1       varchar(255) NOT NULL,
--              tag2       varchar(255) NOT NULL ) ;
--.import fam_annots_ok2.tab fam_annots

DROP TABLE IF EXISTS gene_fams; 
CREATE TABLE gene_fams (
       gid varchar(255) NOT NULL,
       fid varchar(255) NOT NULL );
.separator " "
.import maped_stable_ids_v2 gene_fams

-- Last but not least : init flags 

update pfam set blamed = 0;
update pfam set annot  = "NO_ANNOTS";
update motifs set annot  = "na";
update motifs set blamed = 0;

-- hard coded join
--update motifs set efid = (select fid from motifs as m, gene_fams as e where e.gid = m.gid);

-- in ens_data 
-- DROP TABLE IF EXISTS gene_fams2;
-- CREATE TABLE gene_fams2 (gid varchar(255) NOT NULL,fid varchar(255) NOT NULL);
-- .separator " "
-- .import /Volumes/PROJETS/EUTR/trwww/db3/maped_stable_ids_v3 gene_fams2
-- CREATE INDEX gene_fams2_idx  ON gene_fams2 (gid);