package main

import (
	"database/sql"
	"fmt"
	"html/template"
	"io/ioutil"
	"log"
	"net/http"
	"os"
	"os/exec"
	"strings"
	"sort"

	_ "github.com/mattn/go-sqlite3"
)


var MODEL_DB_PATH = os.Getenv("MODEL_DB_PATH") /*  "/Users/sdescorp/projets/NO_ALTER/scripts/ens_models.db" */
var GEPARD_CP = os.Getenv("GEPARD_CP")         /* "/Users/sdescorp/projets/gepard/src" */

const (
	TR_DB_PATH        = "db3/motifs_v3.db"
	GEPARD_MAIN       = "org.gepard.client.cmdline.CommandLine"
	HTTP_DEFAULT_PORT = 8081
	SRVR_VERSION      = "devel3"
)

type Protein struct {
	Sp    string
	Gid   string
	Tid   string
	Pid   string
	AAseq string
}

type Gene struct {
	Gid    string
	Tid    string
	Fid		 string
	Exons  string
	Seq    string
	AAseq  string
	Blamed bool
}

type TRFam struct {
	Famid  string
	Gid    string
	Exons  string
	Sp     string
	Pos    string
	AaSeq  string
	NuSeq  string
	Vulgar string
	Clade  string
	Efamid  string
}

type IndexPage struct {
	Title string
	Fams  []FamPrototype
}

type IndexPage2 struct {
	Title string
	Fams  []FamNote
}

type FamPrototype struct {
	Famid  string
	Efamid  string
	Count  int
	SpList string
	Blamed bool
	Annot string
	//PFamList string
}

type FamNote struct {
	Famid  string
	EIndex float64
	Notes  string
}

func (g Protein) RunGepard() string {

	// path to current dir
	curdir, err := os.Getwd()
	if err != nil {
		log.Fatal(err)
	}

	// test if dotplot already exists
	if _, err := os.Stat(curdir + "/aa/" + g.Gid + ".png"); err == nil {
		log.Print(fmt.Sprintf("loading dotplot %s from cache", g.Gid))
		return "/aa/" + g.Gid + ".png"
	}

	// create tmp fasta file
	infile, err := ioutil.TempFile("/tmp", "proteic")
	if err != nil {
		log.Fatal(err)
	}
	log.Print("create fasta file for : " + g.Gid)
	defer os.Remove(infile.Name()) // clean up

	if _, err := infile.Write([]byte(">" + g.Gid + "\n" + g.AAseq)); err != nil {
		log.Fatal(err)
	}
	if err := infile.Close(); err != nil {
		log.Fatal(err)
	}

	// run dotplot
	var gepard_cmd = []string{
		"/usr/bin/java",
		"-cp " + GEPARD_CP, GEPARD_MAIN,
		"-seq1 " + infile.Name(),
		"-seq2 " + infile.Name(),
		"-matrix " + GEPARD_CP + "/matrices/blosum62.mat",
		"-outfile " + curdir + "/aa/" + g.Gid + ".png",
	}
	str_cmd := strings.Join(gepard_cmd, " ")
	log.Print("Gepard cmd is " + str_cmd)

	// do ewec
	cmd := exec.Command("bash", "-c", str_cmd)
	err = cmd.Run()
	if err != nil {
		log.Fatal(err)
	}

	return "/aa/" + g.Gid + ".png"
}

func (g Gene) RunGepard() string {

	// path to current dir
	curdir, err := os.Getwd()
	if err != nil {
		log.Fatal(err)
	}

	// test if dotplot already exists
	if _, err := os.Stat(curdir + "/cache/" + g.Gid + ".png"); err == nil {
		log.Print(fmt.Sprintf("loading dotplot %s from cache", g.Gid))
		return "/cache/" + g.Gid + ".png"
	}

	// create tmp fasta file
	infile, err := ioutil.TempFile("/tmp", "genomic")
	if err != nil {
		log.Fatal(err)
	}
	log.Print("create fasta file for : " + g.Gid)
	defer os.Remove(infile.Name()) // clean up

	if _, err := infile.Write([]byte(">" + g.Gid + "\n" + g.Seq)); err != nil {
		log.Fatal(err)
	}
	if err := infile.Close(); err != nil {
		log.Fatal(err)
	}

	// run dotplot
	var gepard_cmd = []string{
		"/usr/bin/java",
		"-cp " + GEPARD_CP, GEPARD_MAIN,
		"-seq1 " + infile.Name(),
		"-seq2 " + infile.Name(),
		"-matrix " + GEPARD_CP + "/matrices/edna.mat",
		"-exons " + g.Exons,
		"-outfile " + curdir + "/cache/" + g.Gid + ".png",
	}
	str_cmd := strings.Join(gepard_cmd, " ")
	log.Print("Gepard cmd is " + str_cmd)

	// do ewec
	cmd := exec.Command("bash", "-c", str_cmd)
	err = cmd.Run()
	if err != nil {
		log.Fatal(err)
	}

	return "/cache/" + g.Gid + ".png"
}

func GetFamBlamedStatus(str_gid string) bool {
	// open db
	db, err := sql.Open("sqlite3", TR_DB_PATH)
	if err != nil {
		log.Fatal("Failed to open db:", err)
	}
	defer db.Close()

	// query db
	//
	rows, err := db.Query("SELECT blamed FROM pfam  WHERE famid == \"" + str_gid + "\";")

	if err != nil {
		log.Fatal("Failed to call db.Query:", err)
	}
	defer rows.Close() // close at the end

	// retreive data & allocation of gene object
	var blamed bool
	for rows.Next() {
		err = rows.Scan(&blamed)

		if err != nil {
			log.Fatal("Failed to db.Query:", err)
		}
	}
	return blamed
}

func GetBlamedStatus(str_gid string) bool {
	// open db
	db, err := sql.Open("sqlite3", TR_DB_PATH)
	if err != nil {
		log.Fatal("Failed to open db:", err)
	}
	defer db.Close()

	// query db
	//
	rows, err := db.Query("SELECT blamed FROM motifs WHERE gid == \"" + str_gid + "\";")

	if err != nil {
		log.Fatal("Failed to call db.Query:", err)
	}
	defer rows.Close() // close at the end

	// retreive data & allocation of gene object
	var blamed bool
	for rows.Next() {
		err = rows.Scan(&blamed)

		if err != nil {
			log.Fatal("Failed to db.Query:", err)
		}
	}
	return blamed
}

// query model db, retunr a gene object
func QueryProtein(str_gid string) (Protein, error) {

	// open db
	db, err := sql.Open("sqlite3", MODEL_DB_PATH)
	if err != nil {
		log.Fatal("Failed to open db:", err)
	}
	defer db.Close()

	// query db
	rows, err := db.Query("SELECT sp, gid, tid, pid, aaseq FROM proteins WHERE gid between \"" + str_gid + "\" AND \"" + str_gid + "{\" ;")
	if err != nil {
		log.Fatal("Failed to call db.Query:", err)
	}
	defer rows.Close() // close at the end

	// retreive data & allocation of gene object
	var sp, gid, tid, pid, aaseq string
	var my_gene Protein

	for rows.Next() {
		err = rows.Scan(&sp, &gid, &tid, &pid, &aaseq)

		if err != nil {
			log.Fatal("Failed to db.Query:", err)
		}

		my_gene = Protein{Sp: gid, Gid: gid, Tid: tid, Pid: pid, AAseq: aaseq }
		log.Printf("__protein : %s ", my_gene)
	}

	return my_gene, err
}

// query model db, retunr a gene object
func QueryGene(str_gid string) (Gene, error) {
	log.Print(str_gid)

	// open db
	db, err := sql.Open("sqlite3", MODEL_DB_PATH)
	if err != nil {
		log.Fatal("Failed to open db:", err)
	}
	defer db.Close()

	// query db
	log.Println("SELECT models.gid, models.tid, gene_fams2.fid, models.exons, models.seq FROM models, gene_fams2  WHERE models.gid == gene_fams2.gid and models.gid == \"" + str_gid + "\";")
	rows, err := db.Query("SELECT models.gid, models.tid, gene_fams2.fid, models.exons, models.seq FROM models, gene_fams2  WHERE models.gid == gene_fams2.gid and models.gid == \"" + str_gid + "\";")
	if err != nil {
		log.Fatal("Failed to call db.Query:", err)
	}
	defer rows.Close() // close at the end

	// retreive data & allocation of gene object
	var gid, tid, fid, exons, seq string
	var my_gene Gene
	var blamed = GetBlamedStatus(str_gid)

	for rows.Next() {
		log.Println("Next row")
		err = rows.Scan(&gid, &tid, &fid, &exons, &seq)

		if err != nil {
			log.Fatal("Failed to db.Query:", err)
		}
		my_gene = Gene{Gid: gid, Tid: tid, Fid: fid, Exons: exons, Seq: seq, Blamed: blamed}
		break // ugly isn't it ? :'( 
	}

	// debug 
	//fmt.Printf("QG>%v",my_gene)
	

	return my_gene, err
}

// query model db, retunr a gene object
func QueryEnsFam(str_fam string, str_ens_fam string) ([]TRFam, error) {
	// open db
	db, err := sql.Open("sqlite3", TR_DB_PATH)
	if err != nil {
		log.Fatal("Failed to open db:", err)
	}
	defer db.Close()

	// query db
	//rows, err := db.Query("SELECT sp, gid, pos, nuseq, aaseq FROM motifs where famid == \"" + str_fam + "\";")
	//log.Print("SELECT motifs.sp,motifs.gid,motifs.pos,motifs.nuseq,motifs.aaseq,sp.vulgar,sp.clade FROM motifs,sp where motifs.famid == \"" + str_fam + "\";")
	rows, err := db.Query("SELECT m.sp,m.gid,m.pos,m.nuseq,m.aaseq,m.blamed, sp.vulgar, sp.clade, m.efid FROM motifs as m,sp WHERE m.famid == \"" + str_fam + "\" and m.sp == sp.spid;")

	if err != nil {
		log.Fatal("Failed to call db.Query:", err)
	}
	defer rows.Close() // close at the end

	// retreive data & allocation of gene object
	var sp, gid, pos, nuseq, aaseq, vulgar, clade, ens_famid string
	var blamed bool
	my_fams := []TRFam{}

	for rows.Next() {
		err = rows.Scan(&sp, &gid, &pos, &nuseq, &aaseq, &blamed, &vulgar, &clade, &ens_famid)

		if err != nil {
			log.Fatal("Failed to db.Query:", err)
		}
		my_fam := TRFam{Famid: str_fam, Gid: gid, Pos: pos, Sp: sp, AaSeq: aaseq, NuSeq: nuseq, Vulgar: vulgar, Clade: clade, Efamid: ens_famid}

		// revese append to blamed
		if (my_fam.Efamid == str_ens_fam && !blamed) {
			log.Printf("Filtered view : adding %v", my_fam)
			my_fams = append(my_fams, my_fam)
		}
	}

	return my_fams, err
}

// query model db, retunr a gene object
func QueryFam(str_fam string, rev bool) ([]TRFam, error) {
	// open db
	db, err := sql.Open("sqlite3", TR_DB_PATH)
	if err != nil {
		log.Fatal("Failed to open db:", err)
	}
	defer db.Close()

	
	// query db
	//rows, err := db.Query("SELECT sp, gid, pos, nuseq, aaseq FROM motifs where famid == \"" + str_fam + "\";")
	//log.Print("SELECT motifs.sp,motifs.gid,motifs.pos,motifs.nuseq,motifs.aaseq,sp.vulgar,sp.clade FROM motifs,sp where motifs.famid == \"" + str_fam + "\";")
	log.Println("SELECT m.sp,m.gid,m.pos,m.nuseq,m.aaseq,m.blamed, sp.vulgar,sp.clade FROM motifs as m,sp WHERE m.famid == \"" + str_fam + "\" and m.sp == sp.spid;")
	rows, err := db.Query("SELECT m.sp,m.gid,m.pos,m.nuseq,m.aaseq,m.blamed, sp.vulgar,sp.clade FROM motifs as m,sp WHERE m.famid == \"" + str_fam + "\" and m.sp == sp.spid;")

	if err != nil {
		log.Fatal("Failed to call db.Query:", err)
	}
	defer rows.Close() // close at the end

	// retreive data & allocation of gene object
	var sp, gid, pos, nuseq, aaseq, vulgar, clade string
	var blamed bool
	my_fams := []TRFam{}

	//log.Println("Start:	DB Q")
	for rows.Next() {
		//log.Println("insider")
		err = rows.Scan(&sp, &gid, &pos, &nuseq, &aaseq, &blamed, &vulgar, &clade)

		if err != nil {
			log.Fatal("Failed to db.Query:", err)
		}
		my_fam := TRFam{Famid: str_fam, Gid: gid, Pos: pos, Sp: sp, AaSeq: aaseq, NuSeq: nuseq, Vulgar: vulgar, Clade: clade}

		//log.Println(my_fam)

		// revese append to blamed
		if !rev {
				blamed = !blamed
		}

		if (!blamed) {
			my_fams = append(my_fams, my_fam)
			sort.SliceStable(my_fams, func(i, j int) bool {
    			return my_fams[i].Gid < my_fams[j].Gid
			})
		}
	}
	//log.Println("End:	DB Q")

	return my_fams, err
}

// query two different DB in order to
// guess number of different gene Fams
// It's ugly -> keeps kids away
// func QueryNbOrthoFams(str_famid string) (int){
//
// 	// open db
// 	db, err := sql.Open("sqlite3", TR_DB_PATH)
// 	if err != nil {
// 		log.Fatal("Failed to open db:", err)
// 	}
// 	defer db.Close()
//
// 	rows, err := db.Query("select distinct(gid) from motifs where famid = "  + str_famid + "\";")
// 	if err != nil {
// 		log.Fatal("Failed to call db.Query:", err)
// 	}
// 	defer rows.Close() // close at the end
//
// 	// TO finish
//
// 	return -1;
// }


// query model db, retunr a gene object
func QueryIndex(str_sel string) ([]FamPrototype, error) {
	// open db
	db, err := sql.Open("sqlite3", TR_DB_PATH)
	if err != nil {
		log.Fatal("Failed to open db:", err)
	}
	defer db.Close()

	// query db
	//  select famid, count(gid), count(sp), GROUP_CONCAT(distinct(sp))  from motifs group by(famid) limit 1;
	// rows, err := db.Query("select famid, count(distinct(gid)), GROUP_CONCAT(distinct(sp)) from motifs group by(famid);")
	rows, err := db.Query("select m.famid, count(distinct(gid)), count(distinct(m.efid)), pfam, p.blamed, p.annot from pfam as p, motifs as m where p.famid = m.famid group by m.famid order by count(distinct(gid)) desc")
	if err != nil {
		log.Fatal("Failed to call db.Query:", err)
	}
	defer rows.Close() // close at the end

	// retreive data & allocation of gene object
	var famid, sp string
	var gid_count int
	var blamedp bool
	var annot string
	var ens_fid string

	all_fams := []FamPrototype{}

	for rows.Next() {
		err = rows.Scan(&famid, &ens_fid, &gid_count, &sp, &blamedp, &annot)

		if err != nil {
			log.Fatal("Failed to db.Query:", err)
		}
		my_fam := FamPrototype{Famid: famid, Efamid: ens_fid, Count:gid_count , SpList: sp, Blamed: blamedp, Annot:annot}

		all_fams = append(all_fams, my_fam)
	}

	// Assoc famid against efam
	//m := make(map[string][]string)
	//for _,fam := range(all_fams){
	//	m[fam.Famid] = append(m[fam.Famid] , fam.Efamid)
	//}
	//log.Println(m)

	return all_fams, err
}

// query model db, retunr a gene object
func QueryIndex2(str_sel string) ([]FamPrototype, error) {
	// open db
	db, err := sql.Open("sqlite3", TR_DB_PATH)
	if err != nil {
		log.Fatal("Failed to open db:", err)
	}
	defer db.Close()

	// query db
	//  select famid, count(gid), count(sp), GROUP_CONCAT(distinct(sp))  from motifs group by(famid) limit 1;
	// rows, err := db.Query("select famid, count(distinct(gid)), GROUP_CONCAT(distinct(sp)) from motifs group by(famid);")
	rows, err := db.Query("select m.famid, m.efid, p.pfam, p.blamed, m.annot, count(distinct(m.gid)) from pfam as p, motifs as m where p.famid = m.famid group by m.famid, m.efid order by m.famid")
	if err != nil {
		log.Fatal("Failed to call db.Query:", err)
	}
	defer rows.Close() // close at the end

	// retreive data & allocation of gene object
	var famid, sp string
	var gid_count int
	var blamedp bool
	var annot string
	var ens_fid string

	all_fams := []FamPrototype{}

	for rows.Next() {
		err = rows.Scan(&famid, &ens_fid, &sp, &blamedp, &annot, &gid_count)

		if err != nil {
			log.Fatal("Failed to db.Query:", err)
		}
		my_fam := FamPrototype{Famid: famid, Efamid: ens_fid, Count: gid_count, SpList: sp, Blamed: blamedp, Annot:annot}
		
		//log.Println(my_fam)

		// skip zf pfam annot are in splist : arg ... 
		// skip singleton 
		if (gid_count > 1) && !strings.Contains(my_fam.SpList, "zf-") {
			all_fams = append(all_fams, my_fam)
		}
	}

	// Assoc famid against efam
	//m := make(map[string][]string)
	//for _,fam := range(all_fams){
	//	m[fam.Famid] = append(m[fam.Famid] , fam.Efamid)
	//}
	//log.Println(m)

	return all_fams, err
}

// query model db, retun list of familly
func indexHandler(w http.ResponseWriter, r *http.Request) {

	// get gid by trimming url
	str_select := strings.TrimPrefix(r.URL.Path, "/index/")

	index, _ := QueryIndex(str_select)

	my_index := IndexPage{Title: "Family Index ( All that glitters is not gold )", Fams: index}

	// new instance of template
	t, err := template.New("Index template").ParseFiles("tmpl/index.html")
	if err != nil {
		log.Fatalf("Template parsing: %s", err)
	}

	// La variable p sera réprésentée par le "." dans le layout
	err = t.ExecuteTemplate(w, "index", my_index)
	if err != nil {
		log.Fatalf("Template execution: %s", err)
	}
}

// Query & Render a gid
func gidHandler(w http.ResponseWriter, r *http.Request) {

	// if coming from post
	if r.Method == "POST" {
		// Form submitted
		r.ParseForm() // Required if you don't call r.FormValue()
		log.Print(fmt.Println(r.Form["curation_status"]))
	}

	// get gid by trimming url
	str_gid := strings.TrimPrefix(r.URL.Path, "/gid/")

	my_gene, _ := QueryGene(str_gid)
	my_prot, _ := QueryProtein(str_gid)

	my_gene.AAseq = my_prot.AAseq

	// dotplot
	_ = my_gene.RunGepard()

	// new instance of template
	t, err := template.New("Gene template").ParseFiles("tmpl/gene.html")
	if err != nil {
		log.Fatalf("Template parsing: %s", err)
	}

	// La variable p sera réprésentée par le "." dans le layout
	err = t.ExecuteTemplate(w, "gene", my_gene)
	if err != nil {
		log.Fatalf("Template execution: %s", err)
	}
}

type FamPage struct {
	Title    string
	IsBlamed bool
	Links    []string
	Genes    []string
	Clades   []string
	Sps      []string
	Fids     []string
	Ufids 	 []string
}

// remove duplicates entries in a string array
func unique(intSlice []string) []string {
    keys := make(map[string]bool)
    list := []string{}
    for _, entry := range intSlice {
        if _, value := keys[entry]; !value {
            keys[entry] = true
            list = append(list, entry)
        }
    }
    return list
}

// Query & Render a fam
func famHandler(w http.ResponseWriter, r *http.Request) {
	
	//log.Println("Start:	famHandler")
	// get gid by trimming url

	str_fam := strings.TrimPrefix(r.URL.Path, "/fam/")
	str_fam  = strings.TrimSuffix(str_fam, "/")

	str_fam_all := strings.Split(str_fam, "/")

	//log.Println("Start:	q")
	// do query db
	fam, _ := QueryFam(str_fam_all[0], true)
	//log.Println("End:	q")

	if len(str_fam_all)> 1{
		fam, _ = QueryFam(str_fam_all[0], false)
	}

	// map gid -> sp name
	genes := make(map[string]string)
	sp_dict := make(map[string]string)

	for _, e := range fam {
		genes[e.Gid] = e.Clade
		sp_dict[e.Gid] = e.Sp
	}

	// uniq gene == keys
	var keys []string
	for key, _ := range genes {
		keys = append(keys, key)
	}

	//log.Println("Start: dotplot	loop")
	// build dotplots
	var links, clades, sps, fids []string
	for _, k := range keys {
		g, _ := QueryGene(k)
		path := g.RunGepard()
		links = append(links, path)
		clades = append(clades, genes[k])
		sps = append(sps, sp_dict[k])
		fids = append(fids, g.Fid)
	}
	//log.Println("Start:	dotplot loop")

	isblamed := GetFamBlamedStatus(str_fam)

	page := FamPage{Title: str_fam, Links: links, Genes: keys, Clades: clades, Sps: sps, Fids: fids,  Ufids: unique(fids), IsBlamed: isblamed}

	// new instance of template
	t, err := template.New("Fam template").ParseFiles("tmpl/fam.html")
	if err != nil {
		log.Fatalf("Template parsing: %s", err)
	}

	// La variable p sera réprésentée par le "." dans le layout
	err = t.ExecuteTemplate(w, "fam", page)
	if err != nil {
		log.Fatalf("Template execution: %s", err)
	}
	//log.Println("End:	famHandler")
}

// Query & Render a fam
func xfamHandler(w http.ResponseWriter, r *http.Request) {
	// get gid by trimming url
	str_fam := strings.TrimPrefix(r.URL.Path, "/xfam/")
	str_both := strings.Split(str_fam, "/")
	log.Printf("New fam view : %v", str_both)

	// do query db
	fam, _ := QueryEnsFam(str_both[0], str_both[1])

	// map gid -> sp name
	genes := make(map[string]string)
	sp_dict := make(map[string]string)

	for _, e := range fam {
		genes[e.Gid] = e.Clade
		sp_dict[e.Gid] = e.Sp
	}

	// uniq gene == keys
	var keys []string
	for key, _ := range genes {
		keys = append(keys, key)
	}

	// build dotplots
	var links, clades, sps, fids []string
	for _, k := range keys {
		g, _ := QueryGene(k)
		path := g.RunGepard()
		links = append(links, path)
		clades = append(clades, genes[k])
		sps = append(sps, sp_dict[k])
		fids = append(fids, g.Fid)
	}

	isblamed := GetFamBlamedStatus(str_fam)

	page := FamPage{Title: str_fam, Links: links, Genes: keys, Clades: clades, Sps: sps, Fids: fids,  Ufids: unique(fids), IsBlamed: isblamed}

	// new instance of template
	t, err := template.New("Fam template").ParseFiles("tmpl/fam.html")
	if err != nil {
		log.Fatalf("Template parsing: %s", err)
	}

	// La variable p sera réprésentée par le "." dans le layout
	err = t.ExecuteTemplate(w, "fam", page)
	if err != nil {
		log.Fatalf("Template execution: %s", err)
	}
}

// Query & Render a fam (prot's way)
func pamHandler(w http.ResponseWriter, r *http.Request) {
	// get gid by trimming url
	str_fam := strings.TrimPrefix(r.URL.Path, "/pam/")
	// do query db
	fam, _ := QueryFam(str_fam, true)

	// map gid -> sp name
	genes := make(map[string]string)
	sp_dict := make(map[string]string)

	for _, e := range fam {
		genes[e.Gid] = e.Clade
		sp_dict[e.Gid] = e.Sp
	}

	// uniq gene == keys
	var keys []string
	for key, _ := range genes {
		keys = append(keys, key)
	}

	// build dotplots
	var links, clades, sps []string
	for _, k := range keys {
		g, _ := QueryProtein(k)

		// no protein entry in db ?
		if g.Gid == "" {
			continue
		}

		path := g.RunGepard()
		links = append(links, path)
		clades = append(clades, genes[k])
		sps = append(sps, sp_dict[k])
	}

	page := FamPage{Title: str_fam, Links: links, Genes: keys, Clades: clades, Sps: sps}

	// new instance of template
	t, err := template.New("xFam template").ParseFiles("tmpl/xfam.html")
	if err != nil {
		log.Fatalf("Template parsing: %s", err)
	}

	// La variable p sera réprésentée par le "." dans le layout
	err = t.ExecuteTemplate(w, "fam", page)
	if err != nil {
		log.Fatalf("Template execution: %s", err)
	}
}

// query model db, retun list of annoted familly
func index2Handler(w http.ResponseWriter, r *http.Request) {

	//debug 
	//cpuProfile, _ := os.Create("cpuprofile")
	//pprof.StartCPUProfile(cpuProfile)

	// get gid by trimming url
	str_select := strings.TrimPrefix(r.URL.Path, "/index2/")

	index, _ := QueryIndex2(str_select)

	my_index := IndexPage{Title: "Megasat bestiary", Fams: index}

	// new instance of template
	t, err := template.New("Index template").ParseFiles("tmpl/index3.html")
	if err != nil {
		log.Fatalf("Template parsing: %s", err)
	}

	// La variable p sera réprésentée par le "." dans le layout
	err = t.ExecuteTemplate(w, "index", my_index)
	if err != nil {
		log.Fatalf("Template execution: %s", err)
	}

	//pprof.StopCPUProfile()
}

func blameHandler(w http.ResponseWriter, r *http.Request) {

	// get gid by trimming url remove trailling /
	str_gid := strings.TrimPrefix(r.URL.Path, "/blame/")
	str_gid = strings.TrimSuffix(str_gid, "/")

	str_fam_gid := strings.Split(str_gid, "/")
	log.Print(fmt.Sprintf("do blame: %s", str_gid))

	// open db
	db, err := sql.Open("sqlite3", TR_DB_PATH)
	if err != nil {
		log.Fatal("Failed to open db:", err)
	}
	defer db.Close()

	// update db
	//

	blame_status := SwitchBlamedStatus(str_fam_gid[1])
	sql_update := "UPDATE motifs SET blamed =" + blame_status + " WHERE gid == \"" + str_fam_gid[1] + "\";"

	if (len(str_fam_gid)>2) {
		blame_status = SwitchBlamedStatus(str_fam_gid[2])
		sql_update = "UPDATE motifs SET blamed =" + blame_status + " WHERE gid == \"" + str_fam_gid[2] + "\";"
	}

	log.Print(sql_update)
	_, uerr := db.Exec(sql_update)
	if uerr != nil {
		log.Fatal("Failed to update db:", uerr)
	}

	log.Print(fmt.Sprintf("Redirect to : %s", "/fam/"+str_fam_gid[0]))

	if len(str_fam_gid)>2{
		//http.Redirect(w, r, "/xfam/"+str_fam_gid[0]+"/"+str_fam_gid[1], http.StatusFound)
		http.Redirect(w, r, "/fam/"+str_fam_gid[0], http.StatusFound)
	}else {
		http.Redirect(w, r, "/fam/"+str_fam_gid[0], http.StatusFound)
	}
}

func SwitchBlamedStatus(str_gid string) string{
	// open db
	db, err := sql.Open("sqlite3", TR_DB_PATH)
	if err != nil {
		log.Fatal("Failed to open db:", err)
	}
	defer db.Close()

	// query db
	//
	rows, err := db.Query("SELECT blamed FROM motifs WHERE gid == \"" + str_gid + "\";")

	if err != nil {
		log.Fatal("Failed to call db.Query:", err)
	}
	defer rows.Close() // close at the end

	// retreive data & allocation of gene object
	var blamed bool
	for rows.Next() {
		err = rows.Scan(&blamed)

		if err != nil {
			log.Fatal("Failed to db.Query:", err)
		}
	}
	if blamed{
		return "0"
	}
	return "1"
}

func blessHandler(w http.ResponseWriter, r *http.Request) {

	// get gid by trimming url remove trailling /
	str_gid := strings.TrimPrefix(r.URL.Path, "/bless/")
	str_gid = strings.TrimSuffix(str_gid, "/")

	str_fam_gid := strings.Split(str_gid, "/")
	log.Print(fmt.Sprintf("do blame: %s", str_gid))

	// open db
	db, err := sql.Open("sqlite3", TR_DB_PATH)
	if err != nil {
		log.Fatal("Failed to open db:", err)
	}
	defer db.Close()

	// update db
	//
	sql_update := "UPDATE motifs SET blamed =0 WHERE gid == \"" + str_fam_gid[1] + "\";"

	if (len(str_fam_gid)>2) {
		sql_update = "UPDATE motifs SET blamed = 1 WHERE gid == \"" + str_fam_gid[2] + "\";"
	}

	log.Print(sql_update)
	_, uerr := db.Exec(sql_update)
	if uerr != nil {
		log.Fatal("Failed to update db:", uerr)
	}

	if len(str_fam_gid)>2{
		http.Redirect(w, r, "/fam/"+str_fam_gid[0]+"/"+str_fam_gid[1], http.StatusFound)
	}else {
		http.Redirect(w, r, "/fam/"+str_fam_gid[0], http.StatusFound)
	}

}

func fblameHandler(w http.ResponseWriter, r *http.Request) {

	// get gid by trimming url
	str_gid := strings.TrimPrefix(r.URL.Path, "/fblame/")
	str_fam_annot := strings.Split(str_gid, "/")
	log.Print(fmt.Sprintf("do blame: %s", str_fam_annot))

	// open db
	db, err := sql.Open("sqlite3", TR_DB_PATH)
	if err != nil {
		log.Fatal("Failed to open db:", err)
	}
	defer db.Close()

	// update db
	//
	sql_update := "UPDATE pfam SET blamed = 1 WHERE famid == \"" + str_fam_annot[0] + "\";"
	log.Print(sql_update)

	if len(str_fam_annot) > 1 {
					 sql_update = "UPDATE pfam SET blamed = 1  WHERE famid IN (SELECT p.famid  FROM pfam as p, motifs as m WHERE m.efid== \"" + str_fam_annot[1] + "\" AND p.famid == m.famid AND p.famid == \"" + str_fam_annot[0] + "\");"
					 log.Print(fmt.Sprintf("do annot: %s",sql_update))
	 }

	_, uerr := db.Exec(sql_update)
	if uerr != nil {
		log.Fatal("Failed to update db:", uerr)
	}

	http.Redirect(w, r, "/index2/", http.StatusFound)
}

func fblessHandler(w http.ResponseWriter, r *http.Request) {

	// get gid by trimming url
	str_gid := strings.TrimPrefix(r.URL.Path, "/fbless/")
	str_fam_annot := strings.Split(str_gid, "/")
	log.Print(fmt.Sprintf("do blame: %s", str_fam_annot))


	// open db
	db, err := sql.Open("sqlite3", TR_DB_PATH)
	if err != nil {
		log.Fatal("Failed to open db:", err)
	}
	defer db.Close()

	// update db
	//
	sql_update := "UPDATE pfam SET blamed = 0 WHERE famid == \"" + str_gid + "\";"
	if len(str_fam_annot) > 1 {
					 sql_update = "UPDATE pfam SET blamed = 0 WHERE famid IN (SELECT p.famid  FROM pfam as p, motifs as m WHERE m.efid== \"" + str_fam_annot[1] + "\" AND p.famid == m.famid AND p.famid == \"" + str_fam_annot[0] + "\");"
					 log.Print(fmt.Sprintf("do annot: %s",sql_update))
	 }

	log.Print(sql_update)
	_, uerr := db.Exec(sql_update)
	if uerr != nil {
		log.Fatal("Failed to update db:", uerr)
	}

	http.Redirect(w, r, "/index2/", http.StatusFound)
}

func fannotHandler(w http.ResponseWriter, r *http.Request) {

	// get gid by trimming url
	str_gid := strings.TrimPrefix(r.URL.Path, "/seta/")
	str_fam_annot := strings.Split(str_gid, "/")


	log.Print(fmt.Sprintf("do annot: %s", str_fam_annot))
	// open db
	db, err := sql.Open("sqlite3", TR_DB_PATH)
	if err != nil {
		log.Fatal("Failed to open db:", err)
	}
	defer db.Close()

	// update db
	//
	sql_update := "UPDATE motifs SET annot = " + "\"" + str_fam_annot[1] + "\"" + " WHERE famid == \"" + str_fam_annot[0] + "\";"

	if len(str_fam_annot) > 2 {
		sql_update = "UPDATE motifs SET annot = " + "\"" + str_fam_annot[2] + "\"" + " WHERE famid == \"" + str_fam_annot[0] + "\"" + "AND fid == \"" + str_fam_annot[1] + "\";"

		//sql_update = "UPDATE pfam SET annot = \"" + str_fam_annot[2] + "\"" + " WHERE gid IN (SELECT m.gid  FROM pfam as p, motifs as m WHERE m.fid== \"" + str_fam_annot[1] + "\" AND p.famid == m.famid AND p.famid == \"" + str_fam_annot[0] + "\");"
		//sql_update = "UPDATE pfam SET annot = " + "\"" + str_fam_annot[2] + "\"" + " WHERE EXISTS (SELECT * FROM pfam as p, motifs as m WHERE m.fid== \"" + str_fam_annot[1] + "\" AND p.famid = m.famid AND p.famid == \"" + str_fam_annot[0] + ")\";"
		log.Print(fmt.Sprintf("do annot: %s",sql_update))
	}

	log.Print(sql_update)
	_, uerr := db.Exec(sql_update)
	if uerr != nil {
		log.Fatal("Failed to update db:", uerr)
	}

	http.Redirect(w, r, "/index2/", http.StatusFound)
}

// log func
func Log(handler http.Handler) http.Handler {
	return http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		fmt.Printf("%s %s %s\n", r.RemoteAddr, r.Method, r.URL)
		handler.ServeHTTP(w, r)
	})
}

func main() {

	// http server
	log.Print("Starting web-tr")
	log.Print("Version " + fmt.Sprintf(":%s", SRVR_VERSION))
	log.Print("point your browser on http://localhost:8080/index2")

	/* curation */
	http.HandleFunc("/blame/", blameHandler)
	http.HandleFunc("/bless/", blessHandler)

	http.HandleFunc("/fblame/", fblameHandler)
	http.HandleFunc("/fbless/", fblessHandler)

	http.HandleFunc("/seta/", fannotHandler) // for annotations

	/* url handlers */
	http.HandleFunc("/index2/", index2Handler)

	/* url handlers */
	http.HandleFunc("/index/", indexHandler)

	/* url handlers */
	http.HandleFunc("/gid/", gidHandler)

	/* url handlers */
	http.HandleFunc("/fam/", famHandler)

	/* url handlers */
	http.HandleFunc("/xfam/", xfamHandler)

	/* url handlers */
	http.HandleFunc("/pam/", pamHandler)

	/* Static files handlers : js, css, etc. */
	http.Handle("/static/", http.StripPrefix("/static/", http.FileServer(http.Dir("static"))))

	/* dotplots  */
	http.Handle("/cache/", http.StripPrefix("/cache/", http.FileServer(http.Dir("cache"))))
	http.Handle("/aa/", http.StripPrefix("/aa/", http.FileServer(http.Dir("aa"))))

	/* finally starts web server */
	log.Print(fmt.Sprintf("web-tr up & running on port: %d", HTTP_DEFAULT_PORT))
	log.Fatal(http.ListenAndServe(fmt.Sprintf(":%d", HTTP_DEFAULT_PORT), Log(http.DefaultServeMux)))
}
