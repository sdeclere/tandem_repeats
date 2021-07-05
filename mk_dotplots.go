package main

import (
	"database/sql"
	"fmt"
	"io/ioutil"
	"log"
	"os"
	"os/exec"
	"strings"

	_ "github.com/mattn/go-sqlite3"
)

var MODEL_DB_PATH = os.Getenv("MODEL_DB_PATH") /*  "/Users/sdescorp/projets/NO_ALTER/scripts/ens_models.db" */
var GEPARD_CP = os.Getenv("GEPARD_CP")         /*  "/Users/sdescorp/projets/gepard/src" */

const (
	TR_DB_PATH        = "db2/motifs_v2.db"
	GEPARD_MAIN       = "org.gepard.client.cmdline.CommandLine"
	HTTP_DEFAULT_PORT = 8080
	SRVR_VERSION      = 0.1
)

type Gene struct {
	Gid    string
	Tid    string
	Exons  string
	Seq    string
	AAseq  string
	Blamed bool
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

// query model db, retunr a gene object
func QueryGene(str_gid string) (Gene, error) {

	// open db
	db, err := sql.Open("sqlite3", MODEL_DB_PATH)
	if err != nil {
		log.Fatal("Failed to open db:", err)
	}
	defer db.Close()

	// query db
	rows, err := db.Query("SELECT gid, tid, exons, seq FROM models WHERE gid == \"" + str_gid + "\";")
	if err != nil {
		log.Fatal("Failed to call db.Query:", err)
	}
	defer rows.Close() // close at the end

	// retreive data & allocation of gene object
	var gid, tid, exons, seq string
	var my_gene Gene
	var blamed = false

	for rows.Next() {
		err = rows.Scan(&gid, &tid, &exons, &seq)

		if err != nil {
			log.Fatal("Failed to db.Query:", err)
		}
		my_gene = Gene{Gid: gid, Tid: tid, Exons: exons, Seq: seq, Blamed: blamed}
	}

	return my_gene, err
}

func main() {
	// open db
	db, err := sql.Open("sqlite3", TR_DB_PATH)
	if err != nil {
		log.Fatal("Failed to open db:", err)
	}
	defer db.Close()

	// query db
	rows, err := db.Query("SELECT distinct(gid) FROM motifs ;")

	if err != nil {
		log.Fatal("Failed to call db.Query:", err)
	}
	defer rows.Close() // close at the end

	var gid string
	for rows.Next() {
		err = rows.Scan(&gid)

		if err != nil {
			log.Fatal("Failed to db.Query:", err)
		}

		g, _ := QueryGene(gid)
		path := g.RunGepard()
		log.Print("creating dotplot for : " + path)
	}
}
