# parse seed fasta file 
# 1 - create a list of uniq geneid 
# 2 - request ensembl in order to gather all homolgous proteins 
# 3 - dump ortho families 
## usage : script <SEED:fasta formated> 


import requests, sys
import time 

# genomes list 
str_genomes = "ailuropoda_melanoleuca,anas_platyrhynchos,anolis_carolinensis,astyanax_mexicanus,bos_taurus,callithrix_jacchus,canis_familiaris,cavia_porcellus,chlorocebus_sabaeus,choloepus_hoffmanni,dasypus_novemcinctus,dipodomys_ordii,echinops_telfairi,equus_caballus,erinaceus_europaeus,felis_catus,ficedula_albicollis,gadus_morhua,gallus_gallus,gasterosteus_aculeatus,gorilla_gorilla,homo_sapiens,ictidomys_tridecemlineatus,latimeria_chalumnae,lepisosteus_oculatus,loxodonta_africana,macaca_mulatta,macropus_eugenii,meleagris_gallopavo,microcebus_murinus,monodelphis_domestica,mus_musculus,mustela_putorius_furo,myotis_lucifugus,nomascus_leucogenys,ochotona_princeps,oreochromis_niloticus,ornithorhynchus_anatinus,oryctolagus_cuniculus,oryzias_latipes,otolemur_garnettii,ovis_aries,pan_troglodytes,papio_anubis,pelodiscus_sinensis,petromyzon_marinus,poecilia_formosa,pongo_abelii,procavia_capensis,pteropus_vampyrus,rattus_norvegicus,sarcophilus_harrisii,sorex_araneus,sus_scrofa,taeniopygia_guttata,tarsius_syrichta,tetraodon_nigroviridis,tursiops_truncatus,vicugna_pacos,xenopus_tropicalis,saccharomyces_cerevisiae"

genomes = str_genomes.strip().split(',')

# ensembl rest server 
server = "https://rest.ensembl.org"


# simple class used in order to 
# wrap conveniently 
class HomoStruct:
   def __init__(self, **entries):
      self.__dict__.update(entries)
   def to_fasta(self):
      s = str(self.align_seq).translate(None, '-')
      return ''.join( s.strip().split('\n') )
        
# step 1 
def gather_gene_ids(h):
   names = set()
   for head in h: 
      ssh = head.strip().split('_')
      names.add( ssh[3].split('.')[0] )
   return list(sorted(names))
   
# step 2 
def request_ens_homo(ids):
   db = {}   
   req_count = 0 # count nb of requests 
   last_req = time.time() # store time before a 15 cycle

   for i in ids: 
      # check if we need to rate limit ourselves
      if req_count >= 5:
         delta = time.time() - last_req
         if delta < 1:
            time.sleep(1 - delta)
            last_req = time.time()
            req_count = 0
            
      ext = "/homology/id/%s?type=orthologues" % i
      print >>sys.stderr, '* request : %s ( %s;content-type=application/json )'  % (i,server+ext)
      r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
      req_count += 1
      
      if not r.ok:
         r.raise_for_status()
         sys.exit()
      db[i] = r.json()
   return db 


if __name__ == '__main__':
   heads = [x for x in sys.stdin.readlines() if x.startswith('>')]
   ids = gather_gene_ids(heads)
   dbj = request_ens_homo(ids)
   fam_members= set() # list of family members 
     
   # 3 - dump orthologous gene in a fasta file 
   for k in dbj.keys():
      
      # test if 
      try:
         root = dbj[k]['data'][0][u'homologies']
      except IndexError:
         print >>sys.stderr, '#no data for : %s '  % k
         continue
         
      # keeps only one 2 one orthologous proteins 
      froot = [x for x in root if x['type'] == u'ortholog_one2one']

      # alloc object from json anwsers
      s = [ HomoStruct(**x['source']) for x in froot ]
      h = [ HomoStruct(**x['target']) for x in froot ]

      # dump source seqs (in the seeds)
      for ss in s:
         if (ss.id not in fam_members) and (ss.species in genomes):
            fasta_head = ">%s_%s|<SEED>|%s" % (ss.id, ss.protein_id, ss.species)
            print fasta_head + '\n'+ ss.to_fasta()
            fam_members.add(ss.id)

      # then dump orthologous genes
      for hh in h:
         #print dir(hh)
         if (hh.id not in fam_members) and (hh.species in genomes):
            fasta_head = ">%s_%s|%s" % (hh.id, hh.protein_id, hh.species)
            print fasta_head + '\n'+ hh.to_fasta()
            fam_members.add(hh.id)
