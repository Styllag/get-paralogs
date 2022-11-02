#βιβλιοθήκες
from urllib.request import urlopen
import re
from bs4 import BeautifulSoup
import csv
import argparse


my_parser = argparse.ArgumentParser(prog='findparalogs', description='Find paralogs per species for a given KO')
my_parser.add_argument('KO',metavar='ko',type=str,help='the KO to use')                                    
my_parser.add_argument('-a', action='store_const', const="Animals", help="ψάχνει  για παράλογα ζώων")
my_parser.add_argument('-p', action='store_const', const="Plants", help="ψάχνει  για παράλογα φυτών")
my_parser.add_argument('-f', action='store_const', const="Fungi", help="ψάχνει  για παράλογα μυκήτων")
my_parser.add_argument('-r', action='store_const', const="Protists", help="ψάχνει  για παράλογα πρωτίστων")
my_parser.add_argument('-k','--komore', action='store_true', help="λαμβάνει ως παράλογα όλα τα γονίδια ενός είδους που έχουν το ίδιο ΚΟ ανεξαρτήτως identity score")
args = my_parser.parse_args()

taxalist = []
if args.p:
    taxalist.append("Plants")
if args.f:
    taxalist.append("Fungi")
if args.r:
    taxalist.append("Protists")
if args.a:
    taxalist.append("Animals")

    




#αρχική λίστα για ένα κο

gene_line=re.compile(r'[A-Z ]+ (?P<sp>[A-Z]{3,4}):(?P<rest>(?: \w+(?:\(\w+\))?)+)')


#η συνάρτηση αυτή επιστρέφει μια λίστα με όλα τα γονίδια που αντιστοιχούν σε ένα KO στην kegg
def get_gene_list(ko):
    list_of_genes = []
    genes = False
    with urlopen(f"http://rest.kegg.jp/get/{ko}") as response:
        for line in response:
            line = line.decode()
            if line.startswith("GENES"):
                genes = True
            if line.startswith("REFERENCE"):
                genes = False
            if genes: 
                list_of_genes.append(line)
    return list_of_genes            


#η συνάρτηση αυτή βγάζει τις παρενθέσεις
def clean_from_parentheses(rest):
    out_par = 1
    Ret = ""
    for char in rest:
        if char == "(":
            out_par = 0
        elif char == ")":
            out_par = 1
        else:
            if out_par:
                Ret += char
    return Ret

#η συνάρτηση αυτή διαχωρίζει το όνομα του γονδίου από το όνομα του είδους και επιστρέφει και τα δύο
def organize_genes(line):
           if gene_line.match(line):
                    species_code = gene_line.match(line).group(1)
                    rest_gene = gene_line.match(line).group("rest")
                    gene_id = clean_from_parentheses(rest_gene).strip().split(" ")
                    return [species_code,gene_id]        



InitData = {}
a = get_gene_list(args.KO)   
for i in a:
    InitData[organize_genes(i)[0]] = organize_genes(i)[1]


#η συνάρτηση αυτή είναι η πιο σημαντική. Δέχεται ως όρισμα το όνομα ενός γονιδίου και βρίσκει τα παράλογα αυτού του γονιδίου.
#για να βρει τα παράλογα χρησιμοποιεί regular expressions
def find_paralogs_all(gene):
    with urlopen(f'https://www.kegg.jp/ssdb-bin/ssdb_paralog?org_gene={gene}') as response:
        html = response.read().decode()
        soup = BeautifulSoup(html, features = "lxml")
        table_tag = soup.find_all("pre")
        if table_tag:
            table = str(table_tag[0])
            list_with_paralogs = re.findall("<input checked=.+/><a href=.+>([a-z]{3,4}:\w+)</a> .+ (?:<a href=.+>K\d{5}</a>)? +(\d+) +(\d+) \(.+\) +(\d+) +([01]\.\d{3}) +(\d+) ",table)
            return list_with_paralogs
        else:
            p = soup.p
            return p.string


#η συνάρτηση αυτή φιλτράρει τα παράλογα με βάση το identity score
def filter_paralogs_with_id07(parlist):
    if type(parlist) == list:
        filtered_list = [x for x in parlist if float(x[4]) >= 0.7]
        return filtered_list
    else :
        print(parlist)    



#χρησιμοποιείται ένα αρχείο csv που κατηγοριοποιεί το είδος με βάση την ευρύτερη ταξινομική ομάδα (ζώα, φυτά κλπ)
taxon = {}
with open('species_code_per_taxon.csv',"r") as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for row in csv_reader:
        taxon[row[0]] = row[1]



#κύρια διαδικασία του προγράμματος ενοποιεί όλες τις συναρτήσεις και γράφει τα αποτελέσματα του γενικού φιλτραρίσματος σε ένα αρχείο
filteredData = set()
for species in InitData:
        if species.lower() in taxon.keys():
            if taxon[species.lower()] in taxalist:
                print(f"γίνεται αναζήτηση για τα παράλογα του είδους {species}")
                for gene in InitData[species]:
                    gene_name = f'{species.lower()}:{gene}'
                    if len(InitData[species]) > 1 :
                        if args.komore:
                            print(f"το είδος {species} έχει περισσότερα από ένα γονίδια που ανήκουν στο συγκεκριμένο KO ({args.KO}).")
                            filteredData.add(gene_name)
                    parlist_ssdb = filter_paralogs_with_id07(find_paralogs_all(gene_name))
                    if parlist_ssdb :
                        print(f"το γονίδιο {gene_name} έχει στην ssdb τουλάχιστον ένα παράλογο με identity > 0.7")
                        for par in parlist_ssdb :
                            filteredData.add(gene_name)
                            filteredData.add(par[0])

total_number_of_paralogs = len(filteredData)
print(f"Βρέθηκαν συνολικά {total_number_of_paralogs} παράλογα για το {args.KO}")

with open(f"{args.KO}.txt","w") as f:
    for i in sorted(filteredData):
        f.write(i + '\n')        



