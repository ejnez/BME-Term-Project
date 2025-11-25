import requests
import time
import urllib.parse
import subprocess
import tempfile

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys

def get_tax_id(species_name):
    """
    Gets Taxonomy ID from UNIPROT
    """
    url = "https://rest.uniprot.org/taxonomy/search"
    params = {"query": f'"{species_name}"', "fields": "taxon_id,scientific_name,common_name"}
    resp = requests.get(url, params=params)
    data = resp.json()
    if "results" in data and len(data["results"]) > 0:
        return data["results"][0]["taxonId"]
    # fallback fuzzy search
    params = {"query": species_name}
    resp = requests.get(url, params=params)
    data = resp.json()
    if "results" in data and len(data["results"]) > 0:
        return data["results"][0]["taxonId"]
    return None

def fetch_uniprot(species_name):
    """
    Uniprot search, first for quality results.
    """
    tax_id = get_tax_id(species_name)
    if tax_id is None:
        return None

    query = "(" \
            "gene:FOXO3* OR gene:foxo3* OR gene:Foxo3* OR " \
            "protein_name:\"Forkhead box protein O3*\"" \
            f") AND organism_id:{tax_id}"
    encoded_query = urllib.parse.quote(query)
    url = f"https://rest.uniprot.org/uniprotkb/stream?format=fasta&compressed=false&query={encoded_query}"
    headers = {"Accept": "text/x-fasta"}
    resp = requests.get(url, headers=headers)
    if resp.status_code == 200 and resp.text.startswith(">"):
        return resp.text
    return None

def fetch_ensembl(species_name):
    """
    Ensemble search, more quality than NCBI but less that Uniprot
    """
    species = species_name.lower().replace(" ", "_")
    gene = "FOXO3"
    lookup_url = f"https://rest.ensembl.org/lookup/symbol/{species}/{gene}?expand=1"
    headers = {"Content-Type": "application/json"}
    resp = requests.get(lookup_url, headers=headers)
    if resp.status_code != 200:
        return None
    data = resp.json()
    gene_id = data.get("id")
    if not gene_id:
        return None
    seq_url = f"https://rest.ensembl.org/sequence/id/{gene_id}?type=protein"
    seq_resp = requests.get(seq_url, headers={"Content-Type": "text/x-fasta"})
    if seq_resp.status_code == 200 and seq_resp.text.startswith(">"):
        header, *body = seq_resp.text.split("\n")
        clean_header = f">{species_name.replace(' ', '_')}"
        return clean_header + "\n" + "\n".join(body)
    return None

def fetch_ncbi(species_name):
    """
    Final last case scenario NCBI search
    """
    gene = "FOXO3"
    term = f"{gene}[Gene] AND {species_name}[Organism]"
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {"db": "protein", "term": term, "retmax": 5, "retmode": "json"}
    resp = requests.get(url, params=params)
    if resp.status_code != 200:
        return None
    data = resp.json()
    ids = data.get("esearchresult", {}).get("idlist", [])
    if not ids:
        return None
    # Fetch first protein sequence
    fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    fetch_params = {"db": "protein", "id": ids[0], "rettype": "fasta", "retmode": "text"}
    seq_resp = requests.get(fetch_url, params=fetch_params)
    if seq_resp.status_code == 200 and seq_resp.text.startswith(">"):
        header, *body = seq_resp.text.split("\n")
        clean_header = f">{species_name.replace(' ', '_')}"
        return clean_header + "\n" + "\n".join(body)
    return None

def fetch_foxo3_sequence(species_name):
    """
    Fetch sequence function:

    Searches uniprot first, if it failes Ensemble, and then NCBI
    """
    print(f"\n...Searching {species_name} in UniProt...")
    seq = fetch_uniprot(species_name)
    if seq:
        print("DONE:Found in UniProt")
        return seq

    print(f"...Searching {species_name} in Ensembl...")
    seq = fetch_ensembl(species_name)
    if seq:
        print("DONE:Found in Ensembl")
        return seq

    print(f"...Searching {species_name} in NCBI...")
    seq = fetch_ncbi(species_name)
    if seq:
        print("DONE:Found in NCBI")
        return seq

    print("FAILED:not found in any database")
    return None


def is_valid_ncbi_sequence(seq_fasta, blast_db="human_db"):
    """
    BLAST candidate NCBI sequence against human FOXO3.
    Returns True if FOXO3 is the top hit with reasonable coverage.
    """
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as tmp:
        tmp.write(seq_fasta)
        tmp_name = tmp.name

    cmd = [
        "blastp",
        "-query", tmp_name,
        "-db", blast_db,
        "-outfmt", "6 qlen slen qstart qend sstart send pident bitscore sseqid",
        "-max_target_seqs", "1",
        "-evalue", "1e-5"
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    out = result.stdout.strip()

    if not out:
        return False

    fields = out.split("\t")
    qlen, slen = int(fields[0]), int(fields[1])
    qstart, qend = int(fields[2]), int(fields[3])
    sseqid = fields[8]

    coverage = (qend - qstart + 1) / qlen

    if "FOXO3" not in sseqid.upper():
        return False
    if coverage < 0.70:
        return False

    return True

def split_domains(infile, name):
    """
    Split sequences into DBD (157–251) and IDR regions
    """
    DBD_START = 157
    DBD_END = 251

    dbd_records = []
    idr_records = []

    for record in SeqIO.parse(infile, "fasta"):
        seq = str(record.seq)

        # Skip sequences too short for DBD
        if len(seq) < DBD_END:
            continue

        dbd_seq = seq[DBD_START-1:DBD_END]
        idr_seq = seq[:DBD_START-1] + seq[DBD_END:]

        # DBD record
        r_dbd = SeqRecord(
            Seq(dbd_seq),
            id=record.id + "_DBD",
            description=""
        )
        dbd_records.append(r_dbd)

        # IDR record
        r_idr = SeqRecord(
            Seq(idr_seq),
            id=record.id + "_IDR",
            description=""
        )
        idr_records.append(r_idr)

    dbd_out = f"{name}_foxo3_DBD.fasta"
    idr_out = f"{name}_foxo3_IDR.fasta"

    SeqIO.write(dbd_records, dbd_out, "fasta")
    SeqIO.write(idr_records, idr_out, "fasta")

    print(f"DBD sequences written to: {dbd_out}")
    print(f"IDR sequences written to: {idr_out}")

def build_fasta_dataset(species_list, outfile, blast_db="human_db"):
    """
    Builds FOXO3 dataset with:
    - human first
    - priority: UniProt sp| > tr| > NCBI
    - BLAST QC for NCBI sequences
    """
    collected = []
    missing_species = []
    weeded_out = []

    # ensure human is first
    if "Homo sapiens" in species_list:
        species_list = ["Homo sapiens"] + [s for s in species_list if s != "Homo sapiens"]

    for sp in species_list:
        seq_data = fetch_foxo3_sequence(sp)
        if not seq_data:
            missing_species.append(sp)
            continue

        entries = [s for s in seq_data.split("\n>") if s.strip()]
        entries = [(">"+e if not e.startswith(">") else e) for e in entries]

        # 1) Prefer sp|
        selected_seq = next((e for e in entries if e.startswith(">sp|")), None)

        # 2) else tr|
        if not selected_seq:
            selected_seq = next((e for e in entries if e.startswith(">tr|")), None)

        # 3) else: Ensembl or NCBI
        if not selected_seq:
            candidate = entries[0]

            # If this came from Ensembl fetch, the header is exactly ">Genus_species"
            if candidate.split("\n")[0].startswith(f">{sp.replace(' ', '_')}"):
                # GOOD — Ensembl sequence, keep it
                selected_seq = candidate
            else:
                # Probably NCBI → BLAST QC required
                lines = candidate.split("\n")
                seq = "".join(lines[1:])
                fasta = f">tmp\n{seq}\n"

                if is_valid_ncbi_sequence(fasta, blast_db=blast_db):
                    selected_seq = candidate
                else:
                    weeded_out.append(sp)
                    continue

        # Clean header to >Genus_species
        lines = selected_seq.split("\n")
        body = lines[1:]
        clean_header = f">{sp.replace(' ', '_')}"
        collected.append(clean_header + "\n" + "\n".join(body))

        time.sleep(0.25)

    with open(outfile, "w") as f:
        f.write("\n\n".join(collected))

    print(f"\nSaved {len(collected)} sequences to {outfile}")
    if missing_species:
        print(f"Not found: {missing_species}")
    if weeded_out:
        print(f"Weeded out (failed BLAST QC): {weeded_out}")

def main():
    if len(sys.argv) < 2:
        print("Usage: python Dataset_extraction.py [test|full|myname]")
        sys.exit(1)

    prefix = sys.argv[1]

    if prefix[0] == 't':
        name = 'test'

        species_list = [ "Homo sapiens",
                 "Mus musculus",
                 "Gallus gallus",
                 "Danio rerio",
                 "Anolis carolinensis" ]
        
    else:
        name = 'full'

        # REAL DATASET
        species_list = [
            "Balaena mysticetus",
            "Balaenoptera physalus",
            "Orcinus orca",
            "Myotis brandtii",
            "Myotis myotis",
            "Heterocephalus glaber",
            "Homo sapiens",
            "Pongo abelii",
            "Pan troglodytes",
            "Elephas maximus",
            "Equus caballus",
            "Ursus arctos",
            "Cacatua galerita",
            "Psittacus erithacus",
            "Strigops habroptilus",
            "Aquila chrysaetos",
            "Diomedea exulans",
            "Corvus corax",
            "Somniosus microcephalus",
            "Cyprinus rubrofuscus",
            "Hoplostethus atlanticus",
            "Chelonoidis niger",
            "Aldabrachelys gigantea",
            "Sphenodon punctatus",
            "Balaenoptera acutorostrata",
            "Tursiops truncatus",
            "Myotis lucifugus",
            "Mus musculus",
            "Rattus norvegicus",
            "Macaca mulatta",
            "Sus scrofa",
            "Canis lupus familiaris",
            "Felis catus",
            "Gallus gallus",
            "Taeniopygia guttata",
            "Fulmarus glacialis",
            "Corvus brachyrhynchos",
            "Scyliorhinus canicula",
            "Danio rerio",
            "Oryzias latipes",
            "Chrysemys picta",
            "Anolis carolinensis"
        ]

    full_outfile = f"{name}_foxo3_seq.fasta"

    build_fasta_dataset(species_list, full_outfile)
    split_domains(full_outfile, name)



if __name__ == "__main__":
    main()