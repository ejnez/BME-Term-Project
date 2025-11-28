import requests
import time
import urllib.parse
import subprocess
import tempfile
import csv

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

def build_fasta_dataset(species_list, lifespan_dict, common_name_dict, outfile, blast_db="human_db"):
    collected = []
    missing_species = []
    weeded_out = []

    # Ensure human first
    if "Homo sapiens" in species_list:
        species_list = ["Homo sapiens"] + [s for s in species_list if s != "Homo sapiens"]

    for sp in species_list:
        seq_data = fetch_foxo3_sequence(sp)
        if not seq_data:
            missing_species.append(sp)
            continue

        # Split entries if multiple sequences returned
        entries = [s for s in seq_data.split("\n>") if s.strip()]
        entries = [(">" + e if not e.startswith(">") else e) for e in entries]

        # Sequence prioritization: sp| > tr| > first entry
        selected_seq = next((e for e in entries if e.startswith(">sp|")), None)
        if not selected_seq:
            selected_seq = next((e for e in entries if e.startswith(">tr|")), None)
        if not selected_seq:
            selected_seq = entries[0]

        # If sequence looks like NCBI (gi| or ref|), do BLAST QC
        lines = selected_seq.split("\n")
        seq = "".join(lines[1:])
        fasta = f">tmp\n{seq}\n"
        if selected_seq.startswith(">gi") or selected_seq.startswith(">ref"):
            if not is_valid_ncbi_sequence(fasta, blast_db=blast_db):
                weeded_out.append(sp)
                continue

        # Clean header: Genus_Species_CommonName|AvgLifespan
        try:
            genus, species_name = sp.split()  # assumes "Genus species"
        except ValueError:
            genus, species_name = sp, ""  # fallback if single word

        common = common_name_dict.get(sp, "Unknown").replace(" ", "_")
        avg_life = lifespan_dict.get(sp, 0)
        clean_header = f">{genus}_{species_name}_{common}|{avg_life}"

        body = "\n".join(lines[1:])
        collected.append(clean_header + "\n" + body)

        time.sleep(0.25)  # be polite to APIs

    # Write to FASTA
    with open(outfile, "w") as f:
        f.write("\n\n".join(collected))

    print(f"\nSaved {len(collected)} sequences to {outfile}")
    if missing_species:
        print(f"Not found: {missing_species}")
    if weeded_out:
        print(f"Weeded out (failed BLAST QC): {weeded_out}")

def load_species_table(filename):
    """
    Reads species file with columns:
    Species, Order, AvgLifespan, MaxLifespan, Common name
    Handles tab- or comma-separated files.
    """
    species_list = []
    lifespan_dict = {}
    common_name_dict = {}

    with open(filename, "r", newline="", encoding="utf-8-sig") as f:
        first_line = f.readline()
        delimiter = "\t" if "\t" in first_line else ","
        f.seek(0)

        reader = csv.DictReader(f, delimiter=delimiter)
        # normalize headers: strip spaces
        reader.fieldnames = [name.strip() for name in reader.fieldnames]

        for row in reader:
            # strip spaces from keys and values
            row = {k.strip(): v.strip() for k, v in row.items()}

            species = row.get("Species", "")
            if not species:
                continue  # skip empty species rows

            common_name = row.get("Common name", "Unknown")
            try:
                avg_life = float(row.get("AvgLifespan", 0))
            except ValueError:
                avg_life = 0

            species_list.append(species)
            lifespan_dict[species] = avg_life
            common_name_dict[species] = common_name

    return species_list, lifespan_dict, common_name_dict


def extract_dbd_for_all_sequences(fasta_file, pfam_db, pfam_id="PF00250", domain_name="Forkhead", out_file="dbd_sequences.fasta", domE=1e-3):
    """
    Extract domains from FASTA using hmmscan, filtering by Pfam accession or name.
    """
    import subprocess, tempfile, os
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    def clean_id(s):
        return s.split("|")[0].split()[0]

    # Load sequences
    records = {clean_id(r.id): r for r in SeqIO.parse(fasta_file, "fasta")}
    hits_found = set()
    dbd_records = []

    # Temp FASTA
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as tmp_fasta:
        SeqIO.write(records.values(), tmp_fasta, "fasta")
        tmp_name = tmp_fasta.name

    # Temp domtblout
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as tmp_out:
        out_name = tmp_out.name

    # Run HMMer
    subprocess.run([
        "hmmscan", "--domtblout", out_name, "--noali",
        "--domE", str(domE), "--incdomE", str(domE),
        pfam_db, tmp_name
    ], capture_output=True)

    # Parse domtblout
    with open(out_name) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 22:
                continue
            # Filter by name OR accession (starts with pfam_id)
            if not (parts[0] == domain_name or parts[1].startswith(pfam_id)):
                continue

            seq_id = clean_id(parts[3])
            start, end = int(parts[17]), int(parts[18])
            seq = str(records[seq_id].seq)
            dbd_seq = seq[start-1:end]
            dbd_records.append(SeqRecord(Seq(dbd_seq), id=f"{records[seq_id].id}_DBD", description=""))
            hits_found.add(seq_id)

    # Warn missing
    for r in records.values():
        if clean_id(r.id) not in hits_found:
            print(f"WARNING: No PFAM domain {pfam_id} ({domain_name}) found for {r.id}")

    # Write output
    SeqIO.write(dbd_records, out_file, "fasta")
    print(f"Extracted DBDs written to: {out_file}")

    # Clean up
    os.remove(tmp_name)
    os.remove(out_name)
    
def main():
    if len(sys.argv) < 2:
        print("Usage: python Dataset_extraction.py <PrefixName>")
        sys.exit(1)

    prefix = sys.argv[1]

    species_file = f"{prefix}_Species.csv"
    fasta_outfile = f"{prefix}_foxo3_seq.fasta"
    
    # Load species table
    species_list, lifespan_dict, common_name_dict = load_species_table(species_file)

    # Build FASTA dataset
    build_fasta_dataset(
        species_list,
        lifespan_dict,
        common_name_dict,
        fasta_outfile,
        blast_db="human_db"
    )
    PFAM_DB = "Pfam-A.hmm"  # path to your Pfam HMM database
    print("Finished building FASTA dataset:", fasta_outfile)
    extract_dbd_for_all_sequences(fasta_outfile, PFAM_DB, pfam_id="PF00250", out_file=dbd_outfile)
    dbd_outfile = f"{prefix}_foxo3_DBD.fasta"
   


    

if __name__ == "__main__":
    main()
    