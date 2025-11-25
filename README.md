# BME-Term-Project

**Code descriptions (libraries, output, and running)**

**Dataset_extraction.py**

FOXO3 sequences are extracted from Uniprot, Ensembl, AND NCBI, in that order.
If a sequence isn't found in Uniprot it's then searched in Ensembl, and if not found there
it searches NCBI last. All NCBI sequences are inputted in BLAST to check for orthology using human FOXO3 sequence.

The sequences are then split into 3 files, the full file has the entire sequences, the DBD file
has only the DNA binding domain, using the coordinates given by the human FOXO3 sequence in uniprot (157-251).
Finally the IDR file has all of the other (likely disordered regions).

*Libraries:*

- requests
- time
- urllib.parse
- subprocess
- tempfil
- Bio (biopython) 
- sys

*Test run (only 5 species):*

```
python3 Dataset_extraction.py test
```

*Full run (all species):*

```
python3 Dataset_extraction.py full
```

*Output ({name} = full or test):*

- {name}_foxo3_seq.fasta
- {name}_foxo3_DBD.fasta
- {name}_foxo3_IDR.fasta

**MSA.py**

The MSA file takes the sequence files from Data_extraction.py, aligns them using MAFFT and then trims the alignments
using TrimAl (in order to reduce noise). The test run will output every file (for manually checking quality), so there will be the alignment files before trimming, and the full final sequencs after alignment and trimming.
The full run will only output final sequences.

*Libraries*

- os
- shutil
- subprocess

*Test run (only 5 species):*

```
python3 MSA.py test
```

*Full run (all species):*

```
python3 MSA.py full
```

*TEST RUN Output files:*
- test_full_seq_aligned.fasta
- test_foxo3_DBD.fasta
- test_foxo3_IDR.fasta

- test_full_final.fasta
- test_DBD_final.fasta
- test_IDR_final.fasta

*FULL RUN Output files:*

- full_full_final.fasta
- full_DBD_final.fasta
- full_IDR_final.fasta
  
# basic download information #
Biopython:

```pip install biopython```

Requests:

```pip install requests```

MAFFT:

```sudo apt-get install mafft```

Trimal:
```
git clone https://github.com/inab/trimal.git

cd trimal
cd Source

make
```


**other sorting stuff this part is not organized at all**
```sudo apt install -y python3 python3-pip

sudo apt install python3.12-venv

python3 -m venv venv

source venv/bin/activate

sudo apt update
sudo apt install -y build-essential git

git clone https://github.com/inab/trimal.git
cd trimal
cd Source

# compile
make

# move binary to PATH
sudo cp trimal /usr/local/bin/
sudo chmod +x /usr/local/bin/trimal

# test
trimal -h

pip install biopython
```
