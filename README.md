# BME-Term-Project

**Code descriptions (libraries, output, and running)**

**Dataset_extraction.py**

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
