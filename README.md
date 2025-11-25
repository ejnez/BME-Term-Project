# BME-Term-Project

**Libraries used by each code:**

Dataset_extraction.py

- requests
- time
- urllib.parse
- subprocess
- tempfil
- Bio (biopython) 
- sys

MSA.py

- os
- shutil
- subprocess

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


**other sorting stuff this part is not organized at all **
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
