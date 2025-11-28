# BME-Term-Project

**Code descriptions (libraries, output, and running)**

# Dataset_extraction.py

There are 5 csv files named {order}_Species.csv, each file contains a list of species, their order, their average lifespan, their maximum lifespan, and their common name.

From this species name, average lifespan, and common name are used as fasta titles (ex: {Genus}_{species}_{Common}_{name}|AvgLifespan)

FOXO3 sequences are extracted from Uniprot, Ensembl, AND NCBI, in that order.
If a sequence isn't found in Uniprot it's then searched in Ensembl, and if not found there
it searches NCBI last. All NCBI sequences are inputted in BLAST to check for orthology using human FOXO3 sequence.

The full sequence is then used by HMMer alongside the Pfam database in order to find the "Forkhead" domain with the Pfam ID "PF00250". This is a more accurate DNA Binding Domain than just using one species coordinates. Then an {order}_foxo3_DBD.fasta is created with every species DBD sequence.

*Orders:*
- Artiodactyla
- Carnivora
- Chiroptera
- Primate
- Rodentia
  
*Libraries:*

- requests
- time
- urllib.parse
- subprocess
- tempfile
- Bio (biopython) 
- sys


*run 1 order at a time:*

```
python3 Data_extraction.py {order}
```

*Output ({name} = full or test):*

DISPLAYING CHIROPTERA EXAMPLE ONLY

- Chiroptera_foxo3_seq.fasta
  
```
>Eptesicus_fuscus_Big_brown_bat|19.0
MRVQNEGTGKSSWWIINPDGGKSGKAPRRRAVSMDNSNKYTKSRSRAAKKKAALQTAPESADDSPSQLSK
WPGSPTSRSSDELDAWTDFRSRTNSNASTVSGRLSPILASTELDDVQDDDAPLSPMLYSSSASLSPSVSK
PCTVELPRLTDMAGTMNLNDGLADNLMDDLLDNIPLPPSQPSPTAGLMQRSSSFPYTTKGSGLGSPTSSF
NSTVFGPSSLNSLRQSPMQTIQENKPATFSSMSHYGTQTLQDLLTSDSLSHSDVMMTQSDPLMSQASTAV
SAQNSRRNVMLRNDPMMSFAAQPNQGSLVNQNLLHHQHQTQGALSGSRALSNSVSNMGLSDSSSLGSAKH
QQQSLVSQSMQTLSDSLSGSSLYSTSANLSVMGHEKFPSDLDLDMFNGSLECDMESIIRTDLMDADGLDF
NFDSLISTQNVVGLNVGNFTGAKQASSQSWVPG



>Myotis_lucifugus_Little_brown_bat|7.0
MAEAPASPAPLSPLEVELDPEFEPQSRPRSCTWPLQRPELQGSPAKPSGEAAADSMIPEE
EDDEDDEDGGGRAGSAMAIGGGGSSTLGAGMLLEDSARLLAPGGQDPGSGPAPAAGALSG
GTQTPLQPQQPLPPPQPGAAGGSGQPRKCSSRRNAWGNLSYADLITRAIQSSPDQRLTLS
QIYEWMVRCVPYFNDKGDSNSSAGWKNSIRHNLSLHSRFMRVQNEGTGKSSWWIINPDGG
KSGKAPRRRAVSMDNSNKYTKSRSRAAKKKAALQTAPESADDSPSQLSKWPGSPTSRSSD
ELDAWTDFRSRTNSNASTVSGRLSPILAGTELDDVQDDDAPLSPMLYSSSASLSPSVSKP
CTVELPRLTDMAGTMNLNDGLADNLMDDLLDNIPLPPSQPSPTAGLMQRSSSFPYTTKGS
GLGSPTSSFNSTVFGPSSLNSLRQSPMQTIQENKPATFSSMPHYGTQTLQDLLTSDSLSH
SDVMMTQSDPLMSQASTAVSAQNSRRNVMLRNDPMMSFAAQPNQGSLVNQNLLHHQHQTQ
GALSGSRALSNSVSNMGLSDSSSLGSAKHQQQSLVSQSMQTLSDSLSGSSLYSTSANLSV
MGHEKFPSDLDLDMFNGSLECDMESIIRSDLMDADGLDFNFDSLISTQNVVGLNVGNFTG
AKQASSQSWVPG


>Pteropus_alecto_Black_flying_fox|4.5
MAEAPASPAPLSPLEVELDPEFEPQSRPRSCTWPLQRPELQGSPAKPSGEAAADSMIPEE
EDDEDDEDSSSRAGSAMAIGGGGCGTLGAGLLLEDSARLLAPGGQDPGSGPAPAAGAPSG
ATQTPLQPQQPLPPPQPGAAGGSGQPRKCSSRRNAWGNLSYADLITRAIESSPDKRLTLS
QIYEWMVRCVPYFKDKGDSNSSAGWKVRTHLGRAAGPAGPPAQYETRLRFVGALLQAPGR
GVGREPLL

>Rousettus_aegyptiacus_Egyptian_rousette|9.0
MAEALASPAPLSPLEVELDPEFEPQSRPRSCTWPLQRPELQGSPAKPSGEAAADSMIPEEEDDEDDEDSS
SRAGSAMAIGGGGCGTLGAGLLLEDSARLLAPGGQDPGSGPAPAAGAPSGATQTPLQPQQPLPPPQPGAA
GGSGQPRKCSSRRNAWGNLSYADLITRAIESSPDKRLTLSQIYEWMVRCVPYFKDKGDSNSSAGWKNSIR
HNLSLHSRFMRVQNEGTGKSSWWIINPDGGKSGKAPRRRAVSMDNSNKYTKSRNRAAKKKAALQAAPESA
DDSPSQLSKWPGSPTSRSSDELDAWTDFRSRTNSNASTVSGRLSPILASTELDDVQDDDAPLSPMLYSSS
ASLSPSVSKPCTVELPRLTDMTGTMNLNDGLADNLMDDLLDNITLPSSQPSPTGGLMQRSSSFPYTTKGS
GLGSPASSFNSTVFGPSSLNSLRQSPMQTIQENKPATFSSMSHYSTQTLQDLLASDSLSHSDVMMTQSDP
LMSQASTAVSAQNSRRNVMLRNDPMMSFAAQPNQGSLVNQNLLHHQHQTQGALGGSRALSNSVSNMGLSD
SSTLGSAKHQQQSPVSQSMQTLSDSLSGSSLYSTSANLSVMGHEKFPSDLDLDMFNGSLECDMESIIRSE
LMDADGLDFNFDSLISTQNVVGLNVGNFTGAKQASSQSWVPG



>Myotis_brandtii_Brandts_bat|40.0
MVRCVPYFNDKGDSNSSAGWKNSIRHNLSLHSRFMRVQNEGTGKSSWWIINPDGGKSGKAPRRRAVSMDN
SNKYTKSRSRAAKKKAALQTAPESADDSPSQLSKWPGSPTSRSSDELDAWTDFRSRTNSNASTVSGRLSP
ILAGTELDDVQDDDAPLSPMLYSSSASLSPSVSKPCTVELPRLTDMAGTMNLNDGLADNLMDDLLDNIPL
PPSQPSPAAGLMQRSSSFPYTTKGSGLGSPTSSFNSTVFGPSSLNSLRQSPMQTIQENKPATFSSMPHYG
TQTLQDLLTSDSLSHSDVMMTQSDPLMSQASTAVSAQNSRRNVMLRNDPMMSFAAQPNQGSLVNQNLLHH
QHQTQGALSGSRALSNSVSNMGLSDSSSLGSAKHQQQSLVSQSMQTLSDSLSGSSLYSTSANLSVMGHEK
FPSDLDLDMFNGSLECDMESIIRSDLMDADGLDFNFDSLISTQNVVGLNVGNFTGAKQASSQSWVPG



>Myotis_myotis_Mouse-eared_bat|15.5
MAEAPASPAPLSPLEVELDPEFEPQSRPRSCTWPLQRPELQGSPAKPSGEAAADSMIPEEEDDEDDEDGG
GRAGSAMAIGGGGSSTLGAGMLLEDSARLLAPGGQDPGSGPAPAAGAPSGGTQTPLQPQQPLPPPQPGAA
GGSGQPRKCSSRRNAWGNLSYADLITRAIQSSPDQRLTLSQIYEWMVRCVPYFNDKGDSNSSAGWKNSIR
HNLSLHSRFMRVQNEGTGKSSWWIINPDGGKSGKAPRRRAVSMDNSNKYTKSRSRAAKKKAALQTAPESA
DDSPSQLSKWPGSPTSRSSDELDAWTDFRSRTNSNASTVSGRLSPILAGTELDDVQDDDAPLSPMLYSSS
ASLSPSVSKPCTVELPRLTDMAGTMNLNDGLADNLMDDLLDNIPLPPSQPSPTAGLMQRSSSFPYTTKGS
GLGSPTSSFNSTVFGPSSLNSLRQSPMQTIQENKPATFSSMPHYGTQTLQDLLTSDSLSHSDVMMTQSDP
LMSQASTAVSAQNSRRNVMLRNDPMMSFAAQPNQGSLVNQNLLHHQHQTQGALSGSRALSNSVSNMGLSD
SSSLGSAKHQQQSLVSQSMQTLSDSLSGSSLYSTSANLSVMGHEKFPSDLDLDMFNGSLECDMESIIRSD
LMDADGLDFNFDSLISTQNVVGLNVGNFTGAKQASSQSWVPG



>Rhinolophus_ferrumequinum_Greater_horseshoe_bat|30.0
MAEARASPAPLSPLEVELDPEFEPQSRPRSCTWPLQRPELQGSPAKPSGEAAADSMIPEE
EDDEDDEDGGGRASSAMAIGGGGSGTLSAGLLLEDSARLLAPGGQDPGSGPAPAAGALSG
ATQTPLQPQQPLPPPQPGAAGGSGQPRKCSSRRNAWGNLSYADLITRAIESSPDKRLTLS
QIYEWMVRCVPYFKDKGDSNSSAGWKNSIRHNLSLHSRFMRVQNEGTGKSSWWIINPDGG
KSGKAPRRRAVSMDNSNKYTKSRSRAAKKKAALQAAPEAADDSPSQLSKWPGSPTSRSSD
ELDAWTDFRSRTNSNASTVSGRLSPILASTELDDVQDDDAPLSPMLYSSASSLSPSVSKP
CTVELPRLTDMAGTMNLNDGLADNLMDDLLDNIALPSSQPSPTGGLMPRSSSFPYTPKGS
GLGSPTSSFTSTVFGPSSLTSLRQSPMQTIQENKPATFSSMSHYGTQTLQDLLTSDSLSH
SDVMMTQSDPLMSQASTAVSAQSSRRNVMLRSDPMMSFAAQPNQGSLVNQNLLHHQHQSQ
GALGGSRALSNSVSSMGLSDASSLGAAKHQQQSPGSQSMQTLSDTLSGSSLYSTSANLSV
MGHDKFPSDLDLDIFNGSLECDMESIIRSELMDADGLDFNFDSLISTQNVVGLNVGNFTG
AKQASSQSWVPG
```
  
- Chiroptera_foxo3_DBD.fasta

```
>Myotis_lucifugus_Little_brown_bat|7.0_DBD
NLSYADLITRAIQSSPDQRLTLSQIYEWMVRCVPYFNDKGDSNSSAGWKNSIRHNLSLHS
RFMRVQNEGTGKSSWWIINPDG
>Pteropus_alecto_Black_flying_fox|4.5_DBD
NLSYADLITRAIESSPDKRLTLSQIYEWMVRCVPYFKDKGDSNSSAGW
>Rousettus_aegyptiacus_Egyptian_rousette|9.0_DBD
NLSYADLITRAIESSPDKRLTLSQIYEWMVRCVPYFKDKGDSNSSAGWKNSIRHNLSLHS
RFMRVQNEGTGKSSWWIINPDG
>Myotis_brandtii_Brandts_bat|40.0_DBD
NSSAGWKNSIRHNLSLHSRFMRVQNEGTGKSSWWIINPDG
>Myotis_myotis_Mouse-eared_bat|15.5_DBD
NLSYADLITRAIQSSPDQRLTLSQIYEWMVRCVPYFNDKGDSNSSAGWKNSIRHNLSLHS
RFMRVQNEGTGKSSWWIINPDG
>Rhinolophus_ferrumequinum_Greater_horseshoe_bat|30.0_DBD
NLSYADLITRAIESSPDKRLTLSQIYEWMVRCVPYFKDKGDSNSSAGWKNSIRHNLSLHS
RFMRVQNEGTGKSSWWIINPDG
```


# MSA.py

The MSA file takes the {order}_foxo3_DBD.fasta files from Data_extraction.py, and aligns them using MAFFT (Kazutaka and Daron 2013). Then alignments are trimmed using TrimAl (in order to reduce noise) (Capella-Gutiérrez et al 2009). All orders are done at once, and their files are formatted {order}_DBD_trimmed.fasta

MAFFT alignment parameters
- DBD sequence: mafft --localpair --genafpair --maxiterate 1000 --reorder --adjustdirection --thread N infile *L-INS-i

I chose *L-INS-i for the DBD sequences because it's likely the most accurate, and since the DBD region is highly conserved I was certain alignment would be pretty clean (Kazutaka and Daron 2013). I also chose it as it is accuracy oriented, which I prioritized since I have very small sample and sequence sizes. 

*Libraries*

import os
import shutil
import subprocess

*Full run:*

```
python3 MSA.py
```

*CHIROPTERA Output files:*


- Chiroptera_DBD_trimmed.fasta
  
```
>Myotis_lucifugus_Little_brown_bat|7.0_DBD
NSSAGW
>Myotis_brandtii_Brandts_bat|40.0_DBD
NSSAGW
>Myotis_myotis_Mouse-eared_bat|15.5_DBD
NSSAGW
>Pteropus_alecto_Black_flying_fox|4.5_DBD
NSSAGW
>Rousettus_aegyptiacus_Egyptian_rousette|9.0_DBD
NSSAGW
>Rhinolophus_ferrumequinum_Greater_horseshoe_bat|30.0_DBD
NSSAGW
```

# Phyl_tree.py

I reused the code created by Osbourne and others, making adjustments so it applied to my code (2024). To create the phylogenetic tree it uses a distance based method called neighbor joining. I felt this was also good choice considering again, my smaller sequence and sample size.


*Libraries:*
library(ape)   
library(seqinr)
library(Biostrings)
library(phytools)
library(phangorn)

*RUN:*
```
python3 Phyl_tree.py
```

# res.py

This is another code I reused and adapted from Osbourne and others as it was applicable for my case. The code converts fasta files to arrays, then it calculates a position weight matrix from the array of aligned sequences. In the end it calculates background frequency of each residue over the aligned sequences, and the position score using position weight matrix and background frequency.

The code calculates the scores for every single order, and outputs {order}_res_score.csv.

*Libraries:*
import argparse
import math
import pandas as pd
from Bio import SeqIO
import re

*RUN:*
```
python3 res.py
```

# PGLS_model.py

Calculates the Phylogenetic Generalized Least Squares using the Brownian covariance matrix, Residue scores, and phylogenetic tree.

The code outputs the PGLS for every single order in files formatted "{order}_PGLS_results.csv"

*Libraries*
import pandas as pd
import numpy as np
from skbio import TreeNode
import statsmodels.api as sm

*RUN*
```
python3 PGLS_model.py
```

# basic download information #
Biopython:

```
pip install biopython
```

Requests:

```
pip install requests
```

MAFFT:

```
sudo apt-get install mafft
```

Trimal:
```
git clone https://github.com/inab/trimal.git

cd trimal
cd Source

make
```
HMMer:
```
sudo apt install hmmer
```

Phylogenic tree
```
pip install biopython scikit-bio
pip install matplotlib
```

**other sorting stuff this part is not organized at all**
```
sudo apt install -y python3 python3-pip

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
**CITATIONS**

Kazutaka K, Daron M. S. 2013. MAFFT Multiple Sequence Alignment Software Version 7: Improvements in Performance and Usability, Molecular Biology and Evolution. 30:4(772–780). https://doi.org/10.1093/molbev/mst010

Capella-Gutiérrez S, Silla-Martínez M, Gabaldón T. 2009. trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. Bioinformatics. 25:15(1972–1973). https://doi.org/10.1093/bioinformatics/btp348
