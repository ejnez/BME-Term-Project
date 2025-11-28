# BME-Term-Project

**Code descriptions (libraries, output, and running)**

Summary: 
FOXO3 sequences are extracted from Uniprot, Ensembl, AND NCBI, in that order.
If a sequence isn't found in Uniprot it's then searched in Ensembl, and if not found there
it searches NCBI last. All NCBI sequences are inputted in BLAST to check for orthology using human FOXO3 sequence.

The sequences are then split into 3 files, the full file has the entire sequences, the DBD file
has only the DNA binding domain, using the coordinates given by the human FOXO3 sequence in uniprot (157-251).
Finally the IDR file has all of the other (likely disordered regions).

The MSA file takes the sequence files from Data_extraction.py, aligns them using MAFFT (Kazutaka and Daron 2013). Then alignments are trimmed using TrimAl (in order to reduce noise) (Capella-Gutiérrez et al 2009). The test run will output every file (for manually checking quality), so there will be the alignment files before trimming, and the full final sequencs after alignment and trimming. The full run will only output final sequences.

MAFFT alignment parameters
- DBD sequence: mafft --localpair --genafpair --maxiterate 1000 --reorder --adjustdirection --thread N infile *L-INS-i
- IDR sequence and full sequence:  --genafpair --maxiterate 1000 --ep 0 --adjustdirection --thread N infile *E-INS-i

I chose *L-INS-i for the DBD sequences because it's likely the most accurate, and since the DBD region is highly conserved I was certain alignment would be pretty clean (Kazutaka and Daron 2013). For the whole sequence, and IDR sequence I decided to use *E-INS-i with --ep 0 because there are a lot of disordered regions which might contribute to excess gaps (Kazutaka and Daron 2013). Both of these are accuracy oriented, which I prioritized since I only have a sample size of 73 species (Kazutaka and Daron 2013). 

I utilized HMMer just to verify orthology, and if alignments are clean. I ran them for both the DBD and full sequence, keeping in mind that the disordered regions in the full sequence might affect results. I only saved the file for my DBD HMMer and posted it below. From this file I observed that EFFN was very low (2.55) for a sample size of 72, this reveals that the sequences are all very similar. This confirms (I think) that my sequences are very conserved, orthologous and likely correctly aligned. My full HMM_build file output, had an EFFN = 0.681702, for 73 sequences and was probably affected by IDR regions.

I replicated the code created by Osbourne and others in python to create my phylogenetic tree. I use a distance based method called neighbor joining. I felt this was a good choice considering again, my smaller 73 sequence sample size.

The res.py code from the paper on p53 seems applicable to my case, so I used it to create a csv file of residue specific scores.

Then I will get a file of normalized lifespan. Lifespan, RES file, and Phylogenetic tree will be used with a PGLS code for a final analysis. (hopefully theres something).


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
python3 Dataset_extraction.py {order}
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

I replicated the code created by Osbourne and others in python because I'm more comfortable with it. To create my phylogenetic tree I use a distance based method called neighbor joining. I felt this was a good choice considering again, my smaller sequence and sample size.


*Libraries:*
- AlignIO
- DistanceCalculator
- DistanceTreeConstructor
- Phylo
- matplotlib.pyplot

*RUN:*
```
python3 Phyl_tree.py
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
