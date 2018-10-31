# EAGLE - Essential and Advantageous Genes Location Explorer  
Now it works only for bacterial genomes

## Requirements
MUSCLE  
HMMER3  
EMBOSS + PHYLIPNEW (EMBASSY package)  
FastME 2.07  
KaKs Calculator 2.0  
Blast+  
Redis  
Python 2.7  
### Python packages:  
&nbsp; - wget >= 3.2  
&nbsp; - pyaml >= 3.12  
&nbsp; - numpy >= 1.14.3  
&nbsp; - pandas >= 0.22.0  
&nbsp; - scipy >= 1.1.0  
&nbsp; - biopython >= 1.72  
&nbsp; - redis >= 2.10.6  

## How to use

### Install it:
```
pip install git+https://github.com/loven-doo/EAGLE.git --upgrade  
```
from dev branch:
```
pip install git+https://github.com/loven-doo/EAGLE.git@dev --upgrade
```

### Prepare the database
You can (recommended way) download the default database from [here]()  
Other option is to build it from prepared lists of NCBI genomes:
```
EAGLEdb -dbt bacteria
```
  
Also below is the instruction for building a database from NCBI if you do not like to use the default database or prepared lists (another option):  
1. Download assembly summary ([here](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt) is RefSeq assembly summary table for bacteria and 
[here](https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt) is Genbank assembly summary table for bacteria).  
   
2. Prepare genomes lists:
```
EAGLEdb.prepare_ncbi_summary <downloaded/summary/path> <prepared/genomes/list/path>
```
   
3. Build the database
```
EAGLEdb -dbt bacteria -igenbank <prepared/genomes/list/path>
```
  

All this commands can be run as Python functions: see below EAGLEdb package reference  
  
### Run the analysis
```
EAGLE ...
```
or from Python
```
from EAGLE import ...


...(...)
```
  
## Packages reference

### EAGLE

### EAGLEdb

