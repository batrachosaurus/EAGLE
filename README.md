# EAGLE - Essential and Advantageous Genes Location Explorer  
### Now it works only for bacterial genomes

## Requirements
MUSCLE  
HMMER3  
EMBOSS + PHYLIPNEW (EMBASSY package)  
FastME 2.07  
Blast+  
Python 2.7  
### Python packages:  
&nbsp; - wget >= 3.2  
&nbsp; - pyaml >= 3.12  
&nbsp; - numpy >= 1.14.3  
&nbsp; - pandas == 0.22.0  
&nbsp; - matplotlib >= 2.2.3  
&nbsp; - scipy >= 1.1.0  
&nbsp; - DendroPy >= 4.4.0  
&nbsp; - biopython >= 1.72  
&nbsp; - redis >= 2.10.6  

## How to use

### Pull the docker image
```
docker pull loven7doo/eagle
```
The workdir in this image is '/EAGLE' so to run all commands from the image use this way:
```
docker run -v </host/system/workdr/location/>:/EAGLE <any command>
```
Do not forget to replace the path to the host system workdir with '/EAGLE' in commands

If docker is not the appropriate way follow steps below (requirements installation and the package installation)

### Instal the requirements
It can be very difficult to install some requirements on Windows. Linux is recommended to use.

### Install the package:
```
pip install git+https://github.com/loven-doo/EAGLE.git --upgrade  
```
from dev branch:
```
pip install git+https://github.com/loven-doo/EAGLE.git@dev --upgrade
```

### Prepare the database
You can (recommended way) download the default database from [here](http://ma.fbb.msu.ru/loven-doo/EAGLE/EAGLEdb.tar.gz)  
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
NOTE: sequences names in input fasta file cannot be longer than 10 symbols.  
  
Type the command below to start the analysis:
```
EAGLE -i <fasta/path> -db <EAGLEdb/scheme/json/path> -m <run_mode> -nt <threads_number> -o <out/dir/path>
```
for detailed parameters description type:
```
EAGLE -h
```
or from Python
```
from EAGLE import explore_genes

explore_genes(...)
```
  
## Packages reference

### EAGLE

### EAGLEdb

