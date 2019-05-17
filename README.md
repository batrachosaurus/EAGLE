# EAGLE - Essential and Advantageous Genes Location Explorer  
### Now it works only for bacterial genomes

## Requirements
MUSCLE  
MAFFT  
MSAProbs  
HMMER3  
EMBOSS + PHYLIPNEW (EMBASSY package)  
FastME 2.07  
Blast+  
KaKs_Calculator 2.0  
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
&nbsp; - psutil >= 5.6.1  

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
You can (recommended way) download the default database from [here](http://ma.fbb.msu.ru/loven-doo/EAGLE/EAGLEdb.tar.gz). The downloaded database should be placed into the workdir (This is not usable at all - will be uprades for it. To save the diskspace, a symlink to EAGLEdb directory can be created in each workdir with 'ln -s </path/to/extracted/EAGLEdb> </path/to/workdir/EAGLEdb>' command). Each database scheme placed in <db_name>_info.json files that are located in the archive root directory ('EAGLEdb').  
Other option is to build it from prepared lists of NCBI genomes:
```
eagle_db -dbt bacteria
```
The created database scheme (db_info.json) will be located in the databese directory (EAGLEdb/bacteria)
  
Also below is the instruction for building a database from NCBI if you do not like to use the default database or prepared lists (another option):  
1. Download assembly summary ([here](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt) is RefSeq assembly summary table for bacteria and 
[here](https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt) is Genbank assembly summary table for bacteria).  
   
2. Prepare genomes lists:
```
eagle_db.prepare_ncbi_summary <downloaded/summary/path> <prepared/genomes/list/path>
```
   
3. Build the database
```
eagle_db -dbt bacteria -igenbank <prepared/genomes/list/path>
```
  

All this commands can be run as Python functions: see below eagledb package reference  
  
### Run the analysis
WARNING: sequences names in input fasta file longer than 10 symbols may produce errors.  
  
Type the command below to start the analysis:
```
eagle -i <fasta/path> -db <EAGLEdb/scheme/json/path> -nt <threads_number> -o <out/dir/path>
```
for detailed parameters description type:
```
eagle -h
```
Also the analysis can be run from Python:
```
from eagle import explore_orfs

explore_orfs(...)
```
  
## Packages reference


### eagle

### eagledb

