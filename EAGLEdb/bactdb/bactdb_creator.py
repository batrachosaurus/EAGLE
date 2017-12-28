import os
import urllib2
import subprocess

from EAGLEdb.lib import get_links_from_html


def get_bacteria_from_ncbi(refseq_bacteria_link="https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria",
                           genbank_bacteria_link="https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria",
                           bactdb_dir="EAGLEdb/bacteria"):
    subprocess.call("mkdir -p " + bactdb_dir, shell=True)
    bacteria_list = []
    refseq_list = get_links_from_html(urllib2.urlopen(refseq_bacteria_link))
    genbank_list = get_links_from_html(urllib2.urlopen(genbank_bacteria_link))
    i = 0
    j = 0
    while True:
        if genbank_list[j] < refseq_list[i]:
            bacteria_list.append(get_bacterium(genbank_bacteria_link, genbank_list[j], bactdb_dir, "genbank"))
            j += 1
        elif genbank_list[j] == refseq_list[i]:
            bacteria_list.append(get_bacterium(refseq_bacteria_link, refseq_list[i], bactdb_dir, "refseq"))
            i += 1
            j += 1
        else:
            bacteria_list.append(get_bacterium(refseq_bacteria_link, refseq_list[i], bactdb_dir, "refseq"))
            i += 1
        if i == len(refseq_list) and j == len(genbank_list):
            break
    return bacteria_list


def get_bacterium(ncbi_db_link, bacterium_name, db_dir, source_db=None):
    bacterium_info = {"family": None,
                      "genus": None,
                      "species": None,
                      "strain": None,
                      "download_prefix": None,
                      "16S_rRNA_file": None,
                      "source_db": source_db,
                      "repr": False}
    bacterium_link = ncbi_db_link + "/" + bacterium_name
    bacterium_list = get_links_from_html(urllib2.urlopen(bacterium_link))
    if "representative" in bacterium_list:
        next_page = bacterium_link + "/" + "representative"
        bacterium_info["repr"] = True
    else:
        next_page = bacterium_link + "/" + "latest_assembly_versions"
    assemblies_list = get_links_from_html(urllib2.urlopen(next_page))
    bacterium_prefix = (next_page + "/" + assemblies_list[-1] + "/" + assemblies_list[-1]).replace("https", "ftp")
    bacterium_info["download_prefix"] = bacterium_prefix
    download_bacterium_files(bacterium_prefix, ["_wgsmaster.gbff.gz", "_rna_from_genomic.fna.gz"], db_dir)
    bacterium_info["family"], bacterium_info["genus"], bacterium_info["species"], bacterium_info["strain"] = \
        get_taxonomy(assemblies_list[-1] + "_wgsmaster.gbff.gz", db_dir)
    bacterium_info["16S_rRNA_file"] = get_16S_fasta(assemblies_list[-1] + "_rna_from_genomic.fna.gz",
                                                    db_dir,
                                                    bacterium_info["strain"])
    return bacterium_info


def download_bacterium_files(bact_prefix, suffixes, download_dir="./"):
    if type(suffixes) is str:
        suffixes_list = [suffixes]
    else:
        suffixes_list = list(suffixes)
    for suffix in suffixes_list:
        file_link = None
        file_link = bact_prefix + suffix
        subprocess.call("wget -P " + download_dir + "/ " + file_link, shell=True)


def get_taxonomy(f_name, f_dir):
    family = None
    genus = None
    species = None
    strain = None
    org = False
    tax_list = []
    f_path = os.path.join(f_dir, f_name)
    f = open (f_path)
    for line_ in f:
        line = None
        line = line_.strip()
        if not line: continue
        if line[:9] == "REFERENCE":
            family = get_family(tax_list, genus, species, strain)
            break
        if line[:8] == "ORGANISM":
            org = True
            line_list = line.split()
            genus = line_list[1]
            species = genus + "_" + line_list[2]
            strain = "_".join(line_list[1:])
        elif org:
            tax_list += list(prepare_tax_line(line))
    os.remove(f_path)
    return family, genus, species, strain


def get_family(tax_list, g, sp, st):
    fam = None
    n = -1
    while -n <= len(tax_list):
        tax_u = None
        tax_u = tax_list[n].replace(" ", "_")
        n = n - 1
        if tax_u == st or tax_u == sp or tax_u == g:
            continue
        else:
            fam = tax_u
            break
    return fam


def prepare_tax_line(tax_line):
    tax_line_list = tax_line.split(";")
    for elm_ in tax_line_list:
        elm = None
        elm = elm_.strip(" .\t")
        if elm: yield elm


def get_16S_fasta(f_name, f_dir, strain):
    fasta_path = os.path.join(f_dir, strain + ".fasta")
    fasta_f = open(fasta_path, 'w')
    f_path = os.path.join(f_dir, f_name)
    rRNA = False
    seq_list = []
    f = open(f_path)
    for line_ in f:
        line = None
        line = line_.strip()
        if not line: continue
        if line[0] == ">":
            if rRNA:
                fasta_f.write("".join(seq_list) + "\n")
                break
            if "[product=16S ribosomal RNA]" in line:
                fasta_f.write(line + "\n")
        elif rRNA:
            seq_list.append(line)
    if not rRNA:
        print ("No 16 rRNA has been found")
    fasta_f.close()
    os.remove(f_path)
    return fasta_path


def get_families_dict(bacteria_list):
    pass
