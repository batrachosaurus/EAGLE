import os
import sys
import urllib2
import gzip
import wget

from EAGLEdb.lib import get_links_from_html


def get_bacteria_from_ncbi(refseq_bacteria_link="https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria",
                           genbank_bacteria_link="https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria",
                           bactdb_dir="EAGLEdb/bacteria",
                           n_first_bact=None):
    try:
        os.makedirs(bactdb_dir)
    except OSError:
        print "bactdb directory exists"
    bacteria_list = []
    refseq_list = get_links_from_html(urllib2.urlopen(refseq_bacteria_link))
    genbank_list = get_links_from_html(urllib2.urlopen(genbank_bacteria_link))
    n = 0
    i = 0
    j = 0
    while i < len(refseq_list) or j < len(genbank_list):
        if n_first_bact:
            if n >= n_first_bact:
                break
            n += 1
        if genbank_list[j] < refseq_list[i]:
            try:
                bacteria_list.append(get_bacterium(genbank_bacteria_link, genbank_list[j], bactdb_dir, "genbank"))
            except:
                print "%s is not prepared: %s" % (genbank_list[j], sys.exc_info())
            j += 1
        else:
            try:
                bacteria_list.append(get_bacterium(refseq_bacteria_link, refseq_list[i], bactdb_dir, "refseq"))
            except:
                print "%s is not prepared: %s" % (refseq_list[i], sys.exc_info())
            i += 1
        if genbank_list[j] == refseq_list[i-1]:
            j += 1
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
    print bacterium_link
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
    print "got %s taxonomy" % bacterium_info["strain"]
    bacterium_info["16S_rRNA_file"] = get_16S_fasta(assemblies_list[-1] + "_rna_from_genomic.fna.gz",
                                                    db_dir,
                                                    bacterium_info["strain"])
    print "got %s 16S rRNA" % bacterium_info["strain"]
    return bacterium_info


def download_bacterium_files(bact_prefix, suffixes, download_dir="./"):
    if type(suffixes) is str:
        suffixes_list = [suffixes]
    else:
        suffixes_list = list(suffixes)
    for suffix in suffixes_list:
        file_link = None
        file_link = bact_prefix + suffix
        wget.download(file_link, out=download_dir)


def get_taxonomy(f_name, f_dir, remove_tax_f=True):
    family = None
    genus = None
    species = None
    strain = None
    org = False
    tax_list = []
    f_path = os.path.join(f_dir, f_name)
    f = gzip.open(f_path, 'rb')
    for line_ in f:
        line = None
        line = line_.decode("utf-8").strip()
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
    f.close()
    if remove_tax_f:
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


def get_16S_fasta(f_name, f_dir, strain, remove_rna_f=True):
    fasta_path = os.path.join(f_dir, strain + "_16S_rRNA.fasta")
    fasta_f = open(fasta_path, 'w')
    f_path = os.path.join(f_dir, f_name)
    rRNA = False
    seq_list = []
    f = gzip.open(f_path, 'rb')
    for line_ in f:
        line = None
        line = line_.decode("utf-8").strip()
        if not line: continue
        if line[0] == ">":
            if rRNA:
                fasta_f.write("".join(seq_list) + "\n")
                rRNA = False
            if "[product=16S ribosomal RNA]" in line:
                rRNA = True
                fasta_f.write(line + "\n")
        elif rRNA:
            seq_list.append(line)
    if rRNA:
        fasta_f.write("".join(seq_list) + "\n")
        rRNA = False
    fasta_f.close()
    f.close()
    if remove_rna_f:
        os.remove(f_path)
    return fasta_path


def get_families_dict(bacteria_list):
    pass
