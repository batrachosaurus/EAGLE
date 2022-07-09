import os
import shutil
import json
import zlib
import traceback
import urllib.request

import pandas as pd
from Bio import Entrez
from Bio import SearchIO
from Bio.SearchIO.HmmerIO.hmmer3_text import Hmmer3TextParser

from eaglib._utils.logging import eagle_logger
from eaglib.alignment import SeqsProfileInfo, SeqsProfile
from eaglib.seqs import SeqsDict, load_fasta_to_dict


def get_taxonomy(tax_id):
    tax_keys = ["superkingdom", "phylum", "clade", "class", "order", "family", "genus", "species"]
    tax_dict = {tax_key: None for tax_key in tax_keys}
    
    record = Entrez.efetch(db="taxonomy", id=tax_id, retmode='xml')
    tax_info = Entrez.read(record)[0]
    tax_dict["species"] = tax_info['ScientificName']
    
    for lin_tax in tax_info['LineageEx']:
        if lin_tax['Rank'] in tax_dict:
            tax_dict[lin_tax['Rank']] = lin_tax['ScientificName']
            
    return [tax_dict[tax_key] for tax_key in tax_keys]


ssu_rrna_profile = SeqsProfile(SeqsProfileInfo.load_from_dict({
    'name': 'SSU_rRNA_bacteria',
    'path': 'SSU_rRNA_bacteria.cm',
    'type': 'rna',
    'weight': 1.0,
    'method': 'infernal'
}))

hsp70_profile = SeqsProfile(SeqsProfileInfo.load_from_dict({
    'name': 'HSP70',
    'path': 'HSP70.hmm',
    'type': 'protein',
    'weight': 1.0,
    'method': 'hmmer'
}))


assembly_summary_path = "bacteria_assembly_summary.txt"
Entrez.email = "moshenskydenis@gmail.com"
db_dir = "bacteria"

processed_ac = list()
bact_df = pd.read_csv(assembly_summary_path, sep="\t", low_memory=False).query("assembly_level=='Complete Genome' & refseq_category!='representative genome'")
eagle_logger.info(" %s genomes to prepare" % len(bact_df))
for _, row in bact_df.iterrows():
    ac = row['assembly_accession']  # id field in genomes_table
    asm = row['asm_name']
    taxonomy = get_taxonomy(row['species_taxid'])
    name = row['organism_name'] + ("" if pd.isna(row['infraspecific_name']) else " " + row['infraspecific_name'])
    ftp_prefix = (row['ftp_path'] + "/" + ac + "_" + asm).replace(" ", "_")
    fna_seq = [ftp_prefix+"_genomic.fna.gz"]
    btc_seqs = [os.path.join(db_dir, ac+"_btc.fasta")]
    rna_path = os.path.join(db_dir, ac+"_rna_from_genomic.fna")
    tcds_path = os.path.join(db_dir, ac+"_translated_cds.faa")
    btc_seqs_dict = dict()
    
    try:
        with open(rna_path, 'wb') as rna_f:
            rna_f.write(zlib.decompress(urllib.request.urlopen(ftp_prefix+"_rna_from_genomic.fna.gz").read(), 15+32))
        psr_rna_df = ssu_rrna_profile.search(seqdb=rna_path, threads=4)
        if not psr_rna_df.empty:
            max_score_rna = psr_rna_df.loc[psr_rna_df["score"].idxmax()]
            btc_seqs_dict[ssu_rrna_profile.name] = load_fasta_to_dict(fasta_path=rna_path)[max_score_rna["target name"]][max_score_rna["seq from"]-1: max_score_rna["seq to"]]                     
        
        with open(tcds_path, 'wb') as tcds_f:
            tcds_f.write(zlib.decompress(urllib.request.urlopen(ftp_prefix+"_translated_cds.faa.gz").read(), 15+32))
        psr_hsp_df = hsp70_profile.search(seqdb=tcds_path, threads=4)
        if not psr_hsp_df.empty:
            max_score_hsp = psr_hsp_df.loc[psr_hsp_df["domain_score"].idxmax()]
            btc_seqs_dict[hsp70_profile.name] = load_fasta_to_dict(fasta_path=tcds_path)[max_score_hsp["target name"]][max_score_hsp["ali_from"]-1: max_score_hsp["ali_to"]]                     
        
        SeqsDict.load_from_dict(btc_seqs_dict).dump(btc_seqs[0])
        processed_ac.append({"id": ac, "name": name, "taxonomy": taxonomy, "btc_seqs": btc_seqs, "fna_seq": fna_seq})
    except:
        print(traceback.format_exc())
    eagle_logger.info(" " + ac + "\t" + name)
    # if _ >= 10: break ###

with open("processed_bacteria.json", "w") as processed_ac_f:
    json.dump(processed_ac, processed_ac_f, indent=2)
