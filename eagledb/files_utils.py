import json
import os
import sys

import pandas

from eagledb.scheme import GenomeInfo  # may produce import error


def are_bacteria_analyzed(bacteria_list_path=None, analyzed_bacteria_path=None):
    if not bacteria_list_path and not analyzed_bacteria_path:
        try:
            bacteria_list_path = sys.argv[1]
            analyzed_bacteria_path = sys.argv[2]
        except IndexError:
            print("Number of arguments must be 2: 1 - bacteria list json path; 2 - analyzed bacteria json path")
    bacteria_list_f = open(bacteria_list_path)
    bacteria_list = json.load(bacteria_list_f)
    bacteria_list_f.close()
    analyzed_bacteria = dict()
    for bacterium in bacteria_list:
        if os.path.exists(bacterium[u'16S_rRNA_file']):
            analyzed_bacteria[bacterium[u'strain'].replace("_", " ")] = True
    analyzed_bacteria_f = open(analyzed_bacteria_path, 'w')
    json.dump(analyzed_bacteria, analyzed_bacteria_f)
    analyzed_bacteria_f.close()
    return analyzed_bacteria


def prepare_summary_table(summary_table=None, out_table_path=None):
    if not summary_table and not out_table_path:
        try:
            summary_table = sys.argv[1]
            out_table_path = sys.argv[2]
        except IndexError:
            print("Number of arguments must be 2: 1 - summary table path; 2 - output table path")
    summary_df = pandas.read_csv(summary_table, header=1, sep="\t", dtype=str)
    lines_list = summary_df.apply(lambda df_row: {
        "org_name": df_row['organism_name'],
        "ncbi_link": df_row['ftp_path'],
        "repr": True if "repr" in df_row['refseq_category'].lower() else False,
        } if df_row['assembly_level'].lower() == "complete genome" else None, axis=1)
    prepared_df = pandas.DataFrame(filter(None, lines_list))
    prepared_df = prepared_df[["org_name", "ncbi_link", "repr"]]
    prepared_df.to_csv(out_table_path, sep="\t", index=False)


def join_genomes_list_files(genomes_list_1_path=None, genomes_list_2_path=None, joined_genomes_list_path=None):
    if not genomes_list_1_path and not genomes_list_2_path and not joined_genomes_list_path:
        try:
            genomes_list_1_path = sys.argv[1]
            genomes_list_2_path = sys.argv[2]
            joined_genomes_list_path = sys.argv[3]
        except IndexError:
            print("Error: number of arguments must be 3: 1 - genomes list 1 json path; 2 - genomes list 2 json path; "
                  "3 - output joined genomes list json path")

    genomes_list_1_f = open(genomes_list_1_path)
    genomes_list_1 = json.load(genomes_list_1_f)
    genomes_list_1_f.close()
    genomes_list_2_f = open(genomes_list_2_path)
    genomes_list_2 = json.load(genomes_list_2_f)
    genomes_list_2_f.close()

    joined_genomes_list = join_genomes_lists(genomes_list_1, genomes_list_2)
    joined_genomes_f = open(joined_genomes_list_path, "w")
    json.dump(joined_genomes_list, joined_genomes_f)
    joined_genomes_f.close()


def join_genomes_lists(genomes_list_1, genomes_list_2):
    joined_genome_info = set()
    joined_genomes_list = list()

    for genome_info in genomes_list_1:
        if not genome_info:
            continue
        g_info = GenomeInfo.load_from_dict(genome_info)
        joined_genomes_list.append(genome_info)
        joined_genome_info.add(g_info.genome_id)

    for genome_info in genomes_list_2:
        if not genome_info:
            continue
        g_info = GenomeInfo.load_from_dict(genome_info)
        if g_info.org_name not in joined_genome_info:
            joined_genomes_list.append(genome_info)
            joined_genome_info.add(g_info.genome_id)

    return joined_genomes_list
