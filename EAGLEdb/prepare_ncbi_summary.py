import sys
import pandas


def prepare_summary_tables(summary_table=None, out_table_path=None):
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


if __name__ == "__main__":
    prepare_summary_tables(sys.argv[1], sys.argv[2])
