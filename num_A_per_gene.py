import csv
import pandas as pd
from collections import defaultdict

path_symb = "/private8/Projects/itamar/ZNF/Sra_GSE73211_35_samples/REDI/curated/R_results/ZNF_start_end_name.csv"
path_As = "/private8/Projects/itamar/ZNF/Sra_GSE73211_35_samples/REDI/curated/input/Final_file_All_ZNF_A_RefseqCur_knownFormat_sorted.bed"

out_path = "/private8/Projects/itamar/ZNF/Sra_GSE73211_35_samples/REDI/curated/R_results/A_count_per_gene.csv"

genes_symb = pd.read_csv(path_symb, sep="\t")
genes_symb["cdsStart"] = genes_symb["cdsStart"].astype("int64", copy=False)
genes_symb["cdsEnd"] = genes_symb["cdsEnd"].astype("int64", copy=False)

genes_counter = defaultdict(int)

# for every site in file
with open(path_As, "r") as A_sites:
    DictReader_obj = csv.DictReader(
        A_sites, delimiter="\t", fieldnames=["chrom", "positions", "strand"]
    )
    for site in DictReader_obj:
        site["positions"] = int(site["positions"])
        # find genes symbol
        site_genes_symb = genes_symb[
            (genes_symb["chrom"] == site["chrom"])
            & (genes_symb["strand"] == site["strand"])
            & (genes_symb["cdsStart"] <= site["positions"])
            & (genes_symb["cdsEnd"] >= site["positions"])
        ]
        for genesymb in site_genes_symb.name2.unique():
            genes_counter[genesymb] += 1


d = dict(genes_counter)
table = pd.DataFrame.from_dict(d, orient="index", columns=["count"])
table = table.reset_index(level=0)
table.rename(columns={"index": "gene_name"}, inplace=True)
table.to_csv(out_path, index=False, sep="\t")
