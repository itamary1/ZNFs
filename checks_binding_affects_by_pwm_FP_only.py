# this code will get a table of genes and their amino change
# 1. will create a folder with fasta for every gen version and its pwm prdiect
# 2. will save into file differece in pwm prediction between protein versions
# every pwm proccess will be with polenoimial and linear pwm prediction
# note! it will be much faster to run genes in parallell,  but pwm_predict have some bugs which prevent it

import multiprocessing
import pandas as pd
import os
import subprocess


# global constant parametrs
DIFF_THRESH = 0.2
pwm_pol_script_path = "/home/alu/twerski/pwm_predict/run_pwm_one_file_pol.sh"
pwm_lin_script_path = "/home/alu/twerski/pwm_predict/run_pwm_one_file_lin.sh"
parent_dir = (
    "/private8/Projects/itamar/ZNF/Sra_GSE73211_35_samples/REDI/curated/domain_editing/"
)


def runPwmPredict(fasta_path):
    p1 = multiprocessing.Process(target=runPwmPredictPol, args=[fasta_path])
    p1.start()
    p2 = multiprocessing.Process(target=runPwmPredictLin, args=[fasta_path])
    p2.start()
    p1.join()
    p2.join()


def runPwmPredictPol(fasta_path):
    subprocess.run(["sh", pwm_pol_script_path, fasta_path + "_pol.fasta"])


def runPwmPredictLin(fasta_path):
    subprocess.run(["sh", pwm_lin_script_path, fasta_path + "_lin.fasta"])


def compare_one_result(origial_prot_fasta_path, new_prot_fasta_path, suffix):
    with open(origial_prot_fasta_path + suffix, "r") as FDorig_pwm:
        orig_pwm = FDorig_pwm.readlines()
    with open(new_prot_fasta_path + suffix, "r") as FDnew:
        new_pwm = FDnew.readlines()
    assert len(orig_pwm) == len(new_pwm), " diffrent pwn len"
    diffs = [0.0, 0]
    for i in range(len(orig_pwm)):
        # if we come to the binding prediction
        if orig_pwm[i].startswith(">"):
            orig_list = [
                orig_pwm[i + 1].split(" ")[:-1],
                orig_pwm[i + 2].split(" ")[:-1],
                orig_pwm[i + 3].split(" ")[:-1],
                orig_pwm[i + 4].split(" ")[:-1],
            ]
            new_list = [
                new_pwm[i + 1].split(" ")[:-1],
                new_pwm[i + 2].split(" ")[:-1],
                new_pwm[i + 3].split(" ")[:-1],
                new_pwm[i + 4].split(" ")[:-1],
            ]
            # get tuple of sum diff and count in the binding prediction
            res = prob_compare(orig_list, new_list)
            diffs[0] += res[0]
            diffs[1] += res[1]
    return diffs


def prob_compare(orig_list, new_list):
    diff_counts = 0
    total_diff = 0
    # iterate over all nucleotides positions
    for i in range(len(orig_list[0])):
        nuc_diff = 0
        # iterate over 4 nucleotids and sum difference
        for orig_nuc_list, new_nuc_list in zip(orig_list, new_list):
            try:
                float(orig_nuc_list[i])
            except ValueError:
                print("Not a float" + orig_nuc_list[i])
            nuc_diff += abs(float(orig_nuc_list[i]) - float(new_nuc_list[i]))
            total_diff += nuc_diff
        # if sum differance of 4 nucs in same position is more then threshold - increase counter
        if nuc_diff > DIFF_THRESH:
            diff_counts += 1
    return [total_diff, diff_counts]


def compare_pwm_result(origial_prot_fasta_path, new_prot_fasta_path):
    poll_diff = compare_one_result(
        origial_prot_fasta_path, new_prot_fasta_path, "_pol.pwm"
    )
    lin_diff = compare_one_result(
        origial_prot_fasta_path, new_prot_fasta_path, "_lin.pwm"
    )
    return {"pol_diff": poll_diff, "lin_diff": lin_diff}


if __name__ == "__main__":
    genes_result = pd.read_csv(parent_dir + "above_40_genes_result.csv")
    pwm_dir = parent_dir + "pwm_result/"
    # will contain list for of lists for every gene - all_gene_diff[one_gene_diff[one_ver_diff[]]]
    pol_diff = []
    lin_diff = []
    diffs_counters_pol = []
    diffs_counters_lin = []
    num_sequneces = []
    for index, gene in genes_result.iterrows():
        gene_pol_diff = []
        gene_lin_diff = []
        gene_seqs_count = 0
        gene_count_pol = 0
        gene_count_lin = 0
        mut_list = gene["in_fingerprints_list"]
        if pd.isna(mut_list):
            diffs_counters_pol.append(0)
            diffs_counters_lin.append(0)
            num_sequneces.append(0)
            pol_diff.append([])
            lin_diff.append([])
            continue
        gene_dir = gene["gene_name"]
        gene_path = os.path.join(pwm_dir, gene_dir)
        # os.mkdir(gene_path)
        seq = gene["prot_seq"]
        # write first fasta for origina seq
        origial_prot_fasta_path = os.path.join(gene_path, (gene["gene_name"]))
        # with open(origial_prot_fasta_path+'_pol.fasta', 'w') as FPfile:
        #     FPfile.write('>' + str(gene['gene_name']) + '\n' + str(seq) + '\n')
        # with open(origial_prot_fasta_path+'_lin.fasta', 'w') as FLfile:
        #     FLfile.write('>' + str(gene['gene_name']) + '\n' + str(seq) + '\n')
        # runPwmPredict(origial_prot_fasta_path)
        # write fasta for every seq
        # its wrote like Q352R;R358G
        mut_list = mut_list.split(";")
        for mutation in mut_list:
            mutable_seq = list(seq)
            # update amino changew to new seq (-1 beca)
            mutable_seq[int(mutation[1:-1]) - 1] = mutation[-1]
            if "*" in mutation:
                continue
            # make sure the original amino compatible with our change list
            assert seq[int(mutation[1:-1]) - 1] == mutation[0], (
                "ERROR amino place" + mutation + gene["gene_name"]
            )
            gene_ver_name = str(gene["gene_name"]) + "_" + mutation
            new_prot_fasta_path = os.path.join(gene_path, gene_ver_name)
            # with open(new_prot_fasta_path+'_pol.fasta', 'w') as nFPfile:
            #     nFPfile.write('>' + gene_ver_name + '\n' + "".join(mutable_seq) + '\n')
            # with open(new_prot_fasta_path+'_lin.fasta', 'w') as nFLfile:
            #     nFLfile.write('>' + gene_ver_name + '\n' + "".join(mutable_seq) + '\n')
            # runPwmPredict(new_prot_fasta_path)
            # compare binding prediction of original and edited protein
            ver_diff = compare_pwm_result(origial_prot_fasta_path, new_prot_fasta_path)
            # print(ver_diff)
            gene_count_pol += ver_diff["pol_diff"][1]
            gene_count_lin += ver_diff["lin_diff"][1]
            if ver_diff["pol_diff"][1] > 0 or ver_diff["lin_diff"][1] > 0:
                gene_seqs_count += 1
            gene_pol_diff.append("; ".join([str(x) for x in ver_diff["pol_diff"]]))
            gene_lin_diff.append("; ".join([str(x) for x in ver_diff["lin_diff"]]))
        # add gene's diffs to all genes' diffs list
        diffs_counters_pol.append(gene_count_pol)
        diffs_counters_lin.append(gene_count_lin)
        num_sequneces.append(gene_seqs_count)
        pol_diff.append(gene_pol_diff)
        lin_diff.append(gene_lin_diff)
    genes_result["pol_diff"] = pol_diff
    genes_result["pwm_diffs_count_Polynomial"] = diffs_counters_pol
    genes_result["lin_diff"] = lin_diff
    genes_result["pwm_diffs_count_linear"] = diffs_counters_lin
    genes_result["num_sequneces"] = num_sequneces
    out_path = parent_dir + "finel_table_above_40_genes_pwm_result_FP_only.csv"
    genes_result.to_csv(path_or_buf=out_path, index=False)
