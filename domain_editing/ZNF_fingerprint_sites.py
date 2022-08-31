import pandas
# !!!!!!!!!!!! very bad code - you should use another one!!!!!!!!!!!!!
# 
# will take list of site and will write the old and the new protein seq for every protein
if __name__ == "__main__":
    sites = pandas.read_csv("/private8/Projects/itamar/ZNF/orshai_editing_sites/above_10_sites.csv")
    genes_binding = pandas.read_csv("/private8/Projects/itamar/ZNF/orshai_editing_sites/above_10_genes_result.csv")

    fingerprint_sites = sites[sites['gene_name'].isin(genes_binding['gene_name'])]

    change_list = ['S701P', 'L677P', 'L509P', 'Y346H', 'F46S', 'F270L', 'F270S', 'S271P', 'L274P', 'Y454H', 'L426P',
                   'Y342H', 'S454P', 'S342P', 'L340S', 'Q352R', 'R358G', 'K411E', 'N467D', 'R470G', 'K536E', 'I533V',
                   'Q527R']

    fingerprint_sites = fingerprint_sites[fingerprint_sites['amino_change'].isin(change_list)]
    old_prot_seq = []
    new_prot_seq = []
    for index, site in fingerprint_sites.iterrows():
        old_seq = genes_binding[genes_binding['gene_name'] == site['gene_name']]
        old_seq = old_seq['prot_seq'].tolist()[0]
        change = site['amino_change'] # like "L182P"
        orig_amino = change[0]  # L
        new_amino = change[-1]  # P
        # update place to -1 becouse orshai counts from 1 not from 0
        amino_place = int(change[1:-1]) - 1
        new_seq = list(old_seq)
        assert new_seq[amino_place] == orig_amino, "not comaptible"
        new_seq[amino_place] = new_amino
        old_prot_seq.append(old_seq)
        new_prot_seq.append(''.join(new_seq))
    fingerprint_sites['old_prot_seq'] = old_prot_seq
    fingerprint_sites['new_prot_seq'] = new_prot_seq

    print("hey")
    
    #fingerprint_sites.to_csv("/private8/Projects/itamar/ZNF/orshai_editing_sites/in_fingerprints_sites.csv")
