# this code =
# input:
# sites with amino change and gene symbol
# genes list (symbol) with compatible acc_number
# note !! make sure your amino-change is compatible with its gene's acc_number
# output:
# will print and write to file statistics and detailed properties on the domain/fingerprints changes dou to the amino-changes


import os
import pandas as pd 
from Bio import Entrez
import time
import re
import numpy as np

class EntrezClint:

    def __init__(self):
        Entrez.email = "your.email@domain.tld"

    def get_prot_seq(self,accession_n):
        # get protein id by gene acc
        result = Entrez.esearch(db = 'Nucleotide', term=accession_n)
        record = Entrez.read(result)
        gi = record["IdList"][0]
        # get rwcord by id
        gene_record = Entrez.efetch(db="Nucleotide", id=gi, retmode="xml")
        xml = Entrez.read(gene_record)
        asStr = str(xml[0])
        prot_start = '\'translation\', \'GBQualifier_value\': '
        prot_seq = asStr.split(prot_start)[1].split('}')[0]
        assert prot_seq[0] == '\'' and prot_seq[-1] == '\'', accession_n + " protein seq didnt fully accomplished " + prot_seq
        # Sleep for 0.3 seconds, to keep NCBI servers happy
        time.sleep(0.3)
        # remove start an end "\" and add stop codon at the end
        return (prot_seq[1:-1] + 'X')


class ZnfDomainsHits:
    # will find and save domain details for the protein seq
    # will save:
    # fingerprints aminos location
    # group for every 4  fingerprints aminos
    # domains' start/end
    # domain total length
    # count domains
    def __init__(self,seq):
        #C2H2 pattern from prosite - regex
        c2h2_pat = "C.{2,4}C.{3}[LIVMFYWC].{8}H.{3,5}H"
        c2h2_comped = re.compile(c2h2_pat)
        #pattern for the region with all the binding aminos and with constant number of amino acids
        crit_area = "C.{3}[LIVMFYWC].{8}H"
        area_comp = re.compile(crit_area)
        # set of positions of aminos that bind dna
        self.binding_aminos = set()
        # list with groups of 4 positions of aminos that bind dna
        # yes its redundant to have both
        self.four_binding = []
        #wil save domains' start/end
        self.domains_regions = []
        # will sum all domain length
        self.domains_span = 0
        # will count domain
        self.num_domains = 0
        for hit in c2h2_comped.finditer(seq):
            # get hit's sequnce
            domain = hit.group()
            domain_start = hit.start()
            domain_end = domain_start+len(domain)
            #add (start,end) to regions and its size to span
            self.domains_regions.append(tuple([domain_start,domain_end]))
            self.domains_span += len(domain)
            self.num_domains +=1
            # get hit's critical area
            for area in area_comp.finditer(domain):
                # if its the right place dont look for next critical area
                if 2 < area.start() < 6:
                    # save the 'H' place
                    h_place = domain_start + area.start() + 13
                    self.binding_aminos.add(h_place - 1)
                    self.binding_aminos.add(h_place - 4)
                    self.binding_aminos.add(h_place - 5)
                    self.binding_aminos.add(h_place - 7)
                    self.four_binding.append(tuple([h_place - 7,h_place - 5,h_place - 4,h_place - 1]))

    # will return true/false if site is in binding amino site
    def is_binding_amino(self, amino_site):
        int(amino_site)
        if amino_site in self.binding_aminos:
            return True
        else:
            return False
   
    # will return nparray with 1 in the binding one
    def which_of_forth(self, amino_site):
        #will save in which_of_forth the editing is
        sites = np.zeros([4],'int')
        for fourth in self.four_binding:
            for i, k in enumerate(fourth):
                if int(amino_site) == int(k):
                    sites[i] = 1
        return sites
        
    # will return true/false if site is in binding domain
    def is_within_domain(self,amino_site):
        for region in self.domains_regions:
            if region[0] <= amino_site <= region[1]:
                return True
        return False
    def get_statistcs(self):
        return {'num_of_domains':self.num_domains,'num_of_aminos':len(self.binding_aminos),'span':self.domains_span}

# todo check recurrent amino replacment?
if __name__ == "__main__":
    # for writing results
    with open("/private8/Projects/itamar/ZNF/Sra_GSE73211_35_samples/REDI/curated/domain_editing/ZNF_binding_edit.txt", 'w') as outf:
        # list of gens with more the 10 editing sites
        genes = pd.read_csv("/private8/Projects/itamar/ZNF/Sra_GSE73211_35_samples/REDI/curated/domain_editing/above_40_geneID.csv")
        # all editing sites
        sites = pd.read_csv("/private8/Projects/itamar/ZNF/Sra_GSE73211_35_samples/REDI/curated/domain_editing/above_40_sites.csv")
        # remove stop loss/gain sites - - nothing to do with binding here
        sites = sites[~sites['amino_change'].str.contains("*",regex=False)]
        # create EntrezClint for getting proteins seq
        EzC = EntrezClint()
        # array to save number of editing in all 4 DNA binding amino acids - to check enrichment
        total_fourth = np.zeros(4,'int')
        # save lengteh of all proteins - for norelize
        proteins_len_all = 0
        # save editing in  DNA binding amino acids
        total_binds = 0
        # save number of DNA binding amino acids - (for normilizeing - how many binding aminos in all the genes)
        total_binding_sites = 0
        # save editing in ZNF domains
        total_in_domains = 0
        # save total length of all domains
        total_domain_length = 0
        # save list of changing  in ZNF domains - like "L182P"
        total_bind_change_list = []
        # take only none syn sites TODO - remove?
        sites = sites[sites['mut_type_short'] != 'syn']
        df_for_seq_mut_rows = []
        # check results for every gene
        for index, gene in genes.iterrows():
            prot_seq = EzC.get_prot_seq(gene['accession_n'])
            #take all gene's sites
            g_sites = sites[sites['gene_name'] == gene['gene_name']]
            #get set of places with binding-related amino acids:
            # create hits checker
            hits_check = ZnfDomainsHits(prot_seq)
            gene_bind_counter = 0
            gene_in_domain_counter = 0
            fourth_sites_counter_gene = np.zeros(4, 'int')
            # list of changing  in ZNF domains in gene
            gene_bind_change_list = []
            
            # check every editing site for the gene if its on binding-related amino acids
            for index, site in g_sites.iterrows():
                change = site['amino_change']  # like "L182P"
                orig_amino = change[0] # L
                new_amino = change[-1] # P
                #update place to -1 becouse orshai counts from 1 not from 0
                amino_place = int(change[1:-1]) - 1
                # check matching between orshai and our protein seq
                assert orig_amino == prot_seq[amino_place] , "orig: " + orig_amino + "seq :" + prot_seq[amino_place] + '\n' + str(site) + "\n\n" +prot_seq
                # check if the site fall on special position
                if hits_check.is_within_domain(amino_place):
                    gene_in_domain_counter +=1
                if hits_check.is_binding_amino(amino_place):
                    gene_bind_counter +=1
                    gene_bind_change_list.append(change)
                    fourth_sites_counter_gene += hits_check.which_of_forth(amino_place)
            #check if there was enrichment in binding
            # calculate expected rate
            gene_expected_aminos_rate = round(g_sites.shape[0]/len(prot_seq), 4)
            # calculate resulted rate
            gene_binding_aminos_rate = round(gene_bind_counter/hits_check.get_statistcs()['num_of_aminos'] if hits_check.get_statistcs()['num_of_aminos'] > 0 else 0, 4)

            # write output
            outf.write(gene['gene_name']+"\nnon-syn sites: "+ str(g_sites.shape[0]) + "\nbinding change: "+ str(gene_bind_change_list) + "\n")
            print("")
            print(gene['gene_name'], "non-syn sites: " ,g_sites.shape[0], "binding change: ", gene_bind_change_list)
            outf.write("editing binding aminos rate expected: {}, results: {}\n".format(gene_expected_aminos_rate,
                                                                                 gene_binding_aminos_rate))
            print("editing binding aminos rate expected: {}, results: {}".format(gene_expected_aminos_rate,gene_binding_aminos_rate))

            # check if there was enrichment in ZNF domain:
            # calculate expected rate
            gene_expected_in_domain_rate = hits_check.get_statistcs()['span']/len(prot_seq)
            # calculate resulted rate
            gene_in_domain_rate=gene_in_domain_counter/g_sites.shape[0]
            # write output
            outf.write("in ZNF-domain editing rate expected: {}, results: {}\n\n".format(gene_expected_in_domain_rate,gene_in_domain_rate))
            print("in ZNF-domain editing rate expected: {}, results: {}".format(gene_expected_in_domain_rate,gene_in_domain_rate))
            print("")
            # sum for total
            proteins_len_all += len(prot_seq)
            total_binds += gene_bind_counter
            total_binding_sites += hits_check.get_statistcs()['num_of_aminos']
            total_fourth += fourth_sites_counter_gene
            total_in_domains += gene_in_domain_counter
            total_domain_length += hits_check.get_statistcs()['span']
            total_bind_change_list.extend(gene_bind_change_list)



        # check results for all genes:
        expected_aminos_rate_total = sites.shape[0]/proteins_len_all
        binding_aminos_rate_total = total_binds/total_binding_sites
        expected_aminos_count_total = sites.shape[0]*(total_binding_sites/proteins_len_all)
        outf.write("total non-syn sites: {}, total sites in binding sites: {}, expected in binding site: {}\n".format(sites.shape[0], total_binds, round(expected_aminos_count_total, 4)))
        print("total non-syn non sites: {}, total sites in binding sits: {}, expected in binding site: {}".format(sites.shape[0], total_binds, round(expected_aminos_count_total, 4)))
        outf.write("total expected edit rate per point: {},total results: {}\n".format(round(expected_aminos_rate_total,4), round(binding_aminos_rate_total,4)))
        print("total expected edit rate per point: {},total results: {}".format(round(expected_aminos_rate_total,4), round(binding_aminos_rate_total,4)))

        #summerize in domains editing
        expected_in_domain_rate_total = total_domain_length/proteins_len_all
        expected_in_domain_count_total = sites.shape[0]*expected_in_domain_rate_total
        in_domain_rate_total = total_in_domains/sites.shape[0]
        outf.write("total sites in domains: {}, expected in domains: {}\n".format(total_in_domains, round(expected_in_domain_count_total)))
        print("total sites in domains: {}, expected in domains: {}".format(total_in_domains, round(expected_in_domain_count_total)))
        outf.write(("total expected rate of edit in domain(in domains/in all protein): {},total results: {}\n".format(round(expected_in_domain_rate_total, 4),round(in_domain_rate_total, 4))))
        print("total expected rate of edit in domain(in domains/in all protein): {},total results: {}".format(round(expected_in_domain_rate_total, 4),round(in_domain_rate_total, 4)))
        # summerize editing in 4 binding amino acid options
        outf.write("fourth pattern:\nresidues in positions[-1, +2, +3, +6] : count:"+str(total_fourth))
        print("fourth pattern:\nresidues in positions[-1, +2, +3, +6] : count:",total_fourth)
        # write list of all changes
        outf.write("all binding changes: "+ str(total_bind_change_list))
        print("all binding changes: ", total_bind_change_list)


