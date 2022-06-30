import os
import pandas
from Bio import Entrez
import time
import re
import numpy as np

class EntrezClint:

    def __init__(self):
        Entrez.email = "your.email@domain.tld"

    def get_prot_seq(self,accession_n):
        result = Entrez.esearch(db = 'Nucleotide', term=accession_n)
        record = Entrez.read(result)
        gi = record["IdList"][0]
        gene_record = Entrez.efetch(db="Nucleotide", id=gi, retmode="xml")
        xml = Entrez.read(gene_record)
        asStr = str(xml[0])
        prot_start = '\'translation\', \'GBQualifier_value\': '
        prot_seq = asStr.split(prot_start)[1].split('}')[0]
        assert prot_seq[0] == '\'' and prot_seq[-1] == '\'', accession_n + " protein seq didnt fully accomplished " + prot_seq
        # Sleep for 0.3 seconds, to keep NCBI servers happy
        time.sleep(0.3)
        return (prot_seq[1:-1] + 'X')


class ZnfDomainsHits:
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
        self.domains_regions = []
        self.domains_span = 0
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
    # TODO check enricment to certain amino
    def is_binding_amino(self, amino_site):
        int(amino_site)
        if amino_site in self.binding_aminos:
            return True
        else:
            return False
    def which_of_forth(self, amino_site):
        #will save in which_of_forth the editing is
        sites = np.zeros([4],'int')
        for fourth in self.four_binding:
            for i, k in enumerate(fourth):
                if int(amino_site) == int(k):
                    sites[i] = 1
        return sites

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

    # list of gens with more the 10 editing sites
    genes = pandas.read_csv("/private8/Projects/itamar/ZNF/orshai_editing_sites/above_10_geneID.csv")
    # all editing sites
    sites = pandas.read_csv("/private8/Projects/itamar/ZNF/orshai_editing_sites/above_10_sites.csv")
    # create EntrezClint for getting proteins seq
    EzC = EntrezClint()
    # array to save number of editing in all 4 DNA binding amino acids
    total_fourth = np.zeros(4,'int')
    # save lengteh of all proteins
    proteins_len_all = 0
    # save editing in  DNA binding amino acids
    total_binds = 0
    # save number of DNA binding amino acids
    total_binding_sites = 0
    # save editing in ZNF domains
    total_in_domains = 0
    # save total length of all domains
    total_domain_length = 0
    # save list of changing  in ZNF domains - like "L182P"
    total_bind_change_list = []
    # take only none syn sites TODO - remove?
    sites = sites[sites['mut_type_short'] == 'nonsyn']
    # save lists of results for adding cols to df
    non_syn_sites = []
    prots_seq = []
    in_fingerprint = []
    fingerprint_expected_rate = []
    fingerprint_rate = []
    in_domain = []
    domain_expected_rate = []
    domain_rate = []
    aa_change_list = []
    for index, gene in genes.iterrows():

        prot_seq = EzC.get_prot_seq(gene['accession_n'])
        prots_seq.append(prot_seq)
        #take all gene's sites
        g_sites = sites[sites['gene_name'] == gene['gene_name']]
        non_syn_sites.append(g_sites.shape[0])
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
            assert orig_amino == prot_seq[amino_place] , orig_amino + prot_seq[amino_place]
            # check if the site fall on special position
            if hits_check.is_within_domain(amino_place):
                gene_in_domain_counter +=1
            if hits_check.is_binding_amino(amino_place):
                gene_bind_counter +=1
                gene_bind_change_list.append(change)
                fourth_sites_counter_gene += hits_check.which_of_forth(amino_place)

        in_fingerprint.append(gene_bind_counter)
        in_domain.append(gene_in_domain_counter)
        aa_change_list.append(";".join(gene_bind_change_list))
        #check if there was enrichment in binding
        # calculate expected rate
        gene_expected_aminos_rate = round(g_sites.shape[0]/len(prot_seq), 4)
        fingerprint_expected_rate.append(gene_expected_aminos_rate)
        # calculate resulted rate
        gene_binding_aminos_rate = round(gene_bind_counter/hits_check.get_statistcs()['num_of_aminos'] if hits_check.get_statistcs()['num_of_aminos'] > 0 else 0, 4)
        fingerprint_rate.append(gene_binding_aminos_rate)
        # check if there was enrichment in ZNF domain:
        # calculate expected rate
        gene_expected_in_domain_rate = hits_check.get_statistcs()['span']/len(prot_seq)
        domain_expected_rate.append(gene_expected_in_domain_rate)
        # calculate resulted rate
        gene_in_domain_rate=gene_in_domain_counter/g_sites.shape[0]
        domain_rate.append(gene_in_domain_rate)
        # sum for total
        proteins_len_all += len(prot_seq)
        total_binds += gene_bind_counter
        total_binding_sites += hits_check.get_statistcs()['num_of_aminos']
        total_fourth += fourth_sites_counter_gene
        total_in_domains += gene_in_domain_counter
        total_domain_length += hits_check.get_statistcs()['span']
        total_bind_change_list.extend(gene_bind_change_list)



    genes['non-syn-sites'] = non_syn_sites
    genes['prot_seq'] = prots_seq
    genes['in_fingerprint_count'] = in_fingerprint
    genes['in_fingerprint_expected_rate'] = fingerprint_expected_rate
    genes['in_fingerprint_rate'] = fingerprint_rate
    genes['in_domain_count'] = in_domain
    genes['in_domain_expected_rate'] = domain_expected_rate
    genes['in_domain_rate'] = domain_rate
    genes['change_list'] = aa_change_list

    genes.to_csv("/private8/Projects/itamar/ZNF/orshai_editing_sites/above_10_genes_result.csv", index=False)

    # check results for all genes:
    expected_aminos_rate_total = sites.shape[0]/proteins_len_all
    binding_aminos_rate_total = total_binds/total_binding_sites
    expected_aminos_count_total = sites.shape[0]*(total_binding_sites/proteins_len_all)

    #summerize in domains editing
    expected_in_domain_rate_total = total_domain_length/proteins_len_all
    expected_in_domain_count_total = sites.shape[0]*expected_in_domain_rate_total
    in_domain_rate_total = total_in_domains/sites.shape[0]


