"""
Script created Wednesday, October 9th 2019

@uthor: BP Baranowsky
"""
import os

#################################################################################
""" Importing libraries that will be used. """
from collections import Counter
from scipy.stats import wilcoxon as test
import pandas as pd
from statsmodels.stats.multitest import multipletests as correction  #
import datetime
import json
import math
import pdb
from scipy.stats import poisson
import psycopg2

#################################################################################


class NeighborhoodAnalyzerFromGene:

    def __init__(self, user_list_of_genes, user_distance_value, user_database, user_cutoff, user_correction,
                 user_strand_value,
                 user_output):
        self.user_list_of_genes = user_list_of_genes
        self.user_distance_value = user_distance_value
        self.user_database = user_database
        self.user_cutoff = user_cutoff
        self.user_correction = user_correction
        self.user_strand_value = user_strand_value
        self.user_output = user_output
        self.cursor = self.connect_to_database()

    def connect_to_database(self):
        conn = psycopg2.connect(
            host="172.17.0.1",
            database="mydatabase",
            user="postgres",
            password="postgres")
        cursor = conn.cursor()
        return cursor

    def open_single_column_file(self, file_name):
        """ Opens file that contains one column of data """
        plik = []
        with open(file_name) as inputfile:
            for line in inputfile:
                plik.append(line.strip())
        return plik

    def open_multiple_column_file(self, file_name, split_mark=" "):
        """ Opens file that contain more than one column and split it by space. """
        plik = []
        with open(file_name) as inputfile:
            for line in inputfile:
                plik.append(line.strip().split(split_mark))
        return plik

    def open_singleline_num(self, path_to_file):
        """ Opens file that contain 1 column and strip it by space. """

        with open(path_to_file, 'r') as f:
            file_names = [float(line.strip()) for line in f]
        return file_names

    def open_json_file(self, file_name):
        with open(file_name, 'r') as fp:
            data = json.load(fp)
        all_data = []
        for key, values in data.items():
            all_data.append([key, values[0], values[1]])
        return all_data

    def open_query(self, path_to_file):

        with open(path_to_file, 'r') as f:
            file_names = [line.strip() for line in f]
        return file_names

    def save_complete(self, file, complete_data):
        """ Saves all stuff together as one file """
        complete_data.index.name = 'Domain'
        complete_data.to_csv(file, sep='\t', mode='w')

    def collect_pfam_domain_description(self, file):
        """ Opens file only to gather information"""

        data = {}
        with open(file) as inputfile:
            lines  = [x.strip().split('\t') for x in inputfile]
        for line in lines:
            try:

                domena = line[0]
                delete_this = "(" + domena.upper() + ")"
                family = line[1].replace(delete_this, "")
                summary = line[2]
                data[domena] = (family, summary)

            except IndexError:
                continue
        return data

    def create_6_lists(self, file):

        # data = self.open_multiple_column_file(file)
        start_coord = []
        end_coord = []
        orientation = []
        domains = []
        genes = []
        contig = []
        # NOWA BAZA DANYCH, KLASTER
        # with open("nowa_baza_danych/"+file) as inputfile:
        # VULTR
        with open("modules/data/img_ready/" + file) as inputfile:
            for line in inputfile:
                bit = line.strip().split()

                genes.append(int(bit[0]))
                start_coord.append(int(bit[1]))
                end_coord.append(int(bit[2]))
                orientation.append(bit[3])
                domains.append(bit[4])
                contig.append(bit[5])
        return genes, start_coord, end_coord, orientation, domains, contig

    def find_correct_genome(self, gene, reference_data):

        for line in reference_data:
            genome = line[0]
            starts = [x for x in line[1::2]]
            ends = [x for x in line[2::2]]
            counter = 0
            for start in starts:
                end = ends[counter]
                counter += 1
                if int(gene) in range(int(start), int(end) + 1):
                    return genome

    def find_correct_genome_sql(self, gene):
        query = f"SELECT genome, begins,ends FROM genomes_map WHERE {int(gene)} BETWEEN begins AND ends"
        self.cursor.execute(query)
        records = self.cursor.fetchall()
        if len(records) != 0:
            genome = records[0][0]
            return genome



    def get_genome_size_in_gene(self, genome_id, list_of_genome, list_of_size):
        """ Takes genome's id, list of genome and lists with data about how much genes in genome"""
        indeks_genomu = list_of_genome.index(genome_id)
        genome_size_in_gene = list_of_size[indeks_genomu]
        return genome_size_in_gene

    def get_genome_size_in_gene_sql(self,gene):
        query = f" SELECT size_in_genes FROM genome_sizes WHERE img_id = '{gene}' "
        self.cursor.execute(query)
        records = self.cursor.fetchall()
        if len(records) != 0 :
            genome_size_in_gene = records[0][0]
            return genome_size_in_gene


    def search_for_gene_in_genome(self, q_gene, genes, start_coords, end_coords, orients, contig):
        result = []
        main_index = genes.index(int(q_gene))
        result.append(start_coords[main_index])
        result.append(end_coords[main_index])
        result.append(orients[main_index])
        result.append(contig[main_index])
        result.append(q_gene)
        return result

    def get_range_coordinates(self, target_gene, distance):
        """
        Based on how large neighbourhood user wants to analyze creates
        points FROM and TO, additionaly shows where main user's pfam begins and ENDS
        with orientation on strands
        """

        gene_id = target_gene[4]
        gene_beg = target_gene[0]
        gene_end = target_gene[1]
        searched_gene_orient = target_gene[2]
        last_coordinate = gene_end + int(distance)
        first_coordinate = gene_beg - int(distance)

        return last_coordinate, first_coordinate, gene_beg, gene_end, searched_gene_orient

    def get_domain_both_strands(self, start_coord, end_coord, last_coordinate, first_coordinate, gene_beg,
                                gene_end, contig, result):
        pfam_index_to_neigh = []
        s_coord_counter = 0
        for s_coord in start_coord:
            if s_coord <= last_coordinate and s_coord >= gene_beg and contig[s_coord_counter] == result[
                3] and s_coord_counter not in pfam_index_to_neigh:
                pfam_index_to_neigh.append(s_coord_counter)
                s_coord_counter += 1
            else:
                s_coord_counter += 1
                continue
        e_coord_counter = 0
        for e_coord in end_coord:
            if e_coord >= first_coordinate and e_coord <= gene_end and contig[e_coord_counter] == result[
                3] and e_coord_counter not in pfam_index_to_neigh:
                pfam_index_to_neigh.append(e_coord_counter)
                e_coord_counter += 1
            else:
                e_coord_counter += 1
                continue
        return pfam_index_to_neigh

    def get_domain_same_strand(self, start_coord, end_coord, last_coordinate, first_coordinate, gene_beg,
                               gene_end, contig, result):
        pass

    def get_domain_oposite_strand(self):
        pass

    def get_domain_both_zeros(self, ref_gene, q_gene):

        index_list_gene = []
        counter = 0
        for gene in ref_gene:
            if gene == int(q_gene):
                index_list_gene.append(counter)
            counter += 1
        return index_list_gene

    def get_list_of_genes_in_neigh(self, list_of_genes, index_list):
        """ Takes overall gene's list in genome and creates complete list of genes in neighbourhood """

        list_of_genes_in_neigh = []
        for pfam_domain in index_list:
            list_of_genes_in_neigh.append(list_of_genes[pfam_domain])
        return list_of_genes_in_neigh

    def get_size_in_genes(self, genome_or_gene_list):
        """ Takes list of genes and returns size of neighbourhood """

        list_of_genes = []
        for gene in genome_or_gene_list:
            list_of_genes.append(gene)
        return len(set(list_of_genes))

    def get_list_of_domain_in_neigh(self, pfam_index, gene, domains):

        to_counter = []
        party = []
        party.append(gene)
        for part in pfam_index:
            party.append(domains[part])
            to_counter.append(domains[part])
        #        neighbourhood_complete.append(party)
        # pdb.set_trace()
        return to_counter

    def multiple_test_correction(self, some_dataframe, correction_met):

        if correction_met == 'none':
            return some_dataframe
        else:
            try:

                value_to_correct = [float(x) for x in some_dataframe.Poisson.tolist()]
                reject, pvals_corrected, alphaSidak, alphaBonf = correction(pvals=value_to_correct,
                                                                            method=correction_met,
                                                                            is_sorted=False, returnsorted=False)
                pvals = pvals_corrected.tolist()
                some_dataframe.Poisson = pvals
                return some_dataframe
            except ZeroDivisionError:
                pass

    def multiple_test_correction_list(self, some_tuple_list, correction_met):
        pvalues = [float(x) for x in some_tuple_list]

        if correction_met == 'none':
            output = pvalues

        else:
            try:

                reject, pvals_corrected, alphaSidak, alphaBonf = correction(pvals=pvalues,
                                                                            method=correction_met,
                                                                            is_sorted=False, returnsorted=False)
                pvals_after_corr = pvals_corrected.tolist()
                output = pvals_after_corr
            except ZeroDivisionError:
                output = pvalues
        return output

    def cutoff_value(self, some_dataframe, cutoff):

        domain_list = list(some_dataframe.index)
        if cutoff == 'none':
            return some_dataframe
        elif cutoff == '0':
            for i in domain_list:
                diff = float(some_dataframe.loc[i, 'Poisson'])
                if diff < 0 and diff > 0:
                    some_dataframe = some_dataframe.drop([i])
            return some_dataframe
        else:
            cutoff = float(cutoff)
            for i in domain_list:
                diff = float(some_dataframe.loc[i, 'Poisson'])
                if diff > cutoff:
                    some_dataframe = some_dataframe.drop([i])
            return some_dataframe

    def sort_table(self, some_dataframe):

        sorted_data = some_dataframe.sort_values('Poisson', ascending=True)
        return sorted_data

    def add_information(self, some_dataframe, dictionary):

        indeksy = list(some_dataframe.index)
        for i in indeksy:
            try:

                pfam = i[0:2] + i[4:]
                family = dictionary[pfam][0]
                summary = dictionary[pfam][1]

                some_dataframe.at[i, 'Family'] = family
                some_dataframe.at[i, 'Summary'] = summary
            except KeyError:
                continue
        return some_dataframe

    def pfam_for_pf(self, dataframe):
        """Final customization of dataframe
            --> PF02696 instead of pfam02696
            --> PVALUE in 10e-3 format
            --> In what percentage % format
            --> avg occurences and density  have now 3 digits after coma
            --> drop NOam
            -->
        """

        indexy = dataframe.index
        for i in indexy:
            dataframe = dataframe.rename(index={i: i[:2].upper() + i[4:]})
        dataframe['Poisson'] = dataframe['Poisson'].map('{:.3e}'.format)
        dataframe['In what percentage'] = dataframe['In what percentage'].map('{:.3%}'.format)
        dataframe = dataframe.round({"average occurence in neighbourhood": 3, "average occurence in genome": 3,
                                     "Density difference": 3})
        if "NOam" in indexy:
            dataframe = dataframe.drop(index="NOam")
        # dataframe = dataframe.style.format({'In what percentage': '{:,.3%}'.format})
        dataframe.drop('average occurence in neighborhood', axis=1, inplace=True)
        dataframe.drop('average occurence in genome', axis=1, inplace=True)
        dataframe.drop('avg neigh size', axis=1, inplace=True)
        # dataframe.drop('domain in whole database', axis=1, inplace=True)
        # if self.user_correction == 'none':
        #     # dataframe.drop('Poisson_correct', axis=1, inplace=True)
        # elif :
        #     dataframe.drop('Poisson', axis=1, inplace=True)

        return dataframe

    def collect_gene_without_genome(self, gene):
        with open('modules/data/data_to_download.txt',
                  'a+') as output:
            output.write(gene)
            output.write("\n")

    def make_genera_statistic(self, genome_list):

        with open("/home/djangoadmin/final_site-project/important_files/genera_statistics", "r") as handler:
            data = [x.strip().split() for x in handler]  # data to lista list [[genom,genus],[genome,genus]]
        genomes = [x[0] for x in data]
        genera = [x[1] for x in data]

        just_genera = []

        for genome in genome_list:
            indeks = genomes.index(genome)
            just_genera.append(genera[indeks])
        counter_genera = Counter(just_genera)

        zliczanie = []
        for genus, value in counter_genera.items():
            zliczanie.append((genus, value / len(genome_list)))
        return zliczanie

    def poisson_distribution(self, p, n, domain, percentage_counter):

        if int(percentage_counter[domain]) < p * n:
            return 1
        else:
            sums = 0
            k = 0
            number_of_samples = int(percentage_counter[domain])
            while k < number_of_samples:
                sums += poisson.pmf(k=k, mu=n * p)
                k += 1
            if 1 - sums >= 0.0:
                return 1 - sums
            else:
                return 0.0

            # stat = poisson.pmf(k=int(percentage_counter[domain]), mu=n * p)
            # return stat

    def go(self):

        query_data = self.open_query("modules/data/input")

        domain_information = self.collect_pfam_domain_description(
            "modules/data/domain_information")

        genome_gene_reference = self.open_multiple_column_file(
            "modules/data/genomes_map_new")
        genome_size_in_gene = self.open_multiple_column_file(
            "modules/data/GENOME_ID_SIZE_IN_GENE_NEW.txt")
        with open('modules/data/counter_domains_new', 'r') as handler:
            pfam_counter_domains = json.load(handler)

        print("Checkpoint #1 Data loaded", datetime.datetime.now())
        genome_ids = [x[0] for x in genome_size_in_gene]
        genome_size = [x[1] for x in genome_size_in_gene]
        genome_size_in_gene = None
        just_neigh_data = []
        genomes_with_domains = []
        neigh_genome_size = []
        counter = 0
        counter_db = 0
        sum_of_genes_in_neighborhoods = 0
        number_of_neighborhoods = 0
        genome_number_to_stat = 0
        ################################################################################################################
        genome_counter = Counter()
        ################################################################################################################
        print("DLOK TIME START", datetime.datetime.now())
        # Loop through list of genes and do all stuff

        percentage_counter = Counter()
        for gene in query_data:
            counter_db += 1
            counter += 1
            try:

                # correct_genome = self.find_correct_genome(gene, genome_gene_reference)
                correct_genome = self.find_correct_genome_sql(gene)
                correct_size = self.get_genome_size_in_gene(correct_genome, genome_ids, genome_size)
                if correct_genome is not None:
                    # correct_size = self.get_genome_size_in_gene_sql(correct_genome)
                    ref_gene, ref_start_coord, ref_end_coord, ref_orient, ref_domain, ref_contig = self.create_6_lists(
                        correct_genome)
                    ##################################################################
                    genome_domains_counter = Counter(ref_domain)
                    genome_counter += genome_domains_counter
                    #################################################################
                    search_result = self.search_for_gene_in_genome(gene, ref_gene, ref_start_coord, ref_end_coord,
                                                                   ref_orient,
                                                                   ref_contig)
                    if self.user_distance_value != 0:
                        last_coordiate, first_coordinate, gene_beg, gene_end, searched_gene_orientation, = \
                            self.get_range_coordinates(search_result, self.user_distance_value)
                        pfam_index_to_neigh = self.get_domain_both_strands(ref_start_coord, ref_end_coord, last_coordiate,
                                                                           first_coordinate, gene_beg, gene_end, ref_contig,
                                                                           search_result)

                    else:
                        pfam_index_to_neigh = self.get_domain_both_zeros(ref_gene, search_result[4])

                    genes_in_neigh = self.get_list_of_genes_in_neigh(ref_gene, pfam_index_to_neigh)
                    number_of_genes_in_neigh = self.get_size_in_genes(genes_in_neigh)
                    whole_neigh = self.get_list_of_domain_in_neigh(pfam_index_to_neigh, gene, ref_domain)
                    just_neigh_data.append(whole_neigh)
                    genomes_with_domains.append(correct_genome)
                    neigh_genome_size.append((number_of_genes_in_neigh, correct_size))
                    percentage_counter += Counter(whole_neigh)
                    genome_number_to_stat += 1
                    sum_of_genes_in_neighborhoods += int(number_of_genes_in_neigh)
                    number_of_neighborhoods += 1
            except ValueError:
                self.collect_gene_without_genome(gene)
                continue
        print("DLOK TIME OVER", datetime.datetime.now())
        print("concatenate_domain")
        alls = []
        neigh_counter = Counter()

        for i in just_neigh_data:
            for j in i:
                if j not in alls:
                    alls.append(j)
        # genome_counter = Counter()

        print("MATRIX CREATING", datetime.datetime.now())
        counter = 0
        for domain_list, genome_in_domain, ng_size in zip(just_neigh_data, genomes_with_domains, neigh_genome_size):
            counter += 1
            # genes, start_coord, end_coord, orientation, domains, contig = self.create_6_lists(genome_in_domain)
            # genome_domains_counter = Counter(domains)
            # genome_counter += genome_domains_counter
            neigh_domains_counter = Counter(domain_list)
            neigh_counter += neigh_domains_counter

        print("MATRIX OVER", datetime.datetime.now())

        print("Wilcoxon start", datetime.datetime.now())

        najlepsze_dane = pd.DataFrame(columns=['occurence in neighborhoods',
                                               'average occurence in neighborhood', 'occurence genomes',
                                               'average occurence in genome', 'In what percentage', 'Poisson',
                                               'avg neigh size', 'Family', 'Summary'])

        counter = 0
        scores = []
        # percentage = []
        scores_poisson = []
        for i in alls:
            counter += 1
            # data_to_calculate = self.open_singleline_num("/home/djangoadmin/final_site-project/scripts/temp_data/" + i)
            # for_percentage = [float(x) for x in data_to_calculate if x > 0]
            # percentage.append(len(for_percentage) / len(just_neigh_data))
            # wynik = test(data_to_calculate, zero_method="wilcox")
            # mean = sum(data_to_calculate) / len(just_neigh_data)
            # scores.append((wynik.pvalue, mean))

            p_KP = ((sum_of_genes_in_neighborhoods / genome_number_to_stat) * pfam_counter_domains[i]) / 326628067

            # najlepsze_dane.at[i, 'P'] = p_KP
            najlepsze_dane.at[i, 'avg neigh size'] = sum_of_genes_in_neighborhoods / genome_number_to_stat
            najlepsze_dane.at[i, 'domain in whole database'] = pfam_counter_domains[i]
            # najlepsze_dane.at[i, 'sum of database len'] = 326628067
            poisson_test_result = self.poisson_distribution(p_KP, number_of_neighborhoods, i, percentage_counter)
            # najlepsze_dane.at[i, 'K'] = int(percentage_counter[i])
            # najlepsze_dane.at[i, 'N'] = number_of_neighborhoods
            # najlepsze_dane.at[i, 'P'] = p_KP
            najlepsze_dane.at[i, 'Poisson'] = poisson_test_result
            scores_poisson.append(poisson_test_result)
        scores_after_correction = self.multiple_test_correction_list(scores_poisson, self.user_correction)
        print("Wilcoxon OVER", datetime.datetime.now())

        print("Collecting data START", datetime.datetime.now())


        for domain, wynik in zip(najlepsze_dane.index.to_list(), scores_after_correction):
            najlepsze_dane.at[domain, 'Poisson_correct'] = wynik
            najlepsze_dane.at[domain, 'In what percentage'] = percentage_counter[domain] / len(just_neigh_data)
            najlepsze_dane.at[domain, 'occurence genomes'] = genome_counter[domain]
            najlepsze_dane.at[domain, 'average occurence in genome'] = genome_counter[domain] / len(
                set(genomes_with_domains))
            najlepsze_dane.at[domain, 'occurence in neighborhoods'] = neigh_counter[domain]
            najlepsze_dane.at[domain, 'average occurence in neighborhood'] = neigh_counter[domain] / len(
                just_neigh_data)

        print("Collecting data OVER", datetime.datetime.now())

        print("Customizing_data START", datetime.datetime.now())
        # after_test = self.multiple_test_correction(najlepsze_dane, self.user_correction)
        after_cutoff = self.cutoff_value(najlepsze_dane, self.user_cutoff)
        sorted_table = self.sort_table(after_cutoff)
        added_information = self.add_information(sorted_table, domain_information)
        final_data = self.pfam_for_pf(added_information)
        print("Customizing_data OVER", datetime.datetime.now())

        self.save_complete(self.user_output, final_data)


class SuperSpeedAnalysisFromDomain:

    def __init__(self, user_pfam, user_distance, user_organisms, user_cutoff, user_correction,
                 user_strand, user_output, skip_negative):
        """ Initializing input data
            Inputs:
                    user_pfam:  str (from pfam00000 to pfam99999)
                    user_distance: str non-negative integer 1-20000
                    user_organisms: str
                    user_cutoff: str (none, 0, 0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.005,
                                 0.01, 0.05, 0.1, 0.2, 0.5)
                    user_correction: str (none, bonferroni, fdr_bh)
                    user_strand: str (both)
                    user_output: str
         """
        self.user_pfam = user_pfam
        self.user_distance = user_distance
        self.user_organisms = user_organisms
        self.user_cutoff = user_cutoff
        self.user_correction = user_correction
        self.user_strand = user_strand
        self.user_output = user_output
        self.skip_negative = skip_negative

    def open_singleline(self, path_to_file):
        """ Opens file that contain 1 column and strip it by space.
            Inputs:
                    path_to_file: str
        """

        with open(path_to_file, 'r') as f:
            file_names = [line.strip() for line in f]
        return file_names

    def open_singleline_num(self, path_to_file):
        """ Opens file that contain 1 column and strip it by space. """

        with open(path_to_file, 'r') as f:
            file_names = [float(line.strip()) for line in f]
        return file_names

    def open_multiple_line(self, path_to_file):
        """ Opens file that contains data about single genome """

        plik = []
        with open(path_to_file) as inputfile:
            for line in inputfile:
                plik.append(line.strip().split())
        return plik

    def save_data(self, file_name, message_for_output, message_for_additional_data, message_down, complete_data):
        """ Saves all stuff together as one file """


        complete_data.index.name = 'Domain'
        complete_data.to_csv(file_name, sep='\t', mode='w')

    def create_six_list(self, file):
        """ Open single genome file chew it and return 6 lists -> GENE, START_COORD,
        END_COORDS ,ORIENTATION, DOMAINS """

        data = self.open_multiple_line("modules/data/img_ready/" + file)
        start_coord = []
        end_coord = []
        orientation = []
        domains = []
        genes = []
        contig = []
        for bit in data:
            genes.append(int(bit[0]))
            start_coord.append(int(bit[1]))
            end_coord.append(int(bit[2]))
            orientation.append(bit[3])
            domains.append(bit[4])
            contig.append(bit[5])
        return genes, start_coord, end_coord, orientation, domains, contig

    def open_database(self, tax):

        tax = tax.replace("_", " ")
        if tax == "all genomes":

            # LOCAL
            # with open("/mnt/d/45.76.38.24/final_site-project/important_files/all_genomes.txt","r") as handler:

            # VULTR
            genomes = os.listdir('modules/data/img_ready')
            return genomes
        elif tax == "entero":
            with open("/home/djangoadmin/final_site-project/important_files/entero", "r") as handler:
                genomes = [x.strip() for x in handler]
            return genomes



        else:
            # LOCAL
            # with open('/mnt/d/45.76.38.24/final_site-project/important_files/database.json', "r") as handler:

            # VULTR
            with open("modules/data/database.json", "r") as handler:
                data = json.load(handler)
            if tax in data.keys():
                genomes = []
                for genus, genome in data[tax].items():
                    genomes += genome
                return genomes
            else:
                for family, genera in data.items():
                    if tax in genera:
                        genomes = data[family][tax]
                        return genomes

    def genome_size_in_gene(self, genome_id, list_of_genome, list_of_size):
        """ Takes genome's id, list of genome and lists with data about how much genes in genome"""

        index_genome = list_of_genome.index(genome_id)
        genome_size_in_gene = list_of_size[index_genome]
        return genome_size_in_gene

    def size_in_genes(self, genome_or_genes):
        """ Takes list of genes and returns size of neighbourhood """

        list_of_genes = []
        for gene in genome_or_genes:
            list_of_genes.append(gene)
        return len(set(list_of_genes))

    def list_of_genes_in_neigh(self, list_of_genes, index_list):
        """ Takes overall gene's list in genome and creates complete list of genes in neighbourhood """

        list_of_genes_in_neigh = []
        for pfam_domain in index_list:
            list_of_genes_in_neigh.append(list_of_genes[pfam_domain])
        return list_of_genes_in_neigh

    def searching_for_domain_in_genome(self, pfam, start_coord, end_coord, orient, domains, contig, genes):
        """
        Takes user's pfam domain searches through lists contains coordinates,
        orientation, domains, contigs and returns complete data about all
        users' domains in genome
        """
        coords = []
        coords_counter = 0
        for domain in domains:
            one_coords_part = []
            if domain == pfam:
                orientation_pfam_domain = orient[coords_counter]
                one_coords_part.append(start_coord[coords_counter])
                one_coords_part.append(end_coord[coords_counter])
                one_coords_part.append(orientation_pfam_domain)
                one_coords_part.append(contig[coords_counter])
                one_coords_part.append(genes[coords_counter])
                coords.append(one_coords_part)
                coords_counter += 1
            else:
                coords_counter += 1
                continue
        return coords

    def presence_confirmation(self, coords, file, pfam):
        """
                Just prints out information how many pfam domains is in each genome
                during analysis propably i will have to remove it before putting it to website ...
                """

        if len(coords) == 0:
            print("In genome  " + file + " pfam domain you have been searching do not exist")
        elif len(coords) == 1:
            print("In genome  " + file + " there is " + str(len(coords)) + " " + pfam + " domain")
        else:
            print("In genome  " + file + " there are " + str(len(coords)) + " " + pfam + " domains")

    def get_range_coordinates(self, target_pfam, distance):
        """
        Based on how large neighbourhood user wants to analyze creates
        points FROM and TO, additionaly shows where main user's pfam begins and ENDS
        with orientation on strands
        """

        last_coordinate = target_pfam[1] + int(distance)
        first_coordinate = target_pfam[0] - int(distance)
        pfam_beg = target_pfam[0]
        pfam_end = target_pfam[1]
        searched_pfam_orientation = target_pfam[2]
        return last_coordinate, first_coordinate, pfam_beg, pfam_end, searched_pfam_orientation

    def get_domain_both_zeros(self, ref_gene, q_gene):

        index_list_gene = []
        counter = 0
        for gene in ref_gene:
            if gene == int(q_gene):
                index_list_gene.append(counter)
            counter += 1
        return index_list_gene

    def get_domains_both_strands(self, start_coord, end_coord, last_coordinate, first_coordinate, pfam_beg,
                                 pfam_end, contig, point):

        pfam_index_to_neigh = []
        s_coord_counter = 0
        for s_coord in start_coord:
            if s_coord <= last_coordinate and s_coord >= pfam_beg and contig[s_coord_counter] == point[
                3] and s_coord_counter not in pfam_index_to_neigh:
                pfam_index_to_neigh.append(s_coord_counter)
                s_coord_counter += 1
            else:
                s_coord_counter += 1
                continue
        e_coord_counter = 0
        for e_coord in end_coord:
            if e_coord >= first_coordinate and e_coord <= pfam_end and contig[e_coord_counter] == point[
                3] and e_coord_counter not in pfam_index_to_neigh:
                pfam_index_to_neigh.append(e_coord_counter)
                e_coord_counter += 1
            else:
                e_coord_counter += 1
                continue
        return pfam_index_to_neigh

    def collect_pfam_domain_description(self, file):
        """ Opens file only to gather information"""

        data = {}
        with open(file) as inputfile:
            for line in inputfile:
                try:

                    one_line = line.strip().split("@")
                    domena = one_line[0]
                    pre_family = one_line[1][8:]
                    delete_this = "(" + domena.upper() + ")"
                    family = pre_family.replace(delete_this, "")
                    summary = one_line[2][9:]
                    data[domena] = (family, summary)

                except IndexError:
                    continue
        return data

    def get_domains_plus_minus_strand(self, start_coord, end_coord, last_coordinate, first_coordinate,
                                      pfam_beg, pfam_end, contig, point, orientation):

        pfam_index_to_neigh_same_strand = []
        pfam_index_to_neigh_oposite_strand = []
        s_coord_counter = 0
        for s_coord in start_coord:
            if s_coord <= last_coordinate and s_coord >= pfam_beg and orientation[s_coord_counter] == point[2] and \
                    contig[s_coord_counter] == point[3] and s_coord_counter not in pfam_index_to_neigh_same_strand:
                pfam_index_to_neigh_same_strand.append(s_coord_counter)
                s_coord_counter += 1
            elif s_coord <= last_coordinate and s_coord >= pfam_beg and orientation[s_coord_counter] != point[2] and \
                    contig[s_coord_counter] == point[
                3] and s_coord_counter not in pfam_index_to_neigh_oposite_strand:
                pfam_index_to_neigh_oposite_strand.append(s_coord_counter)
                s_coord_counter += 1
            else:
                s_coord_counter += 1
                continue
        e_coord_counter = 0
        for e_coord in end_coord:
            if e_coord >= first_coordinate and e_coord <= pfam_end and orientation[e_coord_counter] == point[2] and \
                    contig[e_coord_counter] == point[3] and e_coord_counter not in pfam_index_to_neigh_same_strand:
                pfam_index_to_neigh_same_strand.append(e_coord_counter)
                e_coord_counter += 1
            elif e_coord >= first_coordinate and e_coord <= pfam_end and orientation[e_coord_counter] != point[
                2] and contig[e_coord_counter] == point[
                3] and e_coord_counter not in pfam_index_to_neigh_oposite_strand:
                pfam_index_to_neigh_oposite_strand.append(e_coord_counter)
                e_coord_counter += 1
            else:
                e_coord_counter += 1
                continue
        return pfam_index_to_neigh_same_strand, pfam_index_to_neigh_oposite_strand

    def get_list_of_domain_in_neigh(self, pfam_index, file, domains):

        to_counter = []
        party = []
        party.append(file)
        for part in pfam_index:
            party.append(domains[part])
            to_counter.append(domains[part])
        #        neighbourhood_complete.append(party)

        return to_counter

    def counting_dlok_or_dglob(self, some_counter, genome_neigh_size, mighty_domains):
        dlok_glob = []
        for domain_mighty in mighty_domains:
            if domain_mighty in some_counter.keys():
                lok_glob = some_counter.get(domain_mighty) / int(genome_neigh_size)
                dlok_glob.append(lok_glob)
            else:
                dlok_glob.append(0)
        return dlok_glob

    def count_for_dlok(self, somecounter, neigh_size, tax_name, dlok_dataframe):
        for k, v in somecounter.items():
            dlok_dataframe.at[tax_name, k] = v / neigh_size

    def multiple_test_correction(self, some_dataframe, correction_met):

        if correction_met == 'none':
            return some_dataframe
        else:
            value_to_correct = [float(x) for x in some_dataframe.PVALUE.tolist()]
            reject, pvals_corrected, alphaSidak, alphaBonf = correction(pvals=value_to_correct,
                                                                        method=correction_met,
                                                                        is_sorted=False, returnsorted=False)
            pvals = pvals_corrected.tolist()
            some_dataframe.PVALUE = pvals
            return some_dataframe

    def multiple_test_correction_list(self, some_tuple_list, correction_met):
        pvalues = [float(x) for x in some_tuple_list]

        if correction_met == 'none':
            output = pvalues

        else:
            try:

                reject, pvals_corrected, alphaSidak, alphaBonf = correction(pvals=pvalues,
                                                                            method=correction_met,
                                                                            is_sorted=False, returnsorted=False)
                pvals_after_corr = pvals_corrected.tolist()
                output = pvals_after_corr
            except ZeroDivisionError:
                output = pvalues
        return output

    def cutoff_value(self, some_dataframe, cutoff):

        domain_list = list(some_dataframe.index)
        if cutoff == 'none':
            return some_dataframe
        elif cutoff == '0':
            for i in domain_list:
                diff = float(some_dataframe.loc[i, 'Poisson_correct'])
                if diff < 0 and diff > 0:
                    some_dataframe = some_dataframe.drop([i])
            return some_dataframe
        else:
            cutoff = float(cutoff)
            for i in domain_list:
                diff = float(some_dataframe.loc[i, 'Poisson_correct'])
                if diff > cutoff:
                    some_dataframe = some_dataframe.drop([i])
            return some_dataframe

    def sort_table(self, some_dataframe):

        sorted_data = some_dataframe.sort_values('Poisson', ascending=True)
        return sorted_data

    def add_information(self, some_dataframe, dictionary):

        indeksy = list(some_dataframe.index)
        for i in indeksy:
            try:

                pfam = i[0:2] + i[4:]
                family = dictionary[pfam][0]
                summary = dictionary[pfam][1]

                some_dataframe.at[i, 'Family'] = family
                some_dataframe.at[i, 'Summary'] = summary
            except KeyError:
                continue
        return some_dataframe

    def poisson_distribution_KP(self, k, p, n):
        """
        Return value of poisson distribution mine concept
            Parameters:
                k (int): sum of domain occurrence in all neighborhoods
                n (int): number of neighborhoods
                p (float): ( average lengh of neighborhood * number of domain in all database) / sum of lengh of all database
            How to aggregate data:
                k Counter of all domains
                n - number of neigghborhoods (easy +=)
                p - average lengh of neighborhood ( easy +=) / number of genomes, number of domain in all ( from file), sum of lenght of all database ( const value)
        """
        first_value = (n * p) ** k
        second_value = math.e ** (-1 * n * p)
        third_value = math.factorial(k)
        return (first_value * second_value) / third_value

    def poisson_distribution_KP_mod(self, k, p, n):
        """
        Return value of poisson distribution mine concept
            Parameters:
                k (int): sum of domain occurrence in all neighborhoods
                n (int): number of neighborhoods
                p (float): ( average lengh of neighborhood * number of domain in all database) / sum of lengh of all database
            How to aggregate data:
                k Counter of all domains
                n - number of neigghborhoods (easy +=)
                p - average lengh of neighborhood ( easy +=) / number of genomes, number of domain in all ( from file), sum of lenght of all database ( const value)
        """

        if n * p < k:

            first_value = (n * p) ** k
            second_value = math.e ** (-1 * n * p)
            third_value = math.factorial(k)
            return (first_value * second_value) / third_value
        else:
            return 1.0

    def pfam_for_pf(self, dataframe):
        """Final customization of dataframe
            --> PF02696 instead of pfam02696
            --> PVALUE in 10e-3 format
            --> In what percentage % format
            --> avg occurences and density  have now 3 digits after coma
            --> drop NOam
        """

        indexy = dataframe.index
        for i in indexy:
            dataframe = dataframe.rename(index={i: i[:2].upper() + i[4:]})
        dataframe['Poisson'] = dataframe['Poisson'].map('{:.3e}'.format)
        dataframe['Poisson_correct'] = dataframe['Poisson_correct'].map('{:.3e}'.format)

        dataframe['In what percentage'] = dataframe['In what percentage'].map('{:.3%}'.format)
        dataframe['average occurence in neighborhood'] = dataframe['average occurence in neighborhood'].map(
            '{:.3}'.format)
        dataframe['average occurence in genome'] = dataframe['average occurence in genome'].map('{:.3}'.format)

        if "NOam" in indexy:
            dataframe = dataframe.drop(index="NOam")
        # dataframe = dataframe.style.format({'In what percentage': '{:,.3%}'.format})
        dataframe.drop('average occurence in neighborhood', axis=1, inplace=True)
        dataframe.drop('average occurence in genome', axis=1, inplace=True)
        dataframe.drop('avg neigh size', axis=1, inplace=True)
        dataframe.drop('domain in whole database', axis=1, inplace=True)
        if self.user_correction == 'none':
            dataframe.drop('Poisson_correct', axis=1, inplace=True)
        else:
            dataframe.drop('Poisson', axis=1, inplace=True)

        return dataframe

    def poisson_distribution(self, p, n, domain, percentage_counter):

        if int(percentage_counter[domain]) < p * n:
            return 1
        else:
            sums = 0
            k = 0
            number_of_samples = int(percentage_counter[domain])
            while k < number_of_samples:
                sums += poisson.pmf(k=k, mu=n * p)
                k += 1
            if 1 - sums >= 0.0:
                return 1 - sums
            else:
                return 0.0
            # return 1 - sums

            # return 1 - sums
            # stat = poisson.pmf(k=int(percentage_counter[domain]), mu=n * p)
            # return stat

    def go(self):
        print(self.user_organisms)
        """Zipping all functions and execute them"""

        domain_information = self.collect_pfam_domain_description(
            "modules/data/domain_information")
        genome_id_size_in_gene = self.open_multiple_line(
            "modules/data/GENOME_ID_SIZE_IN_GENE_NEW.txt")
        with open('modules/data/counter_domains_new', 'r') as handler:
            pfam_counter_domains = json.load(handler)
        ################################################################################################################
        print("checkpoint 1")

        genome_id = [x[0] for x in genome_id_size_in_gene]
        size_in_gene = [x[1] for x in genome_id_size_in_gene]

        genome_number_overall = 0
        genome_number_to_stat = 0

        message_for_output = "You have looked for conserved neighbourhood for " + self.user_pfam + \
                             " domain, in range " + str(self.user_distance) + \
                             " bp, in " + self.user_organisms + " database."
        message_for_additional_data = "Pfam domain ,occurence in neighbourhoods, " \
                                      "average occurence in neighbourhood ,occurence genomes, " \
                                      "average occurence in genome,In what percentage?, " \
                                      "average neighborhood size, domain in whole database, " \
                                      "Poisson, Poisson_correct P, N, K, Family, Summary"
        print("checkpoint 2")

        ################################################################################################################
        # # LOCAL
        # file_names = self.open_database(self.user_organisms)
        ################################################################################################################
        ################################################################################################################
        # KLASTER
        # file_names = self.open_database('/home/klaster/ProFaNA v1.0/important_files/database', self.user_organism)
        ################################################################################################################
        ################################################################################################################
        # VULTR
        file_names = self.open_database(self.user_organisms)
        ################################################################################################################

        just_neigh_data = []
        genomes_with_domain = []
        neigh_genome_size = []
        counter = 0
        counter_db = 0
        sum_of_genes_in_neighborhoods = 0
        number_of_neighborhoods = 0
        percentage_counter = Counter()
        print("DLOK TIME START", datetime.datetime.now())

        for file in file_names:
            " OPEN FILE "
            # print(file)
            counter_db += 1
            counter += 1
            try:
                genes, start_coord, end_coord, orientation, domains, contig = self.create_six_list(file)
                " DZIELIMY PLIK NA KAWALKI"

            except FileNotFoundError:
                continue
            except ValueError:
                continue
            file_name_raw = file.split('/')
            tax_name = file_name_raw[-1]
            genome_number_overall += 1
            number_of_genes_in_genome = self.genome_size_in_gene(tax_name, genome_id, size_in_gene)
            coords = self.searching_for_domain_in_genome(self.user_pfam, start_coord, end_coord, orientation,
                                                         domains, contig, genes)
            if len(coords) > 0:
                genome_number_to_stat += 1

            for point in coords:
                number_of_neighborhoods += 1
                if self.user_distance != 0:

                    last_coordinate, first_coordinate, pfam_beg, pfam_end, searched_pfam_orientation = \
                        self.get_range_coordinates(point, self.user_distance)

                    if self.user_strand == 'both':
                        pfam_index_to_neigh = self.get_domains_both_strands(start_coord, end_coord, last_coordinate,
                                                                            first_coordinate, pfam_beg, pfam_end,
                                                                            contig,
                                                                            point)
                    else:
                        pfam_index_to_neigh_same, pfam_index_to_neigh_oposite = self.get_domains_plus_minus_strand(
                            start_coord, end_coord, last_coordinate,
                            first_coordinate, pfam_beg, pfam_end,
                            contig, point, orientation)
                        if self.user_strand == 'same':
                            pfam_index_to_neigh = pfam_index_to_neigh_same
                        elif self.user_strand == 'opposite':
                            pfam_index_to_neigh = pfam_index_to_neigh_oposite
                else:
                    pfam_index_to_neigh = self.get_domain_both_zeros(genes, point[4])

                genes_in_neigh = self.list_of_genes_in_neigh(genes, pfam_index_to_neigh)
                number_of_genes_in_neigh = self.size_in_genes(genes_in_neigh)
                whole = self.get_list_of_domain_in_neigh(pfam_index_to_neigh, file, domains)
                percentage_counter += Counter(set(whole))
                just_neigh_data.append(whole)
                genomes_with_domain.append(file)
                neigh_genome_size.append((number_of_genes_in_neigh, number_of_genes_in_genome))
                sum_of_genes_in_neighborhoods += int(number_of_genes_in_neigh)

        print("DLOK TIME OVER", datetime.datetime.now())
        print("concatenate_domain")

        alls = []
        neigh_counter = Counter()

        for i in just_neigh_data:
            for j in i:
                if j not in alls:
                    alls.append(j)
        genome_counter = Counter()
        print("MATRIX CREATING", datetime.datetime.now())

        counter = 0
        for domain_list, genome_in_domain, ng_size in zip(just_neigh_data, genomes_with_domain, neigh_genome_size):
            counter += 1
            genes, start_coord, end_coord, orientation, domains, contig = self.create_six_list(genome_in_domain)
            genome_domains_counter = Counter(domains)
            genome_counter += genome_domains_counter
            neigh_domains_counter = Counter(domain_list)
            neigh_counter += neigh_domains_counter

        print("MATRIX OVER", datetime.datetime.now())
        print("Wilcoxon start", datetime.datetime.now())

        najlepsze_dane = pd.DataFrame(columns=['occurence in neighborhoods',
                                               'average occurence in neighborhood', 'occurence genomes',
                                               'average occurence in genome',
                                               'In what percentage', 'Poisson', 'Poisson_correct', 'P', 'N', 'K',
                                               'avg neigh size', 'domain in whole database',
                                               'Family', 'Summary', ])

        scores_poisson = []
        counter = 0
        for i in alls:
            counter += 1

            p_KP = ((sum_of_genes_in_neighborhoods / genome_number_to_stat) * pfam_counter_domains[i]) / 326628067

            najlepsze_dane.at[i, 'P'] = p_KP
            najlepsze_dane.at[i, 'avg neigh size'] = sum_of_genes_in_neighborhoods / genome_number_to_stat
            najlepsze_dane.at[i, 'domain in whole database'] = pfam_counter_domains[i]
            # najlepsze_dane.at[i,'sum of database len'] = 326628067
            poisson_test_result = self.poisson_distribution(p_KP, number_of_neighborhoods, i, percentage_counter)
            najlepsze_dane.at[i, 'K'] = int(percentage_counter[i])
            najlepsze_dane.at[i, 'N'] = number_of_neighborhoods
            najlepsze_dane.at[i, 'P'] = p_KP
            najlepsze_dane.at[i, 'Poisson'] = poisson_test_result
            scores_poisson.append(poisson_test_result)

        # scores_after_correction = self.multiple_test_correction_list(scores,self.user_correction)
        poisson_after_correction = self.multiple_test_correction_list(scores_poisson, self.user_correction)
        print("Wilcoxon OVER", datetime.datetime.now())
        print("Collecting data START", datetime.datetime.now())

        for domena, poisson_corr in zip(alls, poisson_after_correction):
            najlepsze_dane.at[domena, 'Poisson_correct'] = poisson_corr
            najlepsze_dane.at[domena, 'In what percentage'] = percentage_counter[domena] / len(just_neigh_data)
            najlepsze_dane.at[domena, 'occurence in neighborhoods'] = neigh_counter[domena]
            najlepsze_dane.at[domena, 'average occurence in neighborhood'] = neigh_counter[domena] / len(
                just_neigh_data)
            najlepsze_dane.at[domena, 'occurence genomes'] = genome_counter[domena]
            najlepsze_dane.at[domena, 'average occurence in genome'] = genome_counter[domena] / len(
                set(genomes_with_domain))

        print("Collecting data OVER", datetime.datetime.now())
        print("Customizing_data START", datetime.datetime.now())
        after_cutoff = self.cutoff_value(najlepsze_dane, self.user_cutoff)
        sorted_table = self.sort_table(after_cutoff)
        added_information = self.add_information(sorted_table, domain_information)
        final_data = self.pfam_for_pf(added_information)

        print("Customizing_data OVER", datetime.datetime.now())
        message_down = " In selected database there is {} genomes and in {} searched domain " \
                       "was found".format(str(counter_db), str(genome_number_to_stat))
        self.save_data(self.user_output, message_for_output, message_for_additional_data, message_down, final_data)


class SuperSpeedAnalysisFromDomainAll:

    def __init__(self, user_pfam, user_distance, user_organism, user_cutoff, user_correction,
                 user_strand, user_output, user_level, skip_negative):
        """ Initializing input data
            Inputs:
                    user_pfam:  str (from pfam00000 to pfam99999)
                    user_distance: str non negative integer 1-20000
                    user_level: str
                    user_cutoff: str (none, 0, 0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.005,
                                 0.01, 0.05, 0.1, 0.2, 0.5)
                    user_correction: str (none, bonferroni, fdr_bh)
                    user_strand: str (both)
                    user_output: str
         """
        self.user_pfam = user_pfam  # pfam02696
        self.user_distance = user_distance  # 5000
        self.user_organism = user_organism  #
        self.user_cutoff = user_cutoff  # none
        self.user_correction = user_correction  # none
        self.user_strand = user_strand  # both
        self.user_output = user_output  # place
        self.user_level = user_level  # order
        self.skip_negative = skip_negative

    def open_singleline(self, path_to_file):
        """ Opens file that contain 1 column and strip it by space.
            Inputs:
                    path_to_file: str
        """

        with open(path_to_file, 'r') as f:
            file_names = [line.strip() for line in f]
        return file_names

    def open_singleline_num(self, path_to_file):
        """ Opens file that contain 1 column and strip it by space. """

        with open(path_to_file, 'r') as f:
            file_names = [float(line.strip()) for line in f]
        return file_names

    def open_multiple_line(self, path_to_file):
        """ Opens file that contains data about single genome """

        plik = []
        with open(path_to_file) as inputfile:
            for line in inputfile:
                plik.append(line.strip().split())
        return plik

    def save_data(self, file_name, message_for_output, message_for_additional_data, message_down, complete_data):
        """ Saves all stuff together as one file """

        # # LOCAL
        # with open('/mnt/d/45.76.38.24/final_site-project' + file_name + '.txt', 'w') as output_file:
        # VULTR
        # with open('/home/djangoadmin/final_site-project/media/results/test/' + file_name, 'w') as output_file:
        #     output_file.write(message_for_output)
        #     output_file.write("\n")
        #     output_file.write(message_down)
        #     output_file.write("\n")
        #     output_file.write(message_for_additional_data)
        #     output_file.write("\n")
        #     # output_file.write(complete_data)
        #     for row in complete_data.iterrows():
        #         index, data = row
        #         pre = data.tolist()
        #         output_file.write(index)
        #         output_file.write("\t")
        #         output_file.write("\t".join([str(i) for i in pre]))
        #         output_file.write("\n")

        complete_data.index.name = 'Domain'
        complete_data.to_csv('/home/djangoadmin/final_site-project' + file_name, sep='\t', mode='w')

    def create_six_list(self, file):
        """ Open single genome file chew it and return 6 lists -> GENE, START_COORD,
        END_COORDS ,ORIENTATION, DOMAINS """

        # KLASTER
        # data = self.open_multiple_line("/home/klaster/ProFaNA v1.0/nowa_baza_danych/"+file)

        # LOCAL
        # data = self.open_multiple_line("/mnt/d/45.76.38.24/final_site-project/important_files/database_file/"+file)

        # VULTR
        data = self.open_multiple_line("/home/djangoadmin/final_site-project/important_files/nowa_baza_danych/" + file)
        start_coord = []
        end_coord = []
        orientation = []
        domains = []
        genes = []
        contig = []
        for bit in data:
            genes.append(int(bit[0]))
            start_coord.append(int(bit[1]))
            end_coord.append(int(bit[2]))
            orientation.append(bit[3])
            domains.append(bit[4])
            contig.append(bit[5])
        return genes, start_coord, end_coord, orientation, domains, contig

    def open_database(self, level, organism):

        print(self.user_pfam, self.user_distance, self.user_organism, self.user_cutoff, self.user_correction,
              self.user_strand, self.user_output, self.user_level, self.skip_negative)
        tax = organism.replace("_", " ")
        # CHECK LEVEL AND GET LIST OF ORGANISMS
        if level == "genus":
            # VULTR
            with open('/home/djangoadmin/final_site-project/important_files/genus_level.json', "r") as handler:
                # LOCAL
                # with open('/mnt/d/45.76.38.24/final_site-project/important_files/genus_level.json', "r") as handler:
                taxonomy = json.load(handler)
        elif level == "family":
            # VULTR
            with open('/home/djangoadmin/final_site-project/important_files/family_level.json', "r") as handler:
                # LOCAL
                # with open('/mnt/d/45.76.38.24/final_site-project/important_files/family_level.json', "r") as handler:
                taxonomy = json.load(handler)
        elif level == "order":
            # VULTR
            with open('/home/djangoadmin/final_site-project/important_files/order_level.json', "r") as handler:
                # LOCAL
                # with open('/mnt/d/45.76.38.24/final_site-project/important_files/order_level.json', "r") as handler:
                taxonomy = json.load(handler)
        elif level == "class":
            # VULTR
            with open('/home/djangoadmin/final_site-project/important_files/class_level.json', "r") as handler:
                # LOCAL
                # with open('/mnt/d/45.76.38.24/final_site-project/important_files/class_level.json', "r") as handler:
                taxonomy = json.load(handler)
        elif level == "phylum":
            # VULTR
            with open('/home/djangoadmin/final_site-project/important_files/phylum_level.json', "r") as handler:
                # LOCAL
                # with open('/mnt/d/45.76.38.24/final_site-project/important_files/phylum_level.json', "r") as handler:
                taxonomy = json.load(handler)
        # elif level =='all':

        # organism_list = taxonomy[tax]

        return taxonomy

    def genome_size_in_gene(self, genome_id, list_of_genome, list_of_size):
        """ Takes genome's id, list of genome and lists with data about how much genes in genome"""

        index_genome = list_of_genome.index(genome_id)
        genome_size_in_gene = list_of_size[index_genome]
        return genome_size_in_gene

    def size_in_genes(self, genome_or_genes):
        """ Takes list of genes and returns size of neighbourhood """

        list_of_genes = []
        for gene in genome_or_genes:
            list_of_genes.append(gene)
        return len(set(list_of_genes))

    def list_of_genes_in_neigh(self, list_of_genes, index_list):
        """ Takes overall gene's list in genome and creates complete list of genes in neighbourhood """

        list_of_genes_in_neigh = []
        for pfam_domain in index_list:
            list_of_genes_in_neigh.append(list_of_genes[pfam_domain])
        return list_of_genes_in_neigh

    def searching_for_domain_in_genome(self, pfam, start_coord, end_coord, orient, domains, contig, genes):
        """
        Takes user's pfam domain searches through lists contains coordinates,
        orientation, domains, contigs and returns complete data about all
        users' domains in genome
        """
        coords = []
        coords_counter = 0
        for domain in domains:
            one_coords_part = []
            if domain == pfam:
                orientation_pfam_domain = orient[coords_counter]
                one_coords_part.append(start_coord[coords_counter])
                one_coords_part.append(end_coord[coords_counter])
                one_coords_part.append(orientation_pfam_domain)
                one_coords_part.append(contig[coords_counter])
                one_coords_part.append(genes[coords_counter])
                coords.append(one_coords_part)
                coords_counter += 1
            else:
                coords_counter += 1
                continue
        return coords

    def presence_confirmation(self, coords, file, pfam):
        """
                Just prints out information how many pfam domains is in each genome
                during analysis propably i will have to remove it before putting it to website ...
                """

        if len(coords) == 0:
            print("In genome  " + file + " pfam domain you have been searching do not exist")
        elif len(coords) == 1:
            print("In genome  " + file + " there is " + str(len(coords)) + " " + pfam + " domain")
        else:
            print("In genome  " + file + " there are " + str(len(coords)) + " " + pfam + " domains")

    def get_range_coordinates(self, target_pfam, distance):
        """
        Based on how large neighbourhood user wants to analyze creates
        points FROM and TO, additionaly shows where main user's pfam begins and ENDS
        with orientation on strands
        """

        last_coordinate = target_pfam[1] + int(distance)
        first_coordinate = target_pfam[0] - int(distance)
        pfam_beg = target_pfam[0]
        pfam_end = target_pfam[1]
        searched_pfam_orientation = target_pfam[2]
        return last_coordinate, first_coordinate, pfam_beg, pfam_end, searched_pfam_orientation

    def get_domain_both_zeros(self, ref_gene, q_gene):

        index_list_gene = []
        counter = 0
        for gene in ref_gene:
            if gene == int(q_gene):
                index_list_gene.append(counter)
            counter += 1
        return index_list_gene

    def get_domains_both_strands(self, start_coord, end_coord, last_coordinate, first_coordinate, pfam_beg, pfam_end,
                                 contig, point):

        pfam_index_to_neigh = []
        s_coord_counter = 0
        for s_coord in start_coord:
            if s_coord <= last_coordinate and s_coord >= pfam_beg and contig[s_coord_counter] == point[
                3] and s_coord_counter not in pfam_index_to_neigh:
                pfam_index_to_neigh.append(s_coord_counter)
                s_coord_counter += 1
            else:
                s_coord_counter += 1
                continue
        e_coord_counter = 0
        for e_coord in end_coord:
            if e_coord >= first_coordinate and e_coord <= pfam_end and contig[e_coord_counter] == point[
                3] and e_coord_counter not in pfam_index_to_neigh:
                pfam_index_to_neigh.append(e_coord_counter)
                e_coord_counter += 1
            else:
                e_coord_counter += 1
                continue
        return pfam_index_to_neigh

    def collect_pfam_domain_description(self, file):
        """ Opens file only to gather information"""

        data = {}
        with open(file) as inputfile:
            for line in inputfile:
                try:

                    one_line = line.strip().split("@")
                    domena = one_line[0]
                    pre_family = one_line[1][8:]
                    delete_this = "(" + domena.upper() + ")"
                    family = pre_family.replace(delete_this, "")
                    summary = one_line[2][9:]
                    data[domena] = (family, summary)

                except IndexError:
                    continue
        return data

    def get_domains_plus_minus_strand(self, start_coord, end_coord, last_coordinate, first_coordinate,
                                      pfam_beg, pfam_end, contig, point, orientation):

        pfam_index_to_neigh_same_strand = []
        pfam_index_to_neigh_oposite_strand = []
        s_coord_counter = 0
        for s_coord in start_coord:
            if s_coord <= last_coordinate and s_coord >= pfam_beg and orientation[s_coord_counter] == point[2] and \
                    contig[s_coord_counter] == point[3] and s_coord_counter not in pfam_index_to_neigh_same_strand:
                pfam_index_to_neigh_same_strand.append(s_coord_counter)
                s_coord_counter += 1
            elif s_coord <= last_coordinate and s_coord >= pfam_beg and orientation[s_coord_counter] != point[2] and \
                    contig[s_coord_counter] == point[
                3] and s_coord_counter not in pfam_index_to_neigh_oposite_strand:
                pfam_index_to_neigh_oposite_strand.append(s_coord_counter)
                s_coord_counter += 1
            else:
                s_coord_counter += 1
                continue
        e_coord_counter = 0
        for e_coord in end_coord:
            if e_coord >= first_coordinate and e_coord <= pfam_end and orientation[e_coord_counter] == point[2] and \
                    contig[e_coord_counter] == point[3] and e_coord_counter not in pfam_index_to_neigh_same_strand:
                pfam_index_to_neigh_same_strand.append(e_coord_counter)
                e_coord_counter += 1
            elif e_coord >= first_coordinate and e_coord <= pfam_end and orientation[e_coord_counter] != point[
                2] and contig[e_coord_counter] == point[
                3] and e_coord_counter not in pfam_index_to_neigh_oposite_strand:
                pfam_index_to_neigh_oposite_strand.append(e_coord_counter)
                e_coord_counter += 1
            else:
                e_coord_counter += 1
                continue
        return pfam_index_to_neigh_same_strand, pfam_index_to_neigh_oposite_strand

    def get_list_of_domain_in_neigh(self, pfam_index, file, domains):

        to_counter = []
        party = []
        party.append(file)
        for part in pfam_index:
            party.append(domains[part])
            to_counter.append(domains[part])
        #        neighbourhood_complete.append(party)

        return to_counter

    def counting_dlok_or_dglob(self, some_counter, genome_neigh_size, mighty_domains):
        dlok_glob = []
        for domain_mighty in mighty_domains:
            if domain_mighty in some_counter.keys():
                lok_glob = some_counter.get(domain_mighty) / int(genome_neigh_size)
                dlok_glob.append(lok_glob)
            else:
                dlok_glob.append(0)
        return dlok_glob

    def count_for_dlok(self, somecounter, neigh_size, tax_name, dlok_dataframe):
        for k, v in somecounter.items():
            dlok_dataframe.at[tax_name, k] = v / neigh_size

    def multiple_test_correction(self, some_dataframe, correction_met):

        if correction_met == 'none':
            return some_dataframe
        else:
            value_to_correct = [float(x) for x in some_dataframe.PVALUE.tolist()]
            reject, pvals_corrected, alphaSidak, alphaBonf = correction(pvals=value_to_correct,
                                                                        method=correction_met,
                                                                        is_sorted=False, returnsorted=False)
            pvals = pvals_corrected.tolist()
            some_dataframe.PVALUE = pvals
            return some_dataframe

    def multiple_test_correction_list(self, some_tuple_list, correction_met):
        pvalues = [float(x) for x in some_tuple_list]

        if correction_met == 'none':
            output = pvalues

        else:
            try:

                reject, pvals_corrected, alphaSidak, alphaBonf = correction(pvals=pvalues,
                                                                            method=correction_met,
                                                                            is_sorted=False, returnsorted=False)
                pvals_after_corr = pvals_corrected.tolist()
                output = pvals_after_corr
            except ZeroDivisionError:
                output = pvalues
        return output

    def cutoff_value(self, some_dataframe, cutoff):

        domain_list = list(some_dataframe.index)
        if cutoff == 'none':
            return some_dataframe
        elif cutoff == '0':
            for i in domain_list:

                diff = float(some_dataframe.loc[i, 'PVALUE'])
                if diff < 0 and diff > 0:
                    some_dataframe = some_dataframe.drop([i])
            return some_dataframe
        else:
            cutoff = float(cutoff)
            for i in domain_list:
                diff = float(some_dataframe.loc[i, 'PVALUE'])
                if diff > cutoff:
                    some_dataframe = some_dataframe.drop([i])
            return some_dataframe

    def sort_table(self, some_dataframe):

        sorted_data = some_dataframe.sort_values('Poisson', ascending=True)
        return sorted_data

    def add_information(self, some_dataframe, dictionary):

        indeksy = list(some_dataframe.index)
        for i in indeksy:
            try:

                pfam = i[0:2] + i[4:]
                family = dictionary[pfam][0]
                summary = dictionary[pfam][1]

                some_dataframe.at[i, 'Family'] = family
                some_dataframe.at[i, 'Summary'] = summary
            except KeyError:
                continue
        return some_dataframe

    def poisson_distribution_mine(self, k, p, n):
        """
        Return value (float) of poisson distribution mine concept.
            Parameters:
                k (int): number of domain in neighborhoods
                p (float): (sum (average) of neighborhoods lenghts * sum of number of domain) / (average) sum of genome sizes
                n (int): number of neighborhoods
            How to agreggate data:
                k - Counter of all domains should be ok and easy to do
                p - sum of neighborhoods (easy += ) , sum of number of domain (z Countera), sum of genome sizes ( easy +=)
                n - number of neighborhoods (easy +=)
        """
        first_value = (n * p) ** k
        second_value = math.e ** (-1 * n * p)
        third_value = math.factorial(k)
        return (first_value * second_value) / third_value

    def poisson_distribution_KP(self, k, p, n):
        """
        Return value of poisson distribution mine concept
            Parameters:
                k (int): sum of domain occurrence in all neighborhoods
                n (int): number of neighborhoods
                p (float): ( average lengh of neighborhood * number of domain in all database) / sum of lengh of all database
            How to aggregate data:
                k Counter of all domains
                n - number of neigghborhoods (easy +=)
                p - average lengh of neighborhood ( easy +=) / number of genomes, number of domain in all ( from file), sum of lenght of all database ( const value)
        """
        first_value = (n * p) ** k
        second_value = math.e ** (-1 * n * p)
        third_value = math.factorial(k)
        return (first_value * second_value) / third_value

    def pfam_for_pf(self, dataframe):
        """Final customization of dataframe
            --> PF02696 instead of pfam02696
            --> PVALUE in 10e-3 format
            --> In what percentage % format
            --> avg occurences and density  have now 3 digits after coma
            --> drop NOam
        """

        indexy = dataframe.index
        for i in indexy:
            dataframe = dataframe.rename(index={i: i[:2].upper() + i[4:]})
        dataframe['Poisson'] = dataframe['Poisson'].map('{:.3e}'.format)
        dataframe['Poisson_correct'] = dataframe['Poisson_correct'].map('{:.3e}'.format)

        dataframe['In what percentage'] = dataframe['In what percentage'].map('{:.3%}'.format)
        dataframe['average occurence in neighborhood'] = dataframe['average occurence in neighborhood'].map(
            '{:.3}'.format)
        dataframe['average occurence in genome'] = dataframe['average occurence in genome'].map('{:.3}'.format)

        # dataframe = dataframe.round({"average occurence in neighborhood": 3,"average occurence in genome": 3,
        #                              "Density difference": 3})
        if "NOam" in indexy:
            dataframe = dataframe.drop(index="NOam")
        # dataframe = dataframe.style.format({'In what percentage': '{:,.3%}'.format})
        dataframe.drop('average occurence in neighborhood', axis=1, inplace=True)
        dataframe.drop('average occurence in genome', axis=1, inplace=True)
        dataframe.drop('avg neigh size', axis=1, inplace=True)
        dataframe.drop('domain in whole database', axis=1, inplace=True)
        if self.user_correction == 'none':
            dataframe.drop('Poisson_correct', axis=1, inplace=True)
        else:
            dataframe.drop('Poisson', axis=1, inplace=True)

        return dataframe

    def poisson_distribution(self, p, n, domain, percentage_counter):

        if int(percentage_counter[domain]) < p * n:
            return 1
        else:
            sums = 0
            k = 0
            number_of_samples = int(percentage_counter[domain])
            while k < number_of_samples:
                sums += poisson.pmf(k=k, mu=n * p)
                k += 1
            if sums >= 1:
                return float(0.0)
            else:
                return 1 - sums

    def go(self):
        print(self.user_organism)
        """Zipping all functions and execute them"""
        ################################################################################################################
        # JUST LOCALLY
        # print("checkpoint #1")
        # genome_id_size_in_gene = self.open_multiple_line(
        #     '/mnt/d/45.76.38.24/final_site-project/important_files/GENOME_ID_SIZE_IN_GENE.txt')
        # domain_information = self.collect_pfam_domain_description(
        #     '/mnt/d/45.76.38.24/final_site-project/important_files/domain_information')
        ################################################################################################################

        ################################################################################################################
        # # KLASTER
        #
        # print("checkpoint #1")
        # # KLASTER
        # domain_information = self.collect_pfam_domain_description(
        #     '/home/klaster/ProFaNA v1.0/important_files/domain_information')
        # genome_id_size_in_gene = self.open_multiple_line(
        #     '/home/klaster/ProFaNA v1.0/important_files/GENOME_ID_SIZE_IN_GENE.txt')
        ################################################################################################################

        ################################################################################################################
        # # VULTR
        domain_information = self.collect_pfam_domain_description(
            "/home/djangoadmin/final_site-project/important_files/domain_information")
        genome_id_size_in_gene = self.open_multiple_line(
            "/home/djangoadmin/final_site-project/important_files/GENOME_ID_SIZE_IN_GENE.txt")
        with open('/home/djangoadmin/final_site-project/important_files/counter_domains', 'r') as handler:
            pfam_counter_domains = json.load(handler)
        ################################################################################################################
        print("checkpoint 1")

        genome_id = [x[0] for x in genome_id_size_in_gene]
        size_in_gene = [x[1] for x in genome_id_size_in_gene]

        genome_number_overall = 0
        genome_number_to_stat = 0

        message_for_output = "You have looked for conserved neighbourhood for " + self.user_pfam + \
                             " domain, in range " + str(self.user_distance) + \
                             " bp, in " + self.user_organism + " database."
        message_for_additional_data = "Pfam domain , PVALUE,  occurence in neighbourhoods, " \
                                      "average occurence in neighbourhood ,occurence genomes, " \
                                      "average occurence in genome, avg DLOK-DGLOB, In what percentage?, " \
                                      "Family, Summary"
        print("checkpoint 2")

        ################################################################################################################
        # # LOCAL
        # file_names = self.open_database(self.user_level,self.user_organism)
        ################################################################################################################
        ################################################################################################################
        # KLASTER
        # file_names = self.open_database('/home/klaster/ProFaNA v1.0/important_files/database', self.user_organism)
        ################################################################################################################
        ################################################################################################################
        # VULTR
        level = self.open_database(self.user_level, self.user_organism)
        ################################################################################################################

        for key in level:
            print(key)

            # Create important variables
            just_neigh_data = []
            genomes_with_domain = []
            neigh_genome_size = []
            counter = 0
            counter_db = 0
            sum_of_genes_in_neighborhoods = 0
            number_of_neighborhoods = 0
            percentage_counter = Counter()

            # Dlok operations
            print("DLOK TIME START", datetime.datetime.now())
            for file in level[key]:
                # print(file)
                counter_db += 1
                counter += 1
                try:
                    genes, start_coord, end_coord, orientation, domains, contig = self.create_six_list(file)
                except FileNotFoundError:
                    continue
                except ValueError:
                    continue
                file_name_raw = file.split('/')
                tax_name = file_name_raw[-1]
                genome_number_overall += 1
                number_of_genes_in_genome = self.genome_size_in_gene(tax_name, genome_id, size_in_gene)
                coords = self.searching_for_domain_in_genome(self.user_pfam, start_coord, end_coord, orientation,
                                                             domains, contig, genes)
                if len(coords) > 0:
                    genome_number_to_stat += 1

                for point in coords:
                    number_of_neighborhoods += 1

                    if self.user_distance != 0:
                        last_coordinate, first_coordinate, pfam_beg, pfam_end, searched_pfam_orientation = \
                            self.get_range_coordinates(point, self.user_distance)
                        pfam_index_to_neigh = self.get_domains_both_strands(start_coord, end_coord, last_coordinate,
                                                                            first_coordinate, pfam_beg, pfam_end,
                                                                            contig,
                                                                            point)
                    else:
                        pfam_index_to_neigh = self.get_domain_both_zeros(genes, point[4])

                    genes_in_neigh = self.list_of_genes_in_neigh(genes, pfam_index_to_neigh)
                    number_of_genes_in_neigh = self.size_in_genes(genes_in_neigh)
                    whole = self.get_list_of_domain_in_neigh(pfam_index_to_neigh, file, domains)
                    percentage_counter += Counter(set(whole))
                    just_neigh_data.append(whole)
                    genomes_with_domain.append(file)
                    neigh_genome_size.append((number_of_genes_in_neigh, number_of_genes_in_genome))
                    sum_of_genes_in_neighborhoods += int(number_of_genes_in_neigh)
                    # pdb.set_trace()
            print("DLOK TIME OVER", datetime.datetime.now())
            print("concatenate_domain")

            alls = []
            neigh_counter = Counter()

            for i in just_neigh_data:
                for j in i:
                    if j not in alls:
                        alls.append(j)
            genome_counter = Counter()
            print("MATRIX CREATING", datetime.datetime.now())

            counter = 0
            for domain_list, genome_in_domain, ng_size in zip(just_neigh_data, genomes_with_domain, neigh_genome_size):
                counter += 1
                genes, start_coord, end_coord, orientation, domains, contig = self.create_six_list(genome_in_domain)
                genome_domains_counter = Counter(domains)
                genome_counter += genome_domains_counter
                neigh_domains_counter = Counter(domain_list)
                neigh_counter += neigh_domains_counter

            print("MATRIX OVER", datetime.datetime.now())
            print("Wilcoxon start", datetime.datetime.now())

            najlepsze_dane = pd.DataFrame(columns=['occurence in neighborhoods',
                                                   'average occurence in neighborhood', 'occurence genomes',
                                                   'average occurence in genome',
                                                   'In what percentage', 'Poisson', 'Poisson_correct', 'P', 'N', 'K',
                                                   'avg neigh size', 'domain in whole database',
                                                   'Family', 'Summary', ])

            scores_poisson = []

            counter = 0
            for i in alls:
                counter += 1

                p_KP = ((sum_of_genes_in_neighborhoods / genome_number_to_stat) * pfam_counter_domains[i]) / 326628067
                najlepsze_dane.at[i, 'P'] = p_KP
                najlepsze_dane.at[i, 'avg neigh size'] = sum_of_genes_in_neighborhoods / genome_number_to_stat
                najlepsze_dane.at[i, 'domain in whole database'] = pfam_counter_domains[i]
                # najlepsze_dane.at[i,'sum of database len'] = 326628067
                poisson_test_result = self.poisson_distribution(p_KP, number_of_neighborhoods, i, percentage_counter)
                najlepsze_dane.at[i, 'K'] = int(percentage_counter[i])
                najlepsze_dane.at[i, 'N'] = number_of_neighborhoods
                najlepsze_dane.at[i, 'P'] = p_KP
                najlepsze_dane.at[i, 'Poisson'] = poisson_test_result
                scores_poisson.append(poisson_test_result)

            poisson_after_correction = self.multiple_test_correction_list(scores_poisson, self.user_correction)
            print("Wilcoxon OVER", datetime.datetime.now())
            print("Collecting data START", datetime.datetime.now())

            for domena, poisson_corr, in zip(alls, poisson_after_correction):
                najlepsze_dane.at[domena, 'Poisson_correct'] = poisson_corr
                najlepsze_dane.at[domena, 'In what percentage'] = percentage_counter[domena] / len(just_neigh_data)
                najlepsze_dane.at[domena, 'occurence in neighborhoods'] = neigh_counter[domena]
                najlepsze_dane.at[domena, 'average occurence in neighborhood'] = neigh_counter[domena] / len(
                    just_neigh_data)
                najlepsze_dane.at[domena, 'occurence genomes'] = genome_counter[domena]
                najlepsze_dane.at[domena, 'average occurence in genome'] = genome_counter[domena] / len(
                    set(genomes_with_domain))

            print("Collecting data OVER", datetime.datetime.now())
            print("Customizing_data START", datetime.datetime.now())
            after_cutoff = self.cutoff_value(najlepsze_dane, self.user_cutoff)
            sorted_table = self.sort_table(after_cutoff)
            added_information = self.add_information(sorted_table, domain_information)
            final_data = self.pfam_for_pf(added_information)

            print("Customizing_data OVER", datetime.datetime.now())
            message_down = " In selected database there is {} genomes and in {} searched domain " \
                           "was found".format(str(counter_db), str(genome_number_to_stat))
            self.save_data('{}_{}'.format(self.user_level, key), message_for_output, message_for_additional_data,
                           message_down, final_data)


class SuperSpeedAnalysisFromDomainAllAvg:

    def __init__(self, user_pfam, user_distance, user_organism, user_cutoff, user_correction,
                 user_strand, user_output, user_level):
        """ Initializing input data
            Inputs:
                    user_pfam:  str (from pfam00000 to pfam99999)
                    user_distance: str non negative integer 1-20000
                    user_level: str
                    user_cutoff: str (none, 0, 0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.005,
                                 0.01, 0.05, 0.1, 0.2, 0.5)
                    user_correction: str (none, bonferroni, fdr_bh)
                    user_strand: str (both)
                    user_output: str
         """
        self.user_pfam = user_pfam  # pfam02696
        self.user_distance = user_distance  # 5000
        self.user_organism = user_organism  # none
        self.user_cutoff = user_cutoff  # none
        self.user_correction = user_correction  # none
        self.user_strand = user_strand  # both
        self.user_output = user_output  # place
        self.user_level = user_level  # order

    def open_singleline(self, path_to_file):
        """ Opens file that contain 1 column and strip it by space.
            Inputs:
                    path_to_file: str
        """

        with open(path_to_file, 'r') as f:
            file_names = [line.strip() for line in f]
        return file_names

    def open_singleline_num(self, path_to_file):
        """ Opens file that contain 1 column and strip it by space. """

        with open(path_to_file, 'r') as f:
            file_names = [float(line.strip()) for line in f]
        return file_names

    def open_multiple_line(self, path_to_file):
        """ Opens file that contains data about single genome """

        plik = []
        with open(path_to_file) as inputfile:
            for line in inputfile:
                plik.append(line.strip().split())
        return plik

    def save_data(self, file_name, message_for_output, message_for_additional_data, message_down, complete_data):
        """ Saves all stuff together as one file """

        # # LOCAL
        # with open('/mnt/d/45.76.38.24/final_site-project' + file_name + '.txt', 'w') as output_file:
        # VULTR
        # with open('/home/djangoadmin/final_site-project/media/results/test/' + file_name, 'w') as output_file:
        #     output_file.write(message_for_output)
        #     output_file.write("\n")
        #     output_file.write(message_down)
        #     output_file.write("\n")
        #     output_file.write(message_for_additional_data)
        #     output_file.write("\n")
        #     # output_file.write(complete_data)
        #     for row in complete_data.iterrows():
        #         index, data = row
        #         pre = data.tolist()
        #         output_file.write(index)
        #         output_file.write("\t")
        #         output_file.write("\t".join([str(i) for i in pre]))
        #         output_file.write("\n")

        complete_data.index.name = 'Domain'
        print('/home/djangoadmin/final_site-project/media/results/temp/' + file_name)
        complete_data.to_csv('/home/djangoadmin/final_site-project/media/results/temp/' + file_name, sep='\t', mode='w')

    def create_six_list(self, file):
        """ Open single genome file chew it and return 6 lists -> GENE, START_COORD,
        END_COORDS ,ORIENTATION, DOMAINS """

        # KLASTER
        # data = self.open_multiple_line("/home/klaster/ProFaNA v1.0/nowa_baza_danych/"+file)

        # LOCAL
        # data = self.open_multiple_line("/mnt/d/45.76.38.24/final_site-project/important_files/database_file/"+file)

        # VULTR
        data = self.open_multiple_line("/home/djangoadmin/final_site-project/important_files/nowa_baza_danych/" + file)
        start_coord = []
        end_coord = []
        orientation = []
        domains = []
        genes = []
        contig = []
        for bit in data:
            genes.append(int(bit[0]))
            start_coord.append(int(bit[1]))
            end_coord.append(int(bit[2]))
            orientation.append(bit[3])
            domains.append(bit[4])
            contig.append(bit[5])
        return genes, start_coord, end_coord, orientation, domains, contig

    def open_database(self, level, organism):

        print(self.user_pfam, self.user_distance, self.user_organism, self.user_cutoff, self.user_correction,
              self.user_strand, self.user_output, self.user_level)
        tax = organism.replace("_", " ")
        # CHECK LEVEL AND GET LIST OF ORGANISMS
        if level == "genus":
            # VULTR
            with open('/home/djangoadmin/final_site-project/important_files/genus_level.json', "r") as handler:
                # LOCAL
                # with open('/mnt/d/45.76.38.24/final_site-project/important_files/genus_level.json', "r") as handler:
                taxonomy = json.load(handler)
        elif level == "family":
            # VULTR
            with open('/home/djangoadmin/final_site-project/important_files/family_level.json', "r") as handler:
                # LOCAL
                # with open('/mnt/d/45.76.38.24/final_site-project/important_files/family_level.json', "r") as handler:
                taxonomy = json.load(handler)
        elif level == "order":
            # VULTR
            with open('/home/djangoadmin/final_site-project/important_files/order_level.json', "r") as handler:
                # LOCAL
                # with open('/mnt/d/45.76.38.24/final_site-project/important_files/order_level.json', "r") as handler:
                taxonomy = json.load(handler)
        elif level == "class":
            # VULTR
            with open('/home/djangoadmin/final_site-project/important_files/class_level.json', "r") as handler:
                # LOCAL
                # with open('/mnt/d/45.76.38.24/final_site-project/important_files/class_level.json', "r") as handler:
                taxonomy = json.load(handler)
        elif level == "phylum":
            # VULTR
            with open('/home/djangoadmin/final_site-project/important_files/phylum_level.json', "r") as handler:
                # LOCAL
                # with open('/mnt/d/45.76.38.24/final_site-project/important_files/phylum_level.json', "r") as handler:
                taxonomy = json.load(handler)

        # organism_list = taxonomy[tax]

        return taxonomy

    def genome_size_in_gene(self, genome_id, list_of_genome, list_of_size):
        """ Takes genome's id, list of genome and lists with data about how much genes in genome"""

        index_genome = list_of_genome.index(genome_id)
        genome_size_in_gene = list_of_size[index_genome]
        return genome_size_in_gene

    def size_in_genes(self, genome_or_genes):
        """ Takes list of genes and returns size of neighbourhood """

        list_of_genes = []
        for gene in genome_or_genes:
            list_of_genes.append(gene)
        return len(set(list_of_genes))

    def list_of_genes_in_neigh(self, list_of_genes, index_list):
        """ Takes overall gene's list in genome and creates complete list of genes in neighbourhood """

        list_of_genes_in_neigh = []
        for pfam_domain in index_list:
            list_of_genes_in_neigh.append(list_of_genes[pfam_domain])
        return list_of_genes_in_neigh

    def searching_for_domain_in_genome(self, pfam, start_coord, end_coord, orient, domains, contig, genes):
        """
        Takes user's pfam domain searches through lists contains coordinates,
        orientation, domains, contigs and returns complete data about all
        users' domains in genome
        """
        coords = []
        coords_counter = 0
        for domain in domains:
            one_coords_part = []
            if domain == pfam:
                orientation_pfam_domain = orient[coords_counter]
                one_coords_part.append(start_coord[coords_counter])
                one_coords_part.append(end_coord[coords_counter])
                one_coords_part.append(orientation_pfam_domain)
                one_coords_part.append(contig[coords_counter])
                one_coords_part.append(genes[coords_counter])
                coords.append(one_coords_part)
                coords_counter += 1
            else:
                coords_counter += 1
                continue
        return coords

    def presence_confirmation(self, coords, file, pfam):
        """
                Just prints out information how many pfam domains is in each genome
                during analysis propably i will have to remove it before putting it to website ...
                """

        if len(coords) == 0:
            print("In genome  " + file + " pfam domain you have been searching do not exist")
        elif len(coords) == 1:
            print("In genome  " + file + " there is " + str(len(coords)) + " " + pfam + " domain")
        else:
            print("In genome  " + file + " there are " + str(len(coords)) + " " + pfam + " domains")

    def get_range_coordinates(self, target_pfam, distance):
        """
        Based on how large neighbourhood user wants to analyze creates
        points FROM and TO, additionaly shows where main user's pfam begins and ENDS
        with orientation on strands
        """

        last_coordinate = target_pfam[1] + int(distance)
        first_coordinate = target_pfam[0] - int(distance)
        pfam_beg = target_pfam[0]
        pfam_end = target_pfam[1]
        searched_pfam_orientation = target_pfam[2]
        return last_coordinate, first_coordinate, pfam_beg, pfam_end, searched_pfam_orientation

    def get_domain_both_zeros(self, ref_gene, q_gene):

        index_list_gene = []
        counter = 0
        for gene in ref_gene:
            if gene == int(q_gene):
                index_list_gene.append(counter)
            counter += 1
        return index_list_gene

    def get_domains_both_strands(self, start_coord, end_coord, last_coordinate, first_coordinate, pfam_beg, pfam_end,
                                 contig, point):

        pfam_index_to_neigh = []
        s_coord_counter = 0
        for s_coord in start_coord:
            if s_coord <= last_coordinate and s_coord >= pfam_beg and contig[s_coord_counter] == point[
                3] and s_coord_counter not in pfam_index_to_neigh:
                pfam_index_to_neigh.append(s_coord_counter)
                s_coord_counter += 1
            else:
                s_coord_counter += 1
                continue
        e_coord_counter = 0
        for e_coord in end_coord:
            if e_coord >= first_coordinate and e_coord <= pfam_end and contig[e_coord_counter] == point[
                3] and e_coord_counter not in pfam_index_to_neigh:
                pfam_index_to_neigh.append(e_coord_counter)
                e_coord_counter += 1
            else:
                e_coord_counter += 1
                continue
        return pfam_index_to_neigh

    def collect_pfam_domain_description(self, file):
        """ Opens file only to gather information"""

        data = {}
        with open(file) as inputfile:
            for line in inputfile:
                try:

                    one_line = line.strip().split("@")
                    domena = one_line[0]
                    pre_family = one_line[1][8:]
                    delete_this = "(" + domena.upper() + ")"
                    family = pre_family.replace(delete_this, "")
                    summary = one_line[2][9:]
                    data[domena] = (family, summary)

                except IndexError:
                    continue
        return data

    def get_domains_plus_minus_strand(self, start_coord, end_coord, last_coordinate, first_coordinate,
                                      pfam_beg, pfam_end, contig, point, orientation):

        pfam_index_to_neigh_same_strand = []
        pfam_index_to_neigh_oposite_strand = []
        s_coord_counter = 0
        for s_coord in start_coord:
            if s_coord <= last_coordinate and s_coord >= pfam_beg and orientation[s_coord_counter] == point[2] and \
                    contig[s_coord_counter] == point[3] and s_coord_counter not in pfam_index_to_neigh_same_strand:
                pfam_index_to_neigh_same_strand.append(s_coord_counter)
                s_coord_counter += 1
            elif s_coord <= last_coordinate and s_coord >= pfam_beg and orientation[s_coord_counter] != point[2] and \
                    contig[s_coord_counter] == point[
                3] and s_coord_counter not in pfam_index_to_neigh_oposite_strand:
                pfam_index_to_neigh_oposite_strand.append(s_coord_counter)
                s_coord_counter += 1
            else:
                s_coord_counter += 1
                continue
        e_coord_counter = 0
        for e_coord in end_coord:
            if e_coord >= first_coordinate and e_coord <= pfam_end and orientation[e_coord_counter] == point[2] and \
                    contig[e_coord_counter] == point[3] and e_coord_counter not in pfam_index_to_neigh_same_strand:
                pfam_index_to_neigh_same_strand.append(e_coord_counter)
                e_coord_counter += 1
            elif e_coord >= first_coordinate and e_coord <= pfam_end and orientation[e_coord_counter] != point[
                2] and contig[e_coord_counter] == point[
                3] and e_coord_counter not in pfam_index_to_neigh_oposite_strand:
                pfam_index_to_neigh_oposite_strand.append(e_coord_counter)
                e_coord_counter += 1
            else:
                e_coord_counter += 1
                continue
        return pfam_index_to_neigh_same_strand, pfam_index_to_neigh_oposite_strand

    def get_list_of_domain_in_neigh(self, pfam_index, file, domains):

        to_counter = []
        party = []
        party.append(file)
        for part in pfam_index:
            party.append(domains[part])
            to_counter.append(domains[part])
        #        neighbourhood_complete.append(party)

        return to_counter

    def counting_dlok_or_dglob(self, some_counter, genome_neigh_size, mighty_domains):
        dlok_glob = []
        for domain_mighty in mighty_domains:
            if domain_mighty in some_counter.keys():
                lok_glob = some_counter.get(domain_mighty) / int(genome_neigh_size)
                dlok_glob.append(lok_glob)
            else:
                dlok_glob.append(0)
        return dlok_glob

    def count_for_dlok(self, somecounter, neigh_size, tax_name, dlok_dataframe):
        for k, v in somecounter.items():
            dlok_dataframe.at[tax_name, k] = v / neigh_size

    def multiple_test_correction(self, some_dataframe, correction_met):

        if correction_met == 'none':
            return some_dataframe
        else:
            value_to_correct = [float(x) for x in some_dataframe.PVALUE.tolist()]
            reject, pvals_corrected, alphaSidak, alphaBonf = correction(pvals=value_to_correct,
                                                                        method=correction_met,
                                                                        is_sorted=False, returnsorted=False)
            pvals = pvals_corrected.tolist()
            some_dataframe.PVALUE = pvals
            return some_dataframe

    def multiple_test_correction_list(self, some_tuple_list, correction_met):
        pvalues = [float(x) for x in some_tuple_list]

        if correction_met == 'none':
            output = pvalues

        else:
            try:

                reject, pvals_corrected, alphaSidak, alphaBonf = correction(pvals=pvalues,
                                                                            method=correction_met,
                                                                            is_sorted=False, returnsorted=False)
                pvals_after_corr = pvals_corrected.tolist()
                output = pvals_after_corr
            except ZeroDivisionError:
                output = pvalues
        return output

    def cutoff_value(self, some_dataframe, cutoff):

        domain_list = list(some_dataframe.index)
        if cutoff == 'none':
            return some_dataframe
        elif cutoff == '0':
            for i in domain_list:

                diff = float(some_dataframe.loc[i, 'PVALUE'])
                if diff < 0 and diff > 0:
                    some_dataframe = some_dataframe.drop([i])
            return some_dataframe
        else:
            cutoff = float(cutoff)
            for i in domain_list:
                diff = float(some_dataframe.loc[i, 'PVALUE'])
                if diff > cutoff:
                    some_dataframe = some_dataframe.drop([i])
            return some_dataframe

    def sort_table(self, some_dataframe):

        sorted_data = some_dataframe.sort_values('Poisson', ascending=True)
        return sorted_data

    def add_information(self, some_dataframe, dictionary):

        indeksy = list(some_dataframe.index)
        for i in indeksy:
            try:

                pfam = i.lower()
                family = dictionary[pfam][0]
                summary = dictionary[pfam][1]

                some_dataframe.at[i, 'Family'] = family
                some_dataframe.at[i, 'Summary'] = summary
            except KeyError:
                continue
        return some_dataframe

    def poisson_distribution_mine(self, k, p, n):
        """
        Return value (float) of poisson distribution mine concept.
            Parameters:
                k (int): number of domain in neighborhoods
                p (float): (sum (average) of neighborhoods lenghts * sum of number of domain) / (average) sum of genome sizes
                n (int): number of neighborhoods
            How to agreggate data:
                k - Counter of all domains should be ok and easy to do
                p - sum of neighborhoods (easy += ) , sum of number of domain (z Countera), sum of genome sizes ( easy +=)
                n - number of neighborhoods (easy +=)
        """
        first_value = (n * p) ** k
        second_value = math.e ** (-1 * n * p)
        third_value = math.factorial(k)
        return (first_value * second_value) / third_value

    def poisson_distribution_KP(self, k, p, n):
        """
        Return value of poisson distribution mine concept
            Parameters:
                k (int): sum of domain occurrence in all neighborhoods
                n (int): number of neighborhoods
                p (float): ( average lengh of neighborhood * number of domain in all database) / sum of lengh of all database
            How to aggregate data:
                k Counter of all domains
                n - number of neigghborhoods (easy +=)
                p - average lengh of neighborhood ( easy +=) / number of genomes, number of domain in all ( from file), sum of lenght of all database ( const value)
        """
        first_value = (n * p) ** k
        second_value = math.e ** (-1 * n * p)
        third_value = math.factorial(k)
        return (first_value * second_value) / third_value

    def pfam_for_pf(self, dataframe):
        """Final customization of dataframe
            --> PF02696 instead of pfam02696
            --> PVALUE in 10e-3 format
            --> In what percentage % format
            --> avg occurences and density  have now 3 digits after coma
            --> drop NOam
        """

        indexy = dataframe.index
        for i in indexy:
            dataframe = dataframe.rename(index={i: i[:2].upper() + i[4:]})

        if "NOam" in indexy:
            dataframe = dataframe.drop(index="NOam")
        # dataframe = dataframe.style.format({'In what percentage': '{:,.3%}'.format})

        return dataframe

    def poisson_distribution(self, p, n, K):

        if int(K) < p * n:
            return 1
        else:
            sums = 0
            k = 0
            number_of_samples = int(K)
            while k < number_of_samples:
                sums += poisson.pmf(k=k, mu=n * p)
                k += 1
            if sums >= 1:
                return float(0.0)
            else:
                return 1 - sums

    def magic_number_creator(self, p0, p1, wie):

        import numpy as np

        np.random.choice([0, 1], p=[p0, p1])

        result = []
        for i in range(wie):
            result.append(np.random.choice([0, 1], p=[p0, p1]))
        return result

    def analysis(self, groups):

        max_neigh = 0

        for org in groups:
            data = pd.read_csv('/home/djangoadmin/final_site-project/media/results/temp/' + org, sep='\t', index_col=0)

            domain_query = 'PF' + self.user_pfam[4:]
            max_neigh = max(max_neigh, data.at[domain_query, 'N'])

        all_dictionary = {}
        all_domains = []
        group = {}
        neigh_c = {}
        genome_c = {}

        for org in groups:

            data = pd.read_csv('/home/djangoadmin/final_site-project/media/results/temp/' + org, sep='\t', index_col=0)

            domains = data.index.tolist()

            neigh_counter = list(data['occurence in neighborhoods'])
            genome_counter = list(data['occurence genomes'])
            k = list(data['K'])
            n = list(data['N'])
            p = list(data['P'])

            # domains
            for domain in domains:
                if domain not in group.keys():
                    group[domain] = 0
                group[domain] += 1
            # K N P
            for i in domains:
                if i not in all_domains:
                    all_domains.append(i)
                    all_dictionary[i] = {'k': [], 'n': [], 'p': [], 'neigh': [], 'genome': []}

            for domain, _k, _n, _p, neigh, genome in zip(domains, k, n, p, neigh_counter, genome_counter):
                all_dictionary[domain]['k'].append(_k)
                all_dictionary[domain]['p'].append(_p)
                all_dictionary[domain]['n'].append(_n)
                all_dictionary[domain]['neigh'].append(neigh)
                all_dictionary[domain]['genome'].append(genome)

        new_dictionary = {}
        for domain in list(all_dictionary.keys()):
            print(domain)
            small_data = all_dictionary[domain]
            new_dictionary[domain] = {'k': [], 'n': [], 'p': [], 'neigh': [], 'genome': []}

            for _k, _n, _p in zip(small_data['k'], small_data['n'], small_data['p']):
                new_dictionary[domain]['p'].append(float(_p))
                rate_of_one = int(_k) / int(_n)
                missing_values = max_neigh - int(_n)
                magic_result = self.magic_number_creator(1 - rate_of_one, rate_of_one, missing_values)
                new_dictionary[domain]['k'].append(int(_k) + sum(magic_result))
                new_dictionary[domain]['n'].append(int(_n) + len(magic_result))
        results = {}
        for domain in list(new_dictionary.keys()):
            k = sum(new_dictionary[domain]['k'])
            n = sum(new_dictionary[domain]['n'])
            p = sum(new_dictionary[domain]['p']) / len(new_dictionary[domain]['p'])
            print(domain)
            results[domain] = self.poisson_distribution(p, n, k)

        return results, all_dictionary

    def merge_all_data_and_save(self, results, all_dict, size_of_group, domain_information, save,
                                overall_num_of_genomes, overall_num_of_neigh):

        complete = pd.DataFrame(columns=['Occurence in genome', 'Avg occurence in genome', 'Occurence in neighborhood',
                                         'Avg occurence in neigh', 'Family percentage',
                                         'Poisson', 'Summary', 'Family'])

        domains = list(results.keys())
        values_to_correct = list(results.values())
        if self.user_correction != 'none':
            reject, pvals_corrected, alphaSidak, alphaBonf = correction(pvals=values_to_correct,
                                                                        method=self.user_correction,
                                                                        is_sorted=False, returnsorted=False)
            pvals = pvals_corrected.tolist()
            for domain, res in zip(domains, pvals):
                complete.at[domain, 'Poisson'] = res
        else:
            for domain, res in zip(domains, values_to_correct):
                complete.at[domain, 'Poisson'] = res

        for domain in domains:
            complete.at[domain, 'Occurence in genome'] = sum(all_dict[domain]['genome'])
            complete.at[domain, 'Avg occurence in genome'] = sum(all_dict[domain]['genome'])

            # complete.at[domain,'Avg occurence in genome'] = sum(all_dict[domain]['genome'])/len(all_dict[domain]['genome'])
            complete.at[domain, 'Occurence in neighborhood'] = sum(all_dict[domain]['neigh']) / overall_num_of_neigh
            # complete.at[domain,'Avg occurence in neigh'] = sum(all_dict[domain]['neigh'])/len(all_dict[domain]['neigh'])
            complete.at[domain, 'Avg occurence in neigh'] = sum(all_dict[domain]['neigh']) / overall_num_of_genomes
            complete.at[domain, 'Family percentage'] = len(all_dict[domain]['neigh']) / size_of_group
        dataframe = self.add_information(complete, domain_information)
        dataframe.index.name = 'Domain'
        dataframe.to_csv('/home/djangoadmin/final_site-project/media/results/' + save, sep='\t', mode='w')

    def go(self):
        groups_calculated = []
        if self.user_level == "genus":

            ####################################################################################################
            # VULTR
            with open('/home/djangoadmin/final_site-project/important_files/genus_level.json', "r") as handler:
                ####################################################################################################
                # LOCAL
                # with open('/mnt/d/45.76.38.24/final_site-project/important_files/genus_level.json', "r") as handler:
                taxonomy = json.load(handler)
            ####################################################################################################

        elif self.user_level == "family":
            ####################################################################################################
            # VULTR
            with open('/home/djangoadmin/final_site-project/important_files/family_level.json', "r") as handler:
                ####################################################################################################
                # LOCAL
                # with open('/mnt/d/45.76.38.24/final_site-project/important_files/family_level.json', "r") as handler:
                taxonomy = json.load(handler)
            ####################################################################################################

        elif self.user_level == "order":
            ####################################################################################################
            # VULTR
            with open('/home/djangoadmin/final_site-project/important_files/order_level.json', "r") as handler:
                ####################################################################################################
                # LOCAL
                # with open('/mnt/d/45.76.38.24/final_site-project/important_files/order_level.json', "r") as handler:
                taxonomy = json.load(handler)
            ####################################################################################################

        elif self.user_level == "class":
            ####################################################################################################
            # VULTR
            with open('/home/djangoadmin/final_site-project/important_files/class_level.json', "r") as handler:
                ####################################################################################################
                # LOCAL
                # with open('/mnt/d/45.76.38.24/final_site-project/important_files/class_level.json', "r") as handler:
                taxonomy = json.load(handler)
            ####################################################################################################

        elif self.user_level == "phylum":
            ####################################################################################################
            # VULTR
            with open('/home/djangoadmin/final_site-project/important_files/phylum_level.json', "r") as handler:
                ####################################################################################################
                # LOCAL
                # with open('/mnt/d/45.76.38.24/final_site-project/important_files/phylum_level.json', "r") as handler:
                taxonomy = json.load(handler)
            ####################################################################################################
        genus_database = list(taxonomy.keys())

        # for self.user_organism in genus_database:

        print(self.user_organism)
        # pdb.set_trace()
        """Zipping all functions and execute them"""
        ################################################################################################################
        # JUST LOCALLY
        # print("checkpoint #1")
        # genome_id_size_in_gene = self.open_multiple_line(
        #     '/mnt/d/45.76.38.24/final_site-project/important_files/GENOME_ID_SIZE_IN_GENE.txt')
        # domain_information = self.collect_pfam_domain_description(
        #     '/mnt/d/45.76.38.24/final_site-project/important_files/domain_information')
        ################################################################################################################

        ################################################################################################################
        # # KLASTER
        #
        # print("checkpoint #1")
        # # KLASTER
        # domain_information = self.collect_pfam_domain_description(
        #     '/home/klaster/ProFaNA v1.0/important_files/domain_information')
        # genome_id_size_in_gene = self.open_multiple_line(
        #     '/home/klaster/ProFaNA v1.0/important_files/GENOME_ID_SIZE_IN_GENE.txt')
        ################################################################################################################
        ################################################################################################################
        # # VULTR
        domain_information = self.collect_pfam_domain_description(
            "/home/djangoadmin/final_site-project/important_files/domain_information")
        genome_id_size_in_gene = self.open_multiple_line(
            "/home/djangoadmin/final_site-project/important_files/GENOME_ID_SIZE_IN_GENE.txt")
        with open('/home/djangoadmin/final_site-project/important_files/counter_domains', 'r') as handler:
            pfam_counter_domains = json.load(handler)
        ################################################################################################################
        print("checkpoint 1")

        genome_id = [x[0] for x in genome_id_size_in_gene]
        size_in_gene = [x[1] for x in genome_id_size_in_gene]

        genome_number_overall = 0
        genome_number_to_stat = 0

        message_for_output = "You have looked for conserved neighbourhood for " + self.user_pfam + \
                             " domain, in range " + str(self.user_distance) + \
                             " bp, in " + self.user_organism + " database."
        message_for_additional_data = "Pfam domain , PVALUE,  occurence in neighbourhoods, " \
                                      "average occurence in neighbourhood ,occurence genomes, " \
                                      "average occurence in genome, avg DLOK-DGLOB, In what percentage?, " \
                                      "Family, Summary"
        print("checkpoint 2")

        ################################################################################################################
        # # LOCAL
        # file_names = self.open_database(self.user_level,self.user_organism)
        ################################################################################################################
        ################################################################################################################
        # KLASTER
        # file_names = self.open_database('/home/klaster/ProFaNA v1.0/important_files/database', self.user_organism)
        ################################################################################################################
        ################################################################################################################
        # VULTR
        level = self.open_database(self.user_level, self.user_organism)
        overall_num_of_genomes = 0
        overall_num_of_neigh = 0

        ################################################################################################################
        for key in list(level.keys()):
            print(key)
            # pdb.set_trace()
            # Create important variables
            just_neigh_data = []
            genomes_with_domain = []
            neigh_genome_size = []
            counter = 0
            counter_db = 0
            sum_of_genes_in_neighborhoods = 0
            number_of_neighborhoods = 0
            percentage_counter = Counter()

            # Dlok operations
            print("DLOK TIME START", datetime.datetime.now())
            for file in level[key]:
                counter_db += 1
                counter += 1
                try:
                    genes, start_coord, end_coord, orientation, domains, contig = self.create_six_list(file)
                except FileNotFoundError:
                    continue
                except ValueError:
                    continue
                file_name_raw = file.split('/')
                tax_name = file_name_raw[-1]
                genome_number_overall += 1
                number_of_genes_in_genome = self.genome_size_in_gene(tax_name, genome_id, size_in_gene)
                coords = self.searching_for_domain_in_genome(self.user_pfam, start_coord, end_coord, orientation,
                                                             domains, contig, genes)
                if len(coords) > 0:
                    genome_number_to_stat += 1
                    overall_num_of_genomes += 1
                for point in coords:
                    number_of_neighborhoods += 1
                    overall_num_of_neigh += 1

                    if self.user_distance != 0:
                        last_coordinate, first_coordinate, pfam_beg, pfam_end, searched_pfam_orientation = \
                            self.get_range_coordinates(point, self.user_distance)
                        pfam_index_to_neigh = self.get_domains_both_strands(start_coord, end_coord, last_coordinate,
                                                                            first_coordinate, pfam_beg, pfam_end,
                                                                            contig,
                                                                            point)
                    else:
                        pfam_index_to_neigh = self.get_domain_both_zeros(genes, point[4])

                    genes_in_neigh = self.list_of_genes_in_neigh(genes, pfam_index_to_neigh)
                    number_of_genes_in_neigh = self.size_in_genes(genes_in_neigh)
                    whole = self.get_list_of_domain_in_neigh(pfam_index_to_neigh, file, domains)
                    percentage_counter += Counter(set(whole))
                    just_neigh_data.append(whole)
                    genomes_with_domain.append(file)
                    neigh_genome_size.append((number_of_genes_in_neigh, number_of_genes_in_genome))
                    sum_of_genes_in_neighborhoods += int(number_of_genes_in_neigh)

            print("DLOK TIME OVER", datetime.datetime.now())
            print("concatenate_domain")

            alls = []
            neigh_counter = Counter()

            for i in just_neigh_data:
                for j in i:
                    if j not in alls:
                        alls.append(j)
            genome_counter = Counter()
            print("MATRIX CREATING", datetime.datetime.now())

            counter = 0
            for domain_list, genome_in_domain, ng_size in zip(just_neigh_data, genomes_with_domain, neigh_genome_size):
                counter += 1
                genes, start_coord, end_coord, orientation, domains, contig = self.create_six_list(genome_in_domain)
                genome_domains_counter = Counter(domains)
                genome_counter += genome_domains_counter
                neigh_domains_counter = Counter(domain_list)
                neigh_counter += neigh_domains_counter

            print("MATRIX OVER", datetime.datetime.now())
            print("Wilcoxon start", datetime.datetime.now())

            najlepsze_dane = pd.DataFrame(columns=['occurence in neighborhoods',
                                                   'occurence genomes',
                                                   'P', 'N', 'K',
                                                   'domain in whole database'])

            counter = 0
            for i in alls:
                counter += 1

                p_KP = ((sum_of_genes_in_neighborhoods / genome_number_to_stat) * pfam_counter_domains[i]) / 326628067

                najlepsze_dane.at[i, 'P'] = p_KP
                najlepsze_dane.at[i, 'domain in whole database'] = pfam_counter_domains[i]
                najlepsze_dane.at[i, 'K'] = int(percentage_counter[i])
                najlepsze_dane.at[i, 'N'] = number_of_neighborhoods
                najlepsze_dane.at[i, 'occurence in neighborhoods'] = neigh_counter[i]
                najlepsze_dane.at[i, 'occurence genomes'] = genome_counter[i]
            print("Wilcoxon OVER", datetime.datetime.now())
            print("Collecting data START", datetime.datetime.now())

            added_information = self.add_information(najlepsze_dane, domain_information)

            print("Collecting data OVER", datetime.datetime.now())
            print("Customizing_data START", datetime.datetime.now())
            final_data = self.pfam_for_pf(added_information)

            print("Customizing_data OVER", datetime.datetime.now())
            message_down = " In selected database there is {} genomes and in {} searched domain " \
                           "was found".format(str(counter_db), str(genome_number_to_stat))
            if len(najlepsze_dane) > 0:
                groups_calculated.append(key)
                self.save_data(key, message_for_output, message_for_additional_data,
                               message_down, final_data)

        print('ANALYSIS NAME')

        complete_results, additional = self.analysis(groups_calculated)

        self.merge_all_data_and_save(complete_results, additional, len(list(level.keys())), domain_information,
                                     self.user_output, overall_num_of_genomes, overall_num_of_neigh)


class SuperSpeedAnalysisFamily:

    def __init__(self, user_pfam, user_distance, user_family, user_cutoff, user_correction,
                 user_strand):
        """ Initializing input data
            Inputs:
                    user_pfam:  str (from pfam00000 to pfam99999)
                    user_distance: str non negative integer 1-20000
                    user_organisms: str
                    user_cutoff: str (none, 0, 0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.005,
                                0.01, 0.05, 0.1, 0.2, 0.5)
                    user_correction: str (none, bonferroni, fdr_bh)
                    user_strand: str (both)
                    user_output: str
        """

        self.user_pfam = user_pfam
        self.user_distance = user_distance
        self.user_family = user_family
        self.user_cutoff = user_cutoff
        self.user_correction = user_correction
        self.user_strand = user_strand
        # self.user_output = user_output

    def open_singleline(self, path_to_file):
        """ Opens file that contain 1 column and strip it by space.
            Inputs:
                    path_to_file: str
        """

        with open(path_to_file, 'r') as f:
            file_names = [line.strip() for line in f]
        return file_names

    def open_singleline_num(self, path_to_file):
        """ Opens file that contain 1 column and strip it by space. """

        with open(path_to_file, 'r') as f:
            file_names = [float(line.strip()) for line in f]
        return file_names

    def open_multiple_line(self, path_to_file):
        """ Opens file that contains data about single genome """

        plik = []
        with open(path_to_file) as inputfile:
            for line in inputfile:
                plik.append(line.strip().split())
        return plik

    def save_data(self, file_name, complete_data):
        """ Saves all stuff together as one file """

        # # LOCAL
        # with open('/mnt/d/45.76.38.24/final_site-project' + file_name + '.txt', 'w') as output_file:
        # VULTR
        # with open('/home/djangoadmin/final_site-project' + file_name, 'w') as output_file:
        #     output_file.write(message_for_output)
        #     output_file.write("\n")
        #     output_file.write(message_down)
        #     output_file.write("\n")
        #     output_file.write(message_for_additional_data)
        #     output_file.write("\n")
        #     # output_file.write(complete_data)
        #     for row in complete_data.iterrows():
        #         index, data = row
        #         pre = data.tolist()
        #         output_file.write(index)
        #         output_file.write("\t")
        #         output_file.write("\t".join([str(i) for i in pre]))
        #         output_file.write("\n")
        complete_data.index.name = 'Domain'
        complete_data.to_csv('/home/djangoadmin/final_site-project' + file_name, sep='\t', mode='w')

    def create_six_list(self, file):
        """ Open single genome file chew it and return 6 lists -> GENE, START_COORD,
        END_COORDS ,ORIENTATION, DOMAINS """
        # KLASTER
        # data = self.open_multiple_line("/home/klaster/ProFaNA v1.0/nowa_baza_danych/"+file)

        # LOCAL
        # data = self.open_multiple_line("/mnt/d/45.76.38.24/final_site-project/important_files/database_file/"+file)

        # VULTR
        data = self.open_multiple_line("/home/djangoadmin/final_site-project/important_files/nowa_baza_danych/" + file)
        start_coord = []
        end_coord = []
        orientation = []
        domains = []
        genes = []
        contig = []
        for bit in data:
            genes.append(int(bit[0]))
            start_coord.append(int(bit[1]))
            end_coord.append(int(bit[2]))
            orientation.append(bit[3])
            domains.append(bit[4])
            contig.append(bit[5])
        return genes, start_coord, end_coord, orientation, domains, contig

    def open_database(self, tax):

        with open('/home/djangoadmin/final_site-project/important_files/database.json', "r") as handler:
            taxonomy = json.load(handler)

        genus = list(taxonomy[tax].keys())
        return genus, taxonomy

    def genome_size_in_gene(self, genome_id, list_of_genome, list_of_size):
        """ Takes genome's id, list of genome and lists with data about how much genes in genome"""

        index_genome = list_of_genome.index(genome_id)
        genome_size_in_gene = list_of_size[index_genome]
        return genome_size_in_gene

    def size_in_genes(self, genome_or_genes):
        """ Takes list of genes and returns size of neighbourhood """

        list_of_genes = []
        for gene in genome_or_genes:
            list_of_genes.append(gene)
        return len(set(list_of_genes))

    def list_of_genes_in_neigh(self, list_of_genes, index_list):
        """ Takes overall gene's list in genome and creates complete list of genes in neighbourhood """

        list_of_genes_in_neigh = []
        for pfam_domain in index_list:
            list_of_genes_in_neigh.append(list_of_genes[pfam_domain])
        return list_of_genes_in_neigh

    def searching_for_domain_in_genome(self, pfam, start_coord, end_coord, orient, domains, contig, genes):
        """
        Takes user's pfam domain searches through lists contains coordinates,
        orientation, domains, contigs and returns complete data about all
        users' domains in genome
        """
        coords = []
        coords_counter = 0
        for domain in domains:
            one_coords_part = []
            if domain == pfam:
                orientation_pfam_domain = orient[coords_counter]
                one_coords_part.append(start_coord[coords_counter])
                one_coords_part.append(end_coord[coords_counter])
                one_coords_part.append(orientation_pfam_domain)
                one_coords_part.append(contig[coords_counter])
                one_coords_part.append(genes[coords_counter])
                coords.append(one_coords_part)
                coords_counter += 1
            else:
                coords_counter += 1
                continue
        return coords

    def presence_confirmation(self, coords, file, pfam):
        """
                Just prints out information how many pfam domains is in each genome
                during analysis propably i will have to remove it before putting it to website ...
                """

        if len(coords) == 0:
            print("In genome  " + file + " pfam domain you have been searching do not exist")
        elif len(coords) == 1:
            print("In genome  " + file + " there is " + str(len(coords)) + " " + pfam + " domain")
        else:
            print("In genome  " + file + " there are " + str(len(coords)) + " " + pfam + " domains")

    def get_range_coordinates(self, target_pfam, distance):
        """
        Based on how large neighbourhood user wants to analyze creates
        points FROM and TO, additionaly shows where main user's pfam begins and ENDS
        with orientation on strands
        """

        last_coordinate = target_pfam[1] + int(distance)
        first_coordinate = target_pfam[0] - int(distance)
        pfam_beg = target_pfam[0]
        pfam_end = target_pfam[1]
        searched_pfam_orientation = target_pfam[2]
        return last_coordinate, first_coordinate, pfam_beg, pfam_end, searched_pfam_orientation

    def get_domain_both_zeros(self, ref_gene, q_gene):

        index_list_gene = []
        counter = 0
        for gene in ref_gene:
            if gene == int(q_gene):
                index_list_gene.append(counter)
            counter += 1
        return index_list_gene

    def get_domains_both_strands(self, start_coord, end_coord, last_coordinate, first_coordinate, pfam_beg,
                                 pfam_end, contig, point):

        pfam_index_to_neigh = []
        s_coord_counter = 0
        for s_coord in start_coord:
            if s_coord <= last_coordinate and s_coord >= pfam_beg and contig[s_coord_counter] == point[
                3] and s_coord_counter not in pfam_index_to_neigh:
                pfam_index_to_neigh.append(s_coord_counter)
                s_coord_counter += 1
            else:
                s_coord_counter += 1
                continue
        e_coord_counter = 0
        for e_coord in end_coord:
            if e_coord >= first_coordinate and e_coord <= pfam_end and contig[e_coord_counter] == point[
                3] and e_coord_counter not in pfam_index_to_neigh:
                pfam_index_to_neigh.append(e_coord_counter)
                e_coord_counter += 1
            else:
                e_coord_counter += 1
                continue
        return pfam_index_to_neigh

    def collect_pfam_domain_description(self, file):
        """ Opens file only to gather information"""

        data = {}
        with open(file) as inputfile:
            for line in inputfile:
                try:

                    one_line = line.strip().split("@")
                    domena = one_line[0]
                    pre_family = one_line[1][8:]
                    delete_this = "(" + domena.upper() + ")"
                    family = pre_family.replace(delete_this, "")
                    summary = one_line[2][9:]
                    data[domena] = (family, summary)

                except IndexError:
                    continue
        return data

    def get_domains_plus_minus_strand(self, start_coord, end_coord, last_coordinate, first_coordinate,
                                      pfam_beg, pfam_end, contig, point, orientation):

        pfam_index_to_neigh_same_strand = []
        pfam_index_to_neigh_oposite_strand = []
        s_coord_counter = 0
        for s_coord in start_coord:
            if s_coord <= last_coordinate and s_coord >= pfam_beg and orientation[s_coord_counter] == point[2] and \
                    contig[s_coord_counter] == point[3] and s_coord_counter not in pfam_index_to_neigh_same_strand:
                pfam_index_to_neigh_same_strand.append(s_coord_counter)
                s_coord_counter += 1
            elif s_coord <= last_coordinate and s_coord >= pfam_beg and orientation[s_coord_counter] != point[2] and \
                    contig[s_coord_counter] == point[
                3] and s_coord_counter not in pfam_index_to_neigh_oposite_strand:
                pfam_index_to_neigh_oposite_strand.append(s_coord_counter)
                s_coord_counter += 1
            else:
                s_coord_counter += 1
                continue
        e_coord_counter = 0
        for e_coord in end_coord:
            if e_coord >= first_coordinate and e_coord <= pfam_end and orientation[e_coord_counter] == point[2] and \
                    contig[e_coord_counter] == point[3] and e_coord_counter not in pfam_index_to_neigh_same_strand:
                pfam_index_to_neigh_same_strand.append(e_coord_counter)
                e_coord_counter += 1
            elif e_coord >= first_coordinate and e_coord <= pfam_end and orientation[e_coord_counter] != point[
                2] and contig[e_coord_counter] == point[
                3] and e_coord_counter not in pfam_index_to_neigh_oposite_strand:
                pfam_index_to_neigh_oposite_strand.append(e_coord_counter)
                e_coord_counter += 1
            else:
                e_coord_counter += 1
                continue
        return pfam_index_to_neigh_same_strand, pfam_index_to_neigh_oposite_strand

    def get_list_of_domain_in_neigh(self, pfam_index, file, domains):

        to_counter = []
        party = []
        party.append(file)
        for part in pfam_index:
            party.append(domains[part])
            to_counter.append(domains[part])
        #        neighbourhood_complete.append(party)

        return to_counter

    def counting_dlok_or_dglob(self, some_counter, genome_neigh_size, mighty_domains):
        dlok_glob = []
        for domain_mighty in mighty_domains:
            if domain_mighty in some_counter.keys():
                lok_glob = some_counter.get(domain_mighty) / int(genome_neigh_size)
                dlok_glob.append(lok_glob)
            else:
                dlok_glob.append(0)
        return dlok_glob

    def count_for_dlok(self, somecounter, neigh_size, tax_name, dlok_dataframe):
        for k, v in somecounter.items():
            dlok_dataframe.at[tax_name, k] = v / neigh_size

    def multiple_test_correction(self, some_dataframe, correction_met):

        if correction_met == 'none':
            return some_dataframe
        else:
            value_to_correct = [float(x) for x in some_dataframe.PVALUE.tolist()]
            reject, pvals_corrected, alphaSidak, alphaBonf = correction(pvals=value_to_correct,
                                                                        method=correction_met,
                                                                        is_sorted=False, returnsorted=False)
            pvals = pvals_corrected.tolist()
            some_dataframe.PVALUE = pvals
            return some_dataframe

    def multiple_test_correction_list(self, some_tuple_list, correction_met):
        pvalues = [float(x) for x in some_tuple_list]

        if correction_met == 'none':
            output = pvalues

        else:
            try:

                reject, pvals_corrected, alphaSidak, alphaBonf = correction(pvals=pvalues,
                                                                            method=correction_met,
                                                                            is_sorted=False, returnsorted=False)
                pvals_after_corr = pvals_corrected.tolist()
                output = pvals_after_corr
            except ZeroDivisionError:
                output = pvalues
        return output

    def cutoff_value(self, some_dataframe, cutoff):

        domain_list = list(some_dataframe.index)
        if cutoff == 'none':
            return some_dataframe
        elif cutoff == '0':
            for i in domain_list:
                diff = float(some_dataframe.loc[i, 'Poisson_correct'])
                if diff < 0 and diff > 0:
                    some_dataframe = some_dataframe.drop([i])
            return some_dataframe
        else:
            cutoff = float(cutoff)
            for i in domain_list:
                diff = float(some_dataframe.loc[i, 'Poisson_correct'])
                if diff > cutoff:
                    some_dataframe = some_dataframe.drop([i])
            return some_dataframe

    def sort_table(self, some_dataframe):

        sorted_data = some_dataframe.sort_values('Poisson', ascending=True)
        return sorted_data

    def add_information(self, some_dataframe, dictionary):

        indeksy = list(some_dataframe.index)
        for i in indeksy:
            try:

                pfam = i[0:2] + i[4:]
                family = dictionary[pfam][0]
                summary = dictionary[pfam][1]

                some_dataframe.at[i, 'Family'] = family
                some_dataframe.at[i, 'Summary'] = summary
            except KeyError:
                continue
        return some_dataframe

    def poisson_distribution_KP(self, k, p, n):
        """
        Return value of poisson distribution mine concept
            Parameters:
                k (int): sum of domain occurrence in all neighborhoods
                n (int): number of neighborhoods
                p (float): ( average lengh of neighborhood * number of domain in all database) / sum of lengh of all database
            How to aggregate data:
                k Counter of all domains
                n - number of neigghborhoods (easy +=)
                p - average lengh of neighborhood ( easy +=) / number of genomes, number of domain in all ( from file), sum of lenght of all database ( const value)
        """
        first_value = (n * p) ** k
        second_value = math.e ** (-1 * n * p)
        third_value = math.factorial(k)
        return (first_value * second_value) / third_value

    def poisson_distribution_KP_mod(self, k, p, n):
        """
        Return value of poisson distribution mine concept
            Parameters:
                k (int): sum of domain occurrence in all neighborhoods
                n (int): number of neighborhoods
                p (float): ( average lengh of neighborhood * number of domain in all database) / sum of lengh of all database
            How to aggregate data:
                k Counter of all domains
                n - number of neigghborhoods (easy +=)
                p - average lengh of neighborhood ( easy +=) / number of genomes, number of domain in all ( from file), sum of lenght of all database ( const value)
        """

        if n * p < k:

            first_value = (n * p) ** k
            second_value = math.e ** (-1 * n * p)
            third_value = math.factorial(k)
            return (first_value * second_value) / third_value
        else:
            return 1.0

    def pfam_for_pf(self, dataframe):
        """Final customization of dataframe
            --> PF02696 instead of pfam02696
            --> PVALUE in 10e-3 format
            --> In what percentage % format
            --> avg occurences and density  have now 3 digits after coma
            --> drop NOam
        """

        indexy = dataframe.index
        for i in indexy:
            dataframe = dataframe.rename(index={i: i[:2].upper() + i[4:]})
        dataframe['Poisson'] = dataframe['Poisson'].map('{:.3e}'.format)
        dataframe['Poisson_correct'] = dataframe['Poisson_correct'].map('{:.3e}'.format)

        dataframe['In what percentage'] = dataframe['In what percentage'].map('{:.3%}'.format)
        dataframe['average occurence in neighborhood'] = dataframe['average occurence in neighborhood'].map(
            '{:.3}'.format)
        dataframe['average occurence in genome'] = dataframe['average occurence in genome'].map('{:.3}'.format)

        if "NOam" in indexy:
            dataframe = dataframe.drop(index="NOam")
        # dataframe = dataframe.style.format({'In what percentage': '{:,.3%}'.format})
        dataframe.drop('average occurence in neighborhood', axis=1, inplace=True)
        dataframe.drop('average occurence in genome', axis=1, inplace=True)
        dataframe.drop('avg neigh size', axis=1, inplace=True)
        dataframe.drop('domain in whole database', axis=1, inplace=True)
        if self.user_correction == 'none':
            dataframe.drop('Poisson_correct', axis=1, inplace=True)
        else:
            dataframe.drop('Poisson', axis=1, inplace=True)

        return dataframe

    def poisson_distribution(self, p, n, domain, percentage_counter):

        if int(percentage_counter[domain]) < p * n:
            return 1
        else:
            sums = 0
            k = 0
            number_of_samples = int(percentage_counter[domain])
            while k < number_of_samples:
                sums += poisson.pmf(k=k, mu=n * p)
                k += 1
            if 1 - sums >= 0.0:
                return 1 - sums
            else:
                return 0.0
            # return 1 - sums

            # return 1 - sums
            # stat = poisson.pmf(k=int(percentage_counter[domain]), mu=n * p)
            # return stat

    def go(self):
        print(self.user_family)
        domain_information = self.collect_pfam_domain_description(
            "/home/djangoadmin/final_site-project/important_files/domain_information")
        genome_id_size_in_gene = self.open_multiple_line(
            "/home/djangoadmin/final_site-project/important_files/GENOME_ID_SIZE_IN_GENE.txt")
        with open('/home/djangoadmin/final_site-project/important_files/counter_domains', 'r') as handler:
            pfam_counter_domains = json.load(handler)

        genome_id = [x[0] for x in genome_id_size_in_gene]
        size_in_gene = [x[1] for x in genome_id_size_in_gene]

        genome_number_overall = 0
        genome_number_to_stat = 0
        print("checkpoint 2")

        genus_database, taxonomy = self.open_database(self.user_family)

        genome_number_overall = 0
        genome_number_to_stat = 0

        for genus in genus_database:
            genomes = taxonomy[self.user_family][genus]
            new_tax = "_".join(genus.split(' '))
            temp_data_link = '/media/results/temp/' + new_tax

            just_neigh_data = []
            genomes_with_domain = []
            neigh_genome_size = []
            counter = 0
            counter_db = 0
            sum_of_genes_in_neighborhoods = 0
            number_of_neighborhoods = 0
            percentage_counter = Counter()
            print("DLOK TIME START", datetime.datetime.now())
            for file in genomes:
                " OPEN FILE "
                # print(file)
                counter_db += 1
                counter += 1
                try:
                    genes, start_coord, end_coord, orientation, domains, contig = self.create_six_list(file)
                    " DZIELIMY PLIK NA KAWALKI"

                except FileNotFoundError:
                    continue
                except ValueError:
                    continue
                file_name_raw = file.split('/')
                tax_name = file_name_raw[-1]
                genome_number_overall += 1
                number_of_genes_in_genome = self.genome_size_in_gene(tax_name, genome_id, size_in_gene)
                coords = self.searching_for_domain_in_genome(self.user_pfam, start_coord, end_coord, orientation,
                                                             domains, contig, genes)
                if len(coords) > 0:
                    genome_number_to_stat += 1

                for point in coords:
                    number_of_neighborhoods += 1
                    if self.user_distance != 0:

                        last_coordinate, first_coordinate, pfam_beg, pfam_end, searched_pfam_orientation = \
                            self.get_range_coordinates(point, self.user_distance)

                        if self.user_strand == 'both':
                            pfam_index_to_neigh = self.get_domains_both_strands(start_coord, end_coord, last_coordinate,
                                                                                first_coordinate, pfam_beg, pfam_end,
                                                                                contig,
                                                                                point)
                        else:
                            pfam_index_to_neigh_same, pfam_index_to_neigh_oposite = self.get_domains_plus_minus_strand(
                                start_coord, end_coord, last_coordinate,
                                first_coordinate, pfam_beg, pfam_end,
                                contig, point, orientation)
                            if self.user_strand == 'same':
                                pfam_index_to_neigh = pfam_index_to_neigh_same
                            elif self.user_strand == 'opposite':
                                pfam_index_to_neigh = pfam_index_to_neigh_oposite
                    else:
                        pfam_index_to_neigh = self.get_domain_both_zeros(genes, point[4])

                    genes_in_neigh = self.list_of_genes_in_neigh(genes, pfam_index_to_neigh)
                    number_of_genes_in_neigh = self.size_in_genes(genes_in_neigh)
                    whole = self.get_list_of_domain_in_neigh(pfam_index_to_neigh, file, domains)
                    percentage_counter += Counter(set(whole))
                    just_neigh_data.append(whole)
                    genomes_with_domain.append(file)
                    neigh_genome_size.append((number_of_genes_in_neigh, number_of_genes_in_genome))
                    sum_of_genes_in_neighborhoods += int(number_of_genes_in_neigh)

            print("DLOK TIME OVER", datetime.datetime.now())
            print("concatenate_domain")

            alls = []
            neigh_counter = Counter()

            for i in just_neigh_data:
                for j in i:
                    if j not in alls:
                        alls.append(j)
            genome_counter = Counter()
            print("MATRIX CREATING", datetime.datetime.now())

            counter = 0
            for domain_list, genome_in_domain, ng_size in zip(just_neigh_data, genomes_with_domain, neigh_genome_size):
                counter += 1
                genes, start_coord, end_coord, orientation, domains, contig = self.create_six_list(genome_in_domain)
                genome_domains_counter = Counter(domains)
                genome_counter += genome_domains_counter
                neigh_domains_counter = Counter(domain_list)
                neigh_counter += neigh_domains_counter

            print("MATRIX OVER", datetime.datetime.now())
            print("Wilcoxon start", datetime.datetime.now())

            najlepsze_dane = pd.DataFrame(columns=['occurence in neighborhoods',
                                                   'average occurence in neighborhood', 'occurence genomes',
                                                   'average occurence in genome',
                                                   'In what percentage', 'Poisson', 'Poisson_correct', 'P', 'N', 'K',
                                                   'avg neigh size', 'domain in whole database',
                                                   'Family', 'Summary', ])

            scores_poisson = []
            counter = 0
            for i in alls:
                counter += 1

                p_KP = ((sum_of_genes_in_neighborhoods / genome_number_to_stat) * pfam_counter_domains[i]) / 326628067

                najlepsze_dane.at[i, 'P'] = p_KP
                najlepsze_dane.at[i, 'avg neigh size'] = sum_of_genes_in_neighborhoods / genome_number_to_stat
                najlepsze_dane.at[i, 'domain in whole database'] = pfam_counter_domains[i]
                # najlepsze_dane.at[i,'sum of database len'] = 326628067
                poisson_test_result = self.poisson_distribution(p_KP, number_of_neighborhoods, i, percentage_counter)
                najlepsze_dane.at[i, 'K'] = int(percentage_counter[i])
                najlepsze_dane.at[i, 'N'] = number_of_neighborhoods
                najlepsze_dane.at[i, 'P'] = p_KP
                najlepsze_dane.at[i, 'Poisson'] = poisson_test_result
                scores_poisson.append(poisson_test_result)

            # scores_after_correction = self.multiple_test_correction_list(scores,self.user_correction)
            poisson_after_correction = self.multiple_test_correction_list(scores_poisson, self.user_correction)
            print("Wilcoxon OVER", datetime.datetime.now())
            print("Collecting data START", datetime.datetime.now())

            for domena, poisson_corr in zip(alls, poisson_after_correction):
                najlepsze_dane.at[domena, 'Poisson_correct'] = poisson_corr
                najlepsze_dane.at[domena, 'In what percentage'] = percentage_counter[domena] / len(just_neigh_data)
                najlepsze_dane.at[domena, 'occurence in neighborhoods'] = neigh_counter[domena]
                najlepsze_dane.at[domena, 'average occurence in neighborhood'] = neigh_counter[domena] / len(
                    just_neigh_data)
                najlepsze_dane.at[domena, 'occurence genomes'] = genome_counter[domena]
                najlepsze_dane.at[domena, 'average occurence in genome'] = genome_counter[domena] / len(
                    set(genomes_with_domain))

            print("Collecting data OVER", datetime.datetime.now())
            print("Customizing_data START", datetime.datetime.now())
            after_cutoff = self.cutoff_value(najlepsze_dane, self.user_cutoff)
            sorted_table = self.sort_table(after_cutoff)
            added_information = self.add_information(sorted_table, domain_information)
            final_data = self.pfam_for_pf(added_information)

            print("Customizing_data OVER", datetime.datetime.now())
            message_down = " In selected database there is {} genomes and in {} searched domain " \
                           "was found".format(str(counter_db), str(genome_number_to_stat))
            self.save_data(temp_data_link, final_data)

