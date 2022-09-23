# import os
# from collections import Counter
# import json
#
# files = os.listdir('data/img_ready')
# all_count = Counter()
# for file in files:
#     print(file)
#     with open(f'data/img_ready/{file}','r') as handler:
#         data = [x.strip().split() for x in handler]
#     domains = [x[4] for x in data if x != 'noPfam']
#     small_counter = Counter(domains)
#     all_count += small_counter
#
# with open('data/counter_domains_new', 'w') as fp:
#     json.dump(all_count, fp)
# ####
import os
files = os.listdir('data/img_ready')
#
# for file in files:
#     print(file)
#     with open(f'data/img_ready/{file}','r') as handler:
#         data = [x.strip().split() for x in handler]
#     contigs = list(set([x[5] for x in data]))
#     result = []
#     result.append(str(file))
#     for contig in contigs:
#         small_coord = []
#         for line in data:
#             if line[5] == contig:
#                 small_coord.append(int(line[0]))
#         genes_start = str(min(small_coord))
#         genes_end = str(max(small_coord))
#         result.append(genes_start)
#         result.append(genes_end)
#
#     with open('data/genomes_map_new','a+') as fs:
#         for i in result:
#             fs.write(i)
#             fs.write(' ')
#         fs.write('\n')

import pdb
# files = [2905846765,2905846765]

for file in files:
    # print(file)
    with open(f'data/img_ready/{file}','r') as handler:
        data = [x.strip().split() for x in handler]
    just_genes = list(set([int(x[0]) for x in data]))
    just_genes.sort()
    # print(just_genes)
    groups = []
    first_gene = int(just_genes[0])
    gene_temp = int(first_gene)
    for line in just_genes:
        if line == first_gene:
            continue
        elif int(gene_temp)+1 == int(line) and line != just_genes[-1]:
            gene_temp += 1
            continue
        elif int(gene_temp)+1 != int(line):
            groups.append((first_gene, gene_temp))
            first_gene = int(line)
            gene_temp = int(line)
        elif line == just_genes[-1]:
            groups.append((first_gene, line))
    # print(first_gene, line,just_genes[-1])
    with open('data/genomes_map_new','a+') as fs:
        fs.write(file)
        fs.write(' ')
        for i in groups:
            for j in i:
                fs.write(str(j))
                fs.write(' ')
        fs.write('\n')
# print(groups)










###
# import os
# files = os.listdir('data/img_ready')
#
# for file in files:
#     print(file)
#     with open(f'data/img_ready/{file}','r') as handler:
#         data = [x.strip().split() for x in handler]
#     genes = [int(x[0]) for x in data]
#     liczba_genow = len(set(genes))
#     with open('data/GENOME_ID_SIZE_IN_GENE_NEW.txt','a+') as fs:
#         fs.write(f"{file} {liczba_genow}\n")