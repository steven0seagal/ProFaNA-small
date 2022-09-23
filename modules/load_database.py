import psycopg2
import streamlit as st
import multiprocessing

def load_database():
    st.header("Load database")

    with st.form("load genome sizes"):
        submit = st.form_submit_button("UPLOAD SIZES")
        if submit:
            load_genome_sizes()
    with st.form("load domain information"):
        submit = st.form_submit_button("UPLOAD INFORMATION")
        if submit:
            load_domain_information()
    with st.form("load map"):
        submit = st.form_submit_button("UPLOAD MAP")
        if submit:
            load_map()

    with st.form("test query"):
        gene  = st.text_input(label="Gene input")
        submit = st.form_submit_button("SEARCH FOR GENOME")
        if submit:
            find_correct_genome(gene)


def find_correct_genome(gene):
    conn = psycopg2.connect(
        host="172.17.0.1",
        database="mydatabase",
        user="postgres",
        password="postgres")
    cursor = conn.cursor()
    query = f"SELECT genome, begins,ends FROM genomes_map WHERE {int(gene)} BETWEEN begins AND ends"

    # postgreSQL_insert_Query = "select * from test"
    cursor.execute(query)
    records = cursor.fetchall()
    print(records)


def load_genome_sizes():

    with open("modules/data/GENOME_ID_SIZE_IN_GENE.txt",'r') as handler:
        data = [x.strip().split() for x in handler]

    conn = psycopg2.connect(
        host="172.17.0.1",
        database="mydatabase",
        user="postgres",
        password="postgres")
    cursor = conn.cursor()

    for line in data:
        print(line[0])

        query = f"INSERT INTO genome_sizes (img_id, size_in_genes) VALUES({line[0]}, {line[1]})"

        # postgreSQL_insert_Query = "select * from test"
        cursor.execute(query)
        conn.commit()

def load_domain_information():

    with open("modules/data/domain_information", 'r') as handler:
        data = [x.strip().split('\t') for x in handler]

    conn = psycopg2.connect(
        host="172.17.0.1",
        database="mydatabase",
        user="postgres",
        password="postgres")
    cursor = conn.cursor()

    for line in data:
        print(line)
        query = f"INSERT INTO domain_information (domains, family, summary) VALUES({line[0]}, {str(line[1])}, {line[2]})"

        # postgreSQL_insert_Query = "select * from test"
        cursor.execute(query)
        conn.commit()

def load_map():
    with open("modules/data/genomes_map_new", 'r') as handler:
        data = [x.strip().split() for x in handler]
    p = multiprocessing.Pool(1)
    p.map(mp_worker, data)


    # conn = psycopg2.connect(
    #     host="172.17.0.1",
    #     database="mydatabase",
    #     user="postgres",
    #     password="postgres")
    # cursor = conn.cursor()
    #
    #
    #
    #
    #
    #
    # for line in data:
    #     genome = line[0]
    #     begins = line[1::2]
    #     ends = line[2::2]
    #     for begin, end in zip(begins,ends):
    #         print(genome, begin, end)
    #         query = f"INSERT INTO genomes_map (genome, begins, stop) VALUES({genome}, {begin}, {end})"
    #
    #         cursor.execute(query)
    #         conn.commit()
    #

def mp_worker(line):
    conn = psycopg2.connect(
        host="172.17.0.1",
        database="mydatabase",
        user="postgres",
        password="postgres")
    cursor = conn.cursor()

    genome = line[0]
    begins = line[1::2]
    ends = line[2::2]
    for begin, end in zip(begins,ends):
        print(genome, begin, end)
        # querySelect = f"SELECT genome, begins, ends FROM public.genomes_map WHERE genome = {genome} AND begins ={begin}  AND ends = {end} "
        # cursor.execute(querySelect)
        # select_records = cursor.fetchall()
        # print(len(select_records))
        query = f"INSERT INTO genomes_map (genome, begins, ends) VALUES({genome}, {begin}, {end})"

        cursor.execute(query)
        conn.commit()



