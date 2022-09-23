import streamlit as st
from io import StringIO
from .utils import NeighborhoodAnalyzerFromGene,SuperSpeedAnalysisFromDomain
import pandas as pd
import subprocess
from statsmodels.stats.multitest import multipletests as correction  #
import datetime
def profana_site():
    """
    Function to write what this app is about
    :return:None
    """

    st.header("Profana test site")

    st.subheader("Analysis from gene list")

    with st.form("gene list analysis"):

        domains = st.file_uploader(label="Please load your genes")
        neigh_size = st.number_input("Neighborhood size up and downstream")



        submit = st.form_submit_button("Start analysis")
        if submit:
            st.write(f" your input is  {int(neigh_size)}")
            stringio = StringIO(domains.getvalue().decode("utf-8"))
            string_data = stringio.read().strip().split()
            with open('modules/data/input','w') as handler:
                for i in string_data:
                    handler.write(i)
                    handler.write('\n')
            NeighborhoodAnalyzerFromGene(user_list_of_genes='modules/data/input',user_distance_value=neigh_size,
                                     user_database=None, user_cutoff='none',user_correction='bonferroni',
                                     user_strand_value='none',user_output='modules/data/out').go()
            data = pd.read_csv('modules/data/out',delimiter='\t')
            st.dataframe(data)

    st.subheader("Analysis from domain - all_database")
    with st.form("domain_analysis"):
        domains = st.text_input(label="domain")
        neigh_size = st.number_input("Neighborhood size up and downstream")
        submit = st.form_submit_button("Start analysis")
        if submit:
            st.write(f'START {datetime.datetime.now()}')
            subprocess.run(f"grep -r {domains} modules/data/img_ready/ > modules/data/tmp_list", shell=True  )
            with open('modules/data/tmp_list','r') as handler:
                data = [x.strip().split() for x in handler]
            genes = [x[0].split(':')[1] for x in data]
            with open('modules/data/input','w') as fs:
                for i in genes:
                    fs.write(i)
                    fs.write('\n')

            NeighborhoodAnalyzerFromGene(user_list_of_genes='modules/data/input',user_distance_value=neigh_size,
                                     user_database=None, user_cutoff='none',user_correction='bonferroni',
                                     user_strand_value='none',user_output='modules/data/out').go()
            st.write(f'END {datetime.datetime.now()}')
            # SuperSpeedAnalysisFromDomain(user_pfam='pfam02696', user_distance=neigh_size, user_correction='bonferroni', user_strand='both',
            #                              user_output='modules/data/out', user_cutoff='none', skip_negative='yes',user_organisms='all genomes').go()

# python3 /usr/src/app/scripts/execute_order_66.py pfam13271 0 all_genomes none bonferroni both hoasxbwpidtfdic.txt yes