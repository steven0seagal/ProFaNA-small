import streamlit as st
from PIL import Image
import psycopg2

import streamlit as st
from modules.intro import intro



pages = {
    "Intro":intro,
}

def test():
    conn = psycopg2.connect(
        host="172.17.0.1",
        database="mydatabase",
        user="postgres",
        password="postgres")
    cursor = conn.cursor()
    postgreSQL_select_Query = "select * from test"
    cursor.execute(postgreSQL_select_Query)
    mobile_records = cursor.fetchall()
    for row in mobile_records:
        print("Id = ", row[0], )
test()



demo_name = st.sidebar.selectbox("Choose an algorithm", pages.keys())
pages[demo_name]()


