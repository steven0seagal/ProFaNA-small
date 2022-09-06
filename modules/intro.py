import streamlit as st
from PIL import Image

def intro():
    """
    Function to write what this app is about
    :return:None
    """
    image = Image.open('drawing.png')
    st.image(image, width=256)

    st.header("Machine Learning project")

    st.markdown(
        """
        This is project to collect and use machine learning algorithms.

        Currently available algorithms:


        **👈 Select a section from the dropdown on the left** !
        """
    )