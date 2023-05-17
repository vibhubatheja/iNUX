import numpy as np
import matplotlib.pyplot as plt
import scipy.special
import os
import csv
import streamlit as st 
from PIL import Image
from io import BytesIO
import io
import base64
import pandas as pd
#import plotly.graph_objects as go

st.set_page_config(
    layout="centered", page_title="Theis  Type Curve", page_icon="📈"
)
st.markdown(
    """
    <style>
    body {
        background-color: ##00c69b;
    }
    </style>
    """,
    unsafe_allow_html=True
)




# DATA FROM FETTER (TABLE 5.1), TRANSFERED FROM FEET TO METER
# Q = 42400 ft3/d, r = 824 ft
m_time = [3,5,8,12,20,24,30,38,47,50,60,70,80,90,100,130,160,200,260,320,380,500] # time in minutes
m_ddown = [0.093,0.216,0.401,0.648,0.987,1.110,1.264,1.449,1.573,1.635,1.758,1.881,1.943,2.066,2.159,2.313,2.560,2.621,2.837,2.991,3.146,3.362]   # drawdown in meters

# DATA FROM SCHWARTZ & ZHANG
# Q_SCHandZH = 500 m3/d, r = 300 m
#m_time = [1,1.27,1.61,2.04,2.59,3.29,4.18,5.3,6.72,8.53,10.83,13.74,17.43,22.12,28.07,35.62,45.2,57.36,72.79,92.37,117.21,148.74,188.74,239.5,303.92,385.66,489.39,621.02,788.05,1000] # time in minutes
#m_ddown = [0.03,0.05,0.09,0.15,0.22,0.31,0.41,0.53,0.66,0.8,0.95,1.11,1.27,1.44,1.61,1.79,1.97,2.15,2.33,2.52,2.7,2.89,3.07,3.26,3.45,3.64,3.83,4.02,4.21,4.39]   # drawdown in meters



#Define a function, class, and object for Theis Well analysis

def well_function(u):
    return scipy.special.exp1(u)

def theis_u(T,S,r,t):
    u = r ** 2 * S / 4. / T / t
    return u

def theis_s(Q, T, u):
    s = Q / 4. / np.pi / T * well_function(u)
    return s

def theis_wu(Q, T, s):
    wu = s * 4. * np.pi * T / Q
    return wu




def intro() :
    
    st.title(":red[Theis Type Curve & Data sheets]")
    """

This notebook is intended to draw the necessary plots for graphical pumping test evaluation methods. Plots will be saved in Cache memory.


## Introduction

### General situation
We consider a confined aquifer with constant transmissivity. If a well is pumping water out of the aquifer, radial flow towards the well is induced. The hydraulics are described by the following simplified flow equation. This equation accounts for 1D radial transient flow towards a fully penetrating well within a confined aquifer without further sinks and sources:
"""
    
    latex=r'''$$ \frac{\partial^2 h}{\partial r^2}+\frac{1}{r}\frac{\partial h}{\partial r}=\frac{S}{T}\frac{\partial h}{\partial t} $$'''
    st.write(latex)
    """
### Solution by Theis
Charles V. Theis presented a solution for this by deriving
"""
    
    latex=(r'''$$ s(r,t)=\frac{Q}{4\pi T}W(u) $$''')
    st.write(latex)
    """
with the well function
"""
    latex=(r'''$$ W(u) = \int_{u }^{+\infty} \frac{e^{-\tilde u}}{\tilde u}d\tilde u $$''')
    st.write(latex)
    """
and the dimensionless variable 
"""
    latex=(r'''$$ u = \frac{Sr^2}{4Tt} $$''')
    st.write(latex)
    """
This equations are not easy to solve. Historically, values for the well function were provided by tables or as so called type-curve. The type-curve matching with experimental data for pumping test analysis can be considered as one of the basic hydrogeological methods.

However, modern computer provide an easier and more convinient way to solve the 1D radial flow equation based on the Theis approach.


"""

def datainput() :
    
    dc=st.selectbox("Select the Data set for computation",("Data from FETTER (TABLE 5.1)"," Data from SCHWARTZ & ZHANG","Data from CSV File"))
    
    if dc=="Data from FETTER (TABLE 5.1)" :
        m_time = [3,5,8,12,20,24,30,38,47,50,60,70,80,90,100,130,160,200,260,320,380,500] # time in minutes
        m_ddown = [0.093,0.216,0.401,0.648,0.987,1.110,1.264,1.449,1.573,1.635,1.758,1.881,1.943,2.066,2.159,2.313,2.560,2.621,2.837,2.991,3.146,3.362]   # drawdown in meters

        
        
        st.session_state.m_time = m_time
        st.session_state.m_ddown =  m_ddown
    
    if dc==" Data from SCHWARTZ & ZHANG" :
        st.write("Question is it working inside 2 ?")
        m_time = [1,1.27,1.61,2.04,2.59,3.29,4.18,5.3,6.72,8.53,10.83,13.74,17.43,22.12,28.07,35.62,45.2,57.36,72.79,92.37,117.21,148.74,188.74,239.5,303.92,385.66,489.39,621.02,788.05,1000] # time in minutes
        m_ddown = [0.03,0.05,0.09,0.15,0.22,0.31,0.41,0.53,0.66,0.8,0.95,1.11,1.27,1.44,1.61,1.79,1.97,2.15,2.33,2.52,2.7,2.89,3.07,3.26,3.45,3.64,3.83,4.02,4.21,4.39]   # drawdown in meters

        
        
        st.session_state.m_time = m_time
        st.session_state.m_ddown =  m_ddown


    
    if dc=="Data from CSV File" :
        st.header("General Instructions to read Data from a CSV File")
        """
- Create a CSV file named with different Columns where Data is stored under each (Make Sure Column Header is in Cell 1 of the CSV). (Figure 1) 
- Add  relevant data under each column header

"""
        image=Image.open('CSV_dataent_f1.png')
        st.image(image,width=500,caption='Figure 1')
        """
- Upload the File Below.

"""
        
        
        uploaded_file = st.file_uploader("Choose a CSV file", type="csv")

        if uploaded_file is not None:
          df = pd.read_csv(uploaded_file)
          st.write(df)

          obsdat1 = df['T'].values
          obsdat2= df['Drawdown'].values
      
          m_time=obsdat1
          m_ddown=obsdat2
          st.session_state.m_time = m_time
          st.session_state.m_ddown =  m_ddown

    tc_gen,curvematch,re_c = st.tabs(["Type Curve Genration", "Curve Matching","Result Computation"])
    with tc_gen :
        st.header(":green[Measured Data]")
        ###### Adding measured Data
        m_time = st.session_state.m_time
        m_ddown = st.session_state.m_ddown

        fig2 = plt.figure(figsize=(9,6))
        ax = fig2.add_subplot(1, 1, 1)
        plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9)  # adjust plot area
        ax.plot(m_time, m_ddown,'bo', markersize=3)
        plt.yscale("log")
        plt.xscale("log")
        plt.axis([0.1,1E4,1E-2,1E+1])
        ax.set(xlabel='time t (min)', ylabel='drawdown s (m)',
               title='Measured data')
        ax.grid(which="both", color='grey',linewidth=0.5)

        # Encoding the figure
        buffer2 = io.BytesIO()
        fig2.savefig(buffer2, format='png')
        buffer2.seek(0)
        data2 = base64.b64encode(buffer2.read()).decode()
        html2 = f'<img src="data:image/png;base64,{data2}" width=500, height=300>'
        st.write(html2, unsafe_allow_html=True)

        ##### Adding Theis Curve

        st.header(":green[Theis Curve Data]")
        menu2 = ["No","Yes"]
        choice2 = st.radio(":red[Do you want to modify the Theis Curve Axes]", menu2)
        r_max=10000
        u_max=100
        match_u_inv=1
        match_wu=1
        if choice2=="Yes" :
            u_max = st.number_input("U Max", value=10.0, max_value=100.0, step=0.1)
            r_max=st.number_input("R Max", value=10000,step=10000)
            match_u_inv=st.number_input("Scale Matching Curve 1/U", value=1,step=1)
            match_wu=st.number_input("Scale Matching Curve R Max", value=1,step=1)
        if choice2=="No" :
            r_max=10000
            u_max=100
            match_u_inv=1
            match_wu=1

        u  = [u_max for x in range(r_max)]
        um = [u_max for x in range(r_max)]
        u_inv  = [r_max/u_max for x in range(r_max)]
        um_inv = [r_max/u_max for x in range(r_max)]
        w_u  = [well_function(u_max/r_max) for x in range(r_max)]
        w_um = [well_function(u_max/r_max) for x in range(r_max)]

        for x in range(1,r_max,1):
            if x>0:
                u[x] = x*u_max/r_max
                u_inv[x] = 1/u[x]
                w_u[x] = well_function(u[x])
        matchgrid_x=[match_u_inv, match_u_inv]
        matchgrid_y=[match_wu, match_wu]
        matchgrid  =[0.001, 1000000]

        fig = plt.figure(figsize=(9,6))
        ax = fig.add_subplot(1, 1, 1)
        plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9)  # adjust plot area

        # Plot the data
        ax.plot(u_inv, w_u, color='black', linewidth=2)
        ax.plot(match_u_inv, match_wu, 'ro', markersize=6)
        ax.plot(matchgrid_x, matchgrid, color='lime', linewidth=1)
        ax.plot(matchgrid, matchgrid_y, color='lime', linewidth=1)

        # Set up the diagram
        plt.yscale('log')
        plt.xscale('log')
        plt.axis([0.1, 1E4, 1E-2, 1E+1])
        ax.set(xlabel='1/u', ylabel='w(u)', title='Theis type curve for manual evaluation')
        ax.grid(which='both', color='whitesmoke', linewidth=0.5)

        # Box around plot
        ax.spines['top'].set_color('lightgrey')
        ax.spines['bottom'].set_color('lightgrey')
        ax.spines['right'].set_color('lightgrey')
        ax.spines['left'].set_color('lightgrey')

        ax.tick_params(which='both', colors='lightgrey')

        # Display the plot on Streamlit
        #st.pyplot(fig)
        buffer1 = io.BytesIO()
        fig.savefig(buffer1, format='png')
        buffer1.seek(0)
        data1 = base64.b64encode(buffer1.read()).decode()
        html = f'<img src="data:image/png;base64,{data1}" width=500, height=300>'
        st.write(html, unsafe_allow_html=True)



    with curvematch :

        # Display the images with adjusted opacity and position using CSS
        opacity = st.number_input("Opacity",value=0.6)
        x_position = st.slider("X Position", -1000, 1000, 0)
        y_position = st.slider("Y Position", -1000, 1000, 0)

        st.markdown(
            f"""
            <div style="position: relative;">
                <img src="data:image/jpg;base64,{data2}" alt="Image 2" width="100%" style="display:block; opacity: {1-opacity};">
                <img src="data:image/png;base64,{data1}" alt="Image 1" width="100%" style="position: absolute; top: {y_position}px; left: {x_position}px; opacity: {opacity};">
            </div>
            """,
            unsafe_allow_html=True
        )
    with re_c :
        
        
        """
## Details of the computation
By fitting the type curve and the measured data and under consideration of the matching point, the values of t and s can be derived from the plots. Those parameters can be used to further compute transmissivity T and storativity S by the following equations
"""
        latex=(r'''$$ T=\frac{Q}{4\pi s}W(u) $$''')
        st.write(latex)   

        """
and subsequently
"""
        latex=(r'''$$ S = \frac{4Tt}{r^2}u $$''')
        st.write(latex)
        """

## Data Processing 
Under consideration of the pumping test data (extraction rate Q and distance r between pumping well and observation, the formation parameter can be computed. The next section is to check the manual parameter estimation. Estimated data (T and S) are set in a common plot. If the manual estimation was successful, measured data and type curve in the automatic plot will fit.**

**MAKE SURE THE UNITS ARE CONSISTENT (e.g. METERS / SECONDS)**
"""
        st.header(" :blue[Enter Data into the Sidebar]")
        r=st.sidebar.number_input("r -Distance between pumping well and observation", value=824*0.3048)
        Q=st.sidebar.number_input("Q -Pumping Rate", value=42400*0.3048**3/86400)
        t_est=st.sidebar.number_input("t -Estimated Transmisivity ", value=4.2*60)
        s_est=st.sidebar.number_input("s -Estimated Storativity ", value=0.81)

        T = Q * match_wu/(4*np.pi*s_est)
        S = 4*t_est*T/r**2/match_u_inv

        a=1
        #if st.button('Compute Result'):
        if a==1:
            choice="Result"
            T = Q * match_wu/(4*np.pi*s_est)
            S = 4*t_est*T/r**2/match_u_inv
            #st.write("Distance       r = ",'{:.2f}'.format(r), "m")
            #st.write("Pumpingrate    Q = ",'{:.2e}'.format(Q), "m3/s")
            st.write("Transmissivity T = ",'{:.2e}'.format(T), "m2/s")
            st.write("Storativity    S = ",'{:.2e}'.format(S), "[-]")
            x = 0
            for t in m_time:
                um[x] = theis_u(T,S,r,t*60)
                um_inv[x] = 1/um[x]
                w_um[x] = theis_wu(Q,T,m_ddown[x])
                x = x+1

            fig = plt.figure(figsize=(9,6))
            ax = fig.add_subplot(1, 1, 1)
            ax.plot(u_inv, w_u)
            ax.plot(um_inv, w_um,'ro')
            plt.yscale("log")
            plt.xscale("log")
            plt.axis([1,1E4,1E-2,1E+1])
            plot_text = "T= " + str('{:.2e}'.format(T))+ " m/s  &  S = " + str('{:.2e}'.format(S))+" [-]"
            plt.text(100, 0.1,plot_text,horizontalalignment='center',verticalalignment='top',multialignment='center', size='12')
            ax.set(xlabel='1/u', ylabel='w(u)',title='Theis drawdown')
            ax.grid(which="both",color='whitesmoke')
            plt.legend(('well function','measured'), loc=4)
            st.pyplot(fig)
        
                
        
menu = ["Introduction","Data Input"]
choice = st.sidebar.radio("Select a page", menu)

if choice == "Introduction" :
    intro()

if choice=="Data Input" :
    datainput()
    
