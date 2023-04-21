import numpy as np
import math
import csv
import matplotlib.pyplot as plt
from ipywidgets import interact, interactive, fixed, interact_manual
import ipywidgets as widgets
import streamlit as st



def curvefitting(mod_factor, cd, test_time, norm_heads):
    
    cdts=[]
    dless_time = np.arange(0, 30.1, 0.1)
    omega=math.sqrt(abs(1 - (cd/2)**2)) 
    omegaplus=(-cd/2) + math.sqrt(abs(1 - (cd/2)**2))
    omegaminus=(-cd/2) - math.sqrt(abs(1 - (cd/2)**2))
    
    for j in dless_time :
        
        if cd > 2:
            numerator = omegaminus * math.exp(omegaplus * j) - omegaplus * math.exp(omegaminus * j)
            denominator = (omegaminus - omegaplus)**(-1)
            result = denominator * numerator
        elif cd == 2:
            result = math.exp(-j) * (1 + j)
        else:
            numerator = math.cos(omega * j) + (cd / (2 * omega)) * math.sin(omega * j)
            denominator = math.exp(-(cd * j) / 2)
            result = denominator * numerator
        cdts.append(result)  
    
    cdts = np.array(cdts)
    adjtimes = mod_factor * np.array(dless_time)
    
    # Create plot
    fig, ax = plt.subplots()
    ax.plot(adjtimes, cdts, label='Cd data')
    ax.plot(test_time, norm_heads, label='field data')
    ax.set_xlabel('time')
    ax.set_ylabel('head')
    ax.set_title('Curve Matching')
    ax.legend()
    
    # Display plot in Streamlit app
    st.pyplot(fig)
    
    return (mod_factor, cd)


static_level = 0 # in meters
og_level = 21.288 #in meters
H_0 = -0.288 #in meters
start_time = 28970 # in seconds
depth_to_bottom_screen = 93.06 #in meters from toc
b_usp = 4 # screen_length (b) in meters 
depth_static_water = 53.1 # in meters  from toc
d_usp = 35.960 #topscreen_to_watertab
rw_usp = 0.250 #Radius of well screen in meters
rnc_usp = 0.050 # Nominal radius of well casing
rtc_usp = 0.003 # Radius of Transducer Cable
rw2_usp = 0.05 # Modified screen radius 
B = 75 # formation thickness in m 


aqt = ["Confined","Unconfined"]
aqt2 = st.sidebar.radio("Select Aquifier Type", aqt)

# Create sidebar with input fields
st.sidebar.title("Input parameters")
static_level = st.sidebar.number_input("Static water level (m)", value=static_level)
og_level = st.sidebar.number_input("OG level (m)", value=og_level)
H_0 = st.sidebar.number_input("Reference elevation (m)", value=H_0)
start_time = st.sidebar.number_input("Start time (s)", value=start_time)
depth_to_bottom_screen = st.sidebar.number_input("Depth to bottom screen (m)", value=depth_to_bottom_screen)
b_usp = st.sidebar.number_input("Screen length (m)", value=b_usp)
depth_static_water = st.sidebar.number_input("Depth of static water level (m)", value=depth_static_water)
d_usp = st.sidebar.number_input("Top screen to water table (m)", value=d_usp)
rw_usp = st.sidebar.number_input("Radius of well screen (m)", value=rw_usp)
rnc_usp = st.sidebar.number_input("Nominal radius of well casing (m)", value=rnc_usp)
rtc_usp = st.sidebar.number_input("Radius of Transducer Cable (m)", value=rtc_usp)
rw2_usp = st.sidebar.number_input("Modified screen radius (m)", value=rw2_usp)
B = st.sidebar.number_input("Formation thickness (m)", value=B)


mod_factor=7  # Modulation factor for Cd time series
cd=25.3

rc=math.sqrt((rnc_usp**2)-(rtc_usp**2)) #Effective Radius of casing in m 
aspect_ratio=b_usp/(rw2_usp)

with open('feilddata.csv', 'r') as infile:
  # read the file as a dictionary for each row ({header : value})
  reader = csv.DictReader(infile)
  data = {}
  for row in reader:
    for header, value in row.items():
      try:
        data[header].append(value)
      except KeyError:
        data[header] = [value]
        


obsdat1= data['T']           # Comment : Need to Modify this in Excel make sure header name is same 
obsdat2= data['Pressure Head in Meters']
times = list(map(float,obsdat1))
presshs = list(map(float,obsdat2))
times= np.array(times)
presshs= np.array(presshs)

test_time=times-start_time
stat_dev=presshs-static_level
norm_heads=stat_dev/H_0

# Define the UI using Streamlit and ipywidgets
st.title("Curve Fitting Parameters")
mod_factor = st.slider("Modifying factor", min_value=0.0, max_value=100.0, step=0.1, value=7.0)
cd = st.slider("Cd", value=25.7,max_value=100.0,step=0.1)

# Call the curve fitting function
curvefitting(mod_factor, cd, test_time, norm_heads)


st.write(cd,mod_factor)

timeCR=1/mod_factor
le_rat= 9.81/(timeCR**2)
le_nom= depth_to_bottom_screen-(b_usp+depth_static_water)+((b_usp/2))*((rnc_usp**2)/(rw_usp**2))
le_diff=abs((le_rat-le_nom)/le_nom)


bracketterm= (aspect_ratio/2 + math.sqrt(1 + (aspect_ratio/2)**2))
Kr_conf_mps= (timeCR * rc**2 * math.log(bracketterm)) / (2 * b_usp * cd)
Kr_conf_mpd= Kr_conf_mps*86400
Kr_conf_cmps= Kr_conf_mps*100
Kr_conf_ftpd= Kr_conf_mpd*3.281


termA = 1.472 + 0.03537*aspect_ratio - 0.00008148*aspect_ratio**2 + 0.0000001028*aspect_ratio**3 - 0.00000000006484*aspect_ratio**4 + 0.00000000000001573*aspect_ratio**5
termB = 0.2372 + 0.005151*aspect_ratio - 0.000002682*aspect_ratio**2 - 0.0000000003491*aspect_ratio**3 + 0.0000000000004738*aspect_ratio**4
termC2= math.log((B - (d_usp + b_usp)) / rw2_usp)



if termC2 > 6:
    termC = 6
else:
    termC = math.log((B - (d_usp + b_usp)) / rw2_usp)


secondterm=(termA + termB * termC) / aspect_ratio

firstterm= 1.1 / math.log((d_usp + b_usp) / rw2_usp)


Rerw=(firstterm + secondterm)**(-1)

Kr_unconf_mps=(timeCR * rc ** 2 * Rerw) / (2 * b_usp * cd)
Kr_unconf_mpd= Kr_unconf_mps*86400
Kr_unconf_cmps= Kr_unconf_mps*100
Kr_unconf_ftpd= Kr_unconf_mpd*3.281

st.title("Results from the Input Data")

if aqt2=="Unconfined" :
    st.write("Kr for Unconfined in meters per second",Kr_unconf_mps)
if aqt2=="Confined" :
    st.write("Kr for Confined in meters per second",Kr_conf_mps)

