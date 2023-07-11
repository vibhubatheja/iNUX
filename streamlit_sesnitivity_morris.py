import streamlit as st
import numpy as np
import scipy as sp
from scipy import special
from scipy import integrate
import matplotlib.pyplot as plt
import openturns as ot
import time
import numpy as np
from SALib.sample import morris as salib_morris
from SALib.analyze import morris as salib_analyze

# Define your function for the arbitrary calculation
def calculate_output(x):
    Co_t, H_t, W_t, v_t, alx_t, aly_t, alz_t, R_t, SD_g_t, Deg_C_t = x
    lmax_transient = 0
    t_transient = 0
    # ... your existing code ...
    ## Input Parameters
    Cthres = 5e-5               #[mg/l]
    # Geometry - centreline
    y = 0                       #[m]                                                fixed (no Input)
    z_1 = 0                     #[m]                                                fixed (no input)
    z_2 = H_t                     #[m]                                               Input slider (lower limit: >0, upper limit:50)
    # Source term
    Co = Co_t    
    z = (z_1+z_2)/2                     #[m]
    W = W_t                      #[m]                                                Input slider (lower limit: >0, upper limit:1000)
    
    # hydraulic & mixing
    v = v_t                      #[m/y]                                              Input Box (lower limit: 10, upper limit: 1000)
    al_x = alx_t                   #[m]                                                Input slider (lower limit: 1, upper limit: 100, default: 10)
    al_y = aly_t                    #[m]                                                Input Box (lower limit: 0.1, upper limit: 10, default: 0.5)
    al_z = alz_t                 #[m]                                                Input Box (lower limit: 0.01, upper limit: 1, default: 0.05)
    Df = 0                      #[m^2/y]                                            Input Box (lower limit: 0 upper limit: 0.1, default: 0)
    
    # reaction terms
    R = R_t                       #[-]                                                Input Box (lower limit: 1, upper limit:100, default: 1)
    ga = SD_g_t                   #[1/y]                                              Input Box (lower limit: 0, upper limit: 1, default: 0)
    la = Deg_C_t                 #[1/y]                                              Input slider (lower limit: 0, upper limit: 1, default: 0.1)
    
    # Gauss points: max 256
    m = 256                     #[-]                                                 Input Box (possible values: 4,5,6,10,15,20,60,104,256; default: 60)
    
    Dx = al_x*v + Df            #[m^2/y]
    Dy = al_y*v + Df            #[m^2/y]
    Dz = al_z*v + Df            #[m^2/y]
    
    # used data
    vr = v/R                    #[m/y]
    Dyr = Dy/R                  #[m^2/y]
    Dxr = Dx/R                  #[m^2/y]
    Dyr = Dy/R                  #[m^2/y]
    Dzr = Dz/R                  #[m^2/y]
    def C_transient(x,t):
        
        # Boundary Condition
        
        if x<=1e-6:
            if y <= W/2 and y >= -W/2 and z <= z_2 and z >= z_1:
                C=Co*np.exp(-ga*t)
            else:
                C=0
        else:
            a = Co*np.exp(-ga*t)*x/(8*np.sqrt(np.pi*Dxr))
            roots = sp.special.roots_legendre(m)[0]
            weights = sp.special.roots_legendre(m)[1]
    
            #scaling
        
            bot = 0
            top = np.sqrt(np.sqrt(t))
            Tau = (roots*(top-bot)+top+bot)/2
            Tau4 = Tau**4
    
            #calculation
                
            xTerm = (np.exp(-(((la-ga)*Tau4)+((x-vr*Tau4)**2)/(4*Dxr*Tau4))))/(Tau**3)
            yTerm = sp.special.erfc((y-W/2)/(2*np.sqrt(Dyr*Tau4))) - sp.special.erfc((y+W/2)/(2*np.sqrt(Dyr*Tau4)))
            zTerm = sp.special.erfc((z-z_2)/(2*np.sqrt(Dzr*Tau4))) - sp.special.erfc((z-z_1)/(2*np.sqrt(Dzr*Tau4)))
            Term = xTerm * yTerm * zTerm
            Integrand = Term*(weights*(top-bot)/2)
            C = a*4*sum(Integrand)
        return C
    
    #function for Lmax depending on time past
    
    def Lmax(t):
        #x_array = np.array([0])                                                     #plotting
        #c_array = np.array([C_transient(0,t)])                                                #plotting
        x_min = 0  # minimum plume length to avoid unnecessary calculation (NOT a fixed value!)
        x_lower = x_min
        x_upper = 10000

        optimal_x = 0

        # Optimization Scheme 1
        while x_upper - x_lower > 1:
            #print(optimal_x)
            x_mid = (x_lower + x_upper) / 2
            C_mid = C_transient(x_mid, t)
            if C_mid >= Cthres:

                optimal_x = x_mid
                x_lower = x_mid
            else:
                x_upper = x_mid
        #print((C_transient(x_upper, t),Cthres,optimal_x,x_upper))
        if (optimal_x < x_upper) and (C_transient(x_upper, t) >= Cthres) :
            #print(C_transient(x_upper, t),Cthres)
            optimal_x = x_upper

        x = optimal_x
        return x
        """
        while C_transient(x,t) >= Cthres and x<100000:
            x = x+1
            x_array = np.append(x_array,x)                                          #plotting
            c_array = np.append(c_array,C_transient(x,t))                                     #plotting
            
        return x
        """
        
    
    t_min = 0                                                                       #starting time for numerical evaluation to avoid unnecessary calculation (NOT a fixed value!)
    t = t_min
    delta_t = 1                                                                     #delta t for defining steady state
    Lmax1 = Lmax(t)                                                                 #must be big otherwise maximum of function Lmax(x) can stop the loop
    Lmax2 = Lmax(t+delta_t)
    Lmax_array = np.array([0])                                                      #array for plotting
    t_array = np.array([0])                                                         #array for plotting
    
    while Lmax2 > Lmax1:                                                            #delta Lmax > 0 -> Criterium for steady state
        t = t+1                                                                     #resolution of time steps
        Lmax1 = Lmax(t)
        Lmax2 = Lmax(t+delta_t)
        Lmax_array = np.append(Lmax_array, Lmax(t))                                 #appending data to array (plotting)
        t_array = np.append(t_array, t)                                             #appending data to array (plotting)
        
    lmax_transient=(Lmax(t))
    t_transient=(t)

    #print(lmax_transient)
    #output = x1 + 2 * x2 + 0.5 * x3**2 + np.sin(x4) + np.exp(-x5)
    return lmax_transient

# Function to run the sensitivity analysis
def run_sensitivity_analysis():
    # Define the bounds for each variable
    problem_bounds = [
        [1, 120],       # Variable Co_t bounds
        [1, 16],        # Variable H_t bounds
        [1, 305],       # Variable W_t bounds
        [1, 60],        # Variable v_t bounds
        [1.5, 30.5],    # Variable alx_t bounds
        [0.5, 10.2],    # Variable aly_t bounds
        [0.0075, 1.5],  # Variable alz_t bounds
        [1, 2],         # Variable R_t bounds
        [0, 1],         # Variable SD_g bounds
        [0.1, 0.45]     # Variable Deg_C bounds
    ]

    # Set the number of trajectories and sample size
    r = 10
    sample_size = 1000

    # Perform the Morris sensitivity analysis
    problem = {
        'num_vars': len(problem_bounds),
        'names': ["Co_t", " H_t", " W_t", " v_t", "alx_t", " aly_t", " alz_t", " R_t", " SD_g_t", "Deg_C_t"],
        'bounds': problem_bounds,
    }

    # Generate the sample using SALib
    sample = salib_morris.sample(problem, sample_size, num_levels=r)

    # Evaluate the function at the sample points
    Y = []
    for i, x in enumerate(sample):
        Y.append(calculate_output(x))
        progress = (i + 1) / sample_size * 100
        st.progress(progress)

    # Convert Y to a NumPy array
    Y = np.array(Y)

    # Perform the Morris sensitivity analysis using SALib
    sensitivity_analysis = salib_analyze.analyze(problem, sample, Y, conf_level=0.95, num_levels=r)

    # Print the elementary effects
    for i in range(len(problem_bounds)):
        st.write(f"Elementary Effects for Variable {i+1}:")
        st.write(sensitivity_analysis['mu_star'][i])


# Streamlit app code
st.title("Sensitivity Analysis")

# Run the sensitivity analysis when the button is clicked
if st.button("Run Sensitivity Analysis"):
    with st.spinner("Running sensitivity analysis..."):
        run_sensitivity_analysis()
    st.success("Sensitivity analysis completed.")

