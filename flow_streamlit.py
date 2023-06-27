import math
import numpy as np
import streamlit as st
import matplotlib.pyplot as plt

def calculate_flow_rates(delta_l2,T_w, V, d, kc,g,tau,delta_h):


    delta_l = np.arange(delta_l2, delta_l2+100, 0.1)
    rho = (999.83952 + 16.945176 * T_w - 7.9870401 * T_w**2 - 3 * T_w**3 - 46.170461 * T_w**4 - 6 * T_w**5 + 105.56302 - 9 * T_w**6 - 280.54253 - 12 * T_w**7) / (1 + 16.879850 - 3 * T_w)
    
    # Calculating dynamic viscosity (mu) for cases where T_w is less than 20
    if T_w < 20:
        mu = (10 * ((1301 - 1.30233) / (998.33 + 8.1855 * (T_w - 20) + 0.00585 * (T_w - 20)**2))) / 1000
    else:
        mu = 10**((1.3272 * (20 - T_w) - 0.001053 * (T_w - 20)**2) / (T_w + 105)) * (1.002 / 1000)

    Q_darcy = []
    grad = []
    Q_hp=[]
    Q_turbulent=[]
    
    for delta_l_val in delta_l:
        # Calculating Reynolds number (Re)
        Re = (rho * V * d) / mu

        # Solving for friction factor (f)
        f_guess = 0.02  # Initial guess for f
        f_tolerance = 1e-5  # Tolerance for f iteration
        f_max_iterations = 100  # Maximum number of iterations for f iteration

        for f_n in range(f_max_iterations):
            f_inv_sqrt = -2 * math.log10((kc / (3.71 * d)) + (2.51 / (Re * math.sqrt(f_guess))))
            f_new = 1 / (f_inv_sqrt**2)

            if abs(f_new - f_guess) < f_tolerance:
                f = f_new
                break
            else:
                f_guess = f_new

        # Calculating head loss (∆h or hL)
        #delta_h = f * (delta_l_val / d) * (V**2 / (2 * g))

        # Calculating Darcy flow rate (Q_darcy)
        Q_darcy_val = ((3.14 * d**2) / 4) * ((delta_h * d * 2 * g) / (f * delta_l_val))**0.5

        # Calculating the gradient (delta_h / delta_l)
        grad_val = delta_h / delta_l_val

        Q_darcy.append(Q_darcy_val)
        grad.append(grad_val)
        # Calculating flow rate using Hagen-Poiseuille equation (for laminar flow)
        Q_hp_val=(- (3.14 * d**4 * rho * g * delta_h) / (128 * mu * delta_l_val * tau))


        # Calculating turbulent flow rate using Horlacher and Ludecke (1992) equation
        Q_turbulent_val=(math.sqrt(abs(delta_h) * g * d**5 * math.pi**2 / (2 * delta_l_val * tau)) * math.log((2.51 * V) / math.sqrt((2 * abs(delta_h) * g * d**3) / (delta_l_val * tau)) + kc / (3.71 * d)) * delta_h / abs(delta_h))

        Q_hp.append(Q_hp_val)
        Q_turbulent.append(Q_turbulent_val)

    return Q_darcy,Q_hp,Q_turbulent, grad

# Streamlit app code
st.title("Flow Rate Calculation")

# Input parameters
delta_l2 = st.sidebar.number_input("Length of the pipe segment (∆l)", value=10.0)
d = st.sidebar.number_input("Pipe diameter (d)", value=5.0)
V = st.sidebar.number_input("Mean velocity (V)", value=20.0)
g=9.81
#g = st.number_input("Gravitational acceleration constant (g)", value=9.8)
tau = st.sidebar.number_input("Shear stress (τ)", value=1.0)
kc = st.sidebar.number_input("Equivalent sand roughness (kc)", value=1.0)
T_w = st.sidebar.number_input("Temperature (T_w)", value=15.0)
delta_h = st.sidebar.number_input("Delta_h", value=0.6)
# Calculate flow rates
Q_darcy,Q_hp,Q_turbulent,grad = calculate_flow_rates(delta_l2,T_w, V, d, kc,g,tau,delta_h)

# Plotting Q_darcy vs grad
plt.scatter(Q_darcy, grad, label='Q_darcy')
plt.plot(Q_turbulent, grad,color='red', label='Q_turbulent')
plt.xlabel('Flow')
plt.ylabel('gradent')
plt.title('Flow vs Gradient')
plt.legend()
st.pyplot(plt)


