import time
from scipy.optimize import newton
#from test import expected_stoch_fit_grad

import os

os.system('sage --preparse f1.sage')
os.system('mv f1.sage.py f1.py')
from f1 import expected_stoch_fit_grad, instant_det_grwth_rate, stoch_approximations_symbolic, det_evol_sing_strategy

os.system('sage --preparse default_det_pars.sage')
os.system('mv default_det_pars.sage.py default_det_pars.py')
from default_det_pars import pars

# -----

# Create our stochastic parameter values with our parameter of interest set as a variable
stoch_pars = {
    'G_epsilon': 0,
    'G_rho': 0,
    'C_omega_epsilon_omega_epsilon': 10,
    #'C_omega_epsilon_omega_epsilon': 0,
    'C_omega_rho_omega_rho':  100,
    #'C_omega_rho_omega_rho':  0,
    'C_omega_epsilon_omega_rho': 0,
} 

# There is a limit on the size of the correlation
par_name = 'C_omega_epsilon_omega_rho'
par_max = float(sqrt(stoch_pars['C_omega_epsilon_omega_epsilon'])*sqrt(stoch_pars['C_omega_rho_omega_rho']))
par_min = -par_max
par_n = 21

stoch_pars[par_name] = var(par_name)

# Other variables
state_vars = [var('n')]
stoch_vars = [var('epsilon'), var('rho')]

# Deterministic evolutionarily singular strategy
x_ss = det_evol_sing_strategy(pars)

w = instant_det_grwth_rate(pars, x_d = var('x'), x = var('x'))
n_ss = float(solve(w==1,n)[0].rhs().subs({x: x_ss}))

det_subs = {x: var('x_val'), n: n_ss, epsilon: 0, rho: 0} 
sa = stoch_approximations_symbolic(pars, stoch_pars, state_vars, stoch_vars, det_subs)

# This is a symbolic expression, with variables e.g.
# sage: E_g.variables() = (C_omega_rho_omega_rho, n_s, x_val)
E_g = expected_stoch_fit_grad(pars, stoch_pars, state_vars, stoch_vars, sa) 

# == Prepare output file and write preamble ==

# Prepare output file
fName = 'vary_' + par_name + '.dat'
f = open(fName, 'w')

# Preamble
f.write('# Created with varypar.sage ' + time.strftime("%Y-%m-%d %H:%M:%S") + '\n')
f.write('#\n')

# Write in the deterministic parameter values
f.write('# pars = {\n')
for k,v in pars.iteritems():
    f.write('#\t\'' + k + '\' : ' + str(v) + ',\n')
f.write('# }\n')
f.write('#\n')

# Write in the stochastic parameter values
f.write('# stoch_pars = {\n')
for k,v in stoch_pars.iteritems():
    f.write('#\t\'' + k + '\' : ' + str(v) + ',\n')
f.write('# }\n')
f.write('# The variable parameter is the one named above\n')
f.write('#\n')
f.write('#\n')

# Write in headers
f.write('# varypar \t x_stoch-x_ss \t\t E_v \n')

# equiv to Octave linspace(par_min,par_max,par_n)
parList = [ par_min + i*(par_max-par_min)/(par_n-1) for i in range(par_n) ]
for parVal in parList:

    # Substitute in our value of the stochastic parameter
    E_g_func(x_val) = E_g.subs({stoch_pars[par_name]: parVal, var('n_s'): n_ss}) # Expected fitness gradient

    # Solve for stochastic fitness steady state
    x_stoch = newton(lambda x_val: float(E_g_func(x_val)),x_ss) # Stochastic ESS trait value
    E_v_val = float(sa.E_v[0].subs({stoch_pars[par_name]: parVal, var('x_val'): x_stoch})) # Resulting deviation from n*

    # Write result to file for plotting
    fstr = '  {0:.4e} \t {1:.4e} \t {2:.4e}\n'.format(float(parVal), x_stoch-x_ss, E_v_val)
    f.write(fstr)

f.close()

# [x] NEXT: set some sensible values for the default parameter values, 
#   Find them and their effects in default_pars.dat
# [ ] then take a look at how the stochastic peak varies with the stochastic parameter values
