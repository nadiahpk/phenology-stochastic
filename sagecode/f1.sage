from collections import namedtuple

# Fitness components
# ------------------

def Pr_factory(pars, x_d=var('x_d'), rho=var('rho')):
    # Probability of recruitment

    # Assign needed parameter values
    h_r = pars['h_r']
    u_r = pars['u_r']

    # Find P_r
    P_r = 1/(exp(h_r*(x_d-(u_r+rho)))+1)

    return P_r


def Ps_factory(pars, x_d=var('x_d'), epsilon=var('epsilon')):
    # Probability of adult early-season survival

    # Assign needed parameter values
    h_s = pars['h_s']
    u_s = pars['u_s']

    # Find P_s
    P_s = 1/(exp(h_s*(x_d-(u_s+epsilon)))+1)

    return P_s


def Pe_factory(pars, x_d=var('x_d'), x=var('x'), n=var('n')):
    # Probability of territory acquisition

    # Assign needed parameter values
    m = pars['m']
    K = pars['K']

    # Find P_e
    if x == x_d:
        P_e = K/n # Just comes out cleaner this way
    else:
        t = (1/m)*log((K/n)/(1-(K/n))) # Half-probability point
        P_e = exp(m*(x+t))/(exp(m*(x+t))+exp(m*x_d))
    end

    return P_e


# Deterministic growth rates, fitness, and fitness gradients
# ----------------------------------------------------------

def instant_det_grwth_rate(pars, x_d=var('x_d'), x=var('x'), n=var('n')):
    # Instantaneous deterministic growth rate, Equation 6 14-12-31

    # Assign needed parameter values
    s = pars['s']
    a = pars['a']

    # Get each fitness component, full specification
    P_r = Pr_factory(pars, x_d=x_d, rho=0)
    P_s = Ps_factory(pars, x_d=x_d, epsilon=0)
    P_e = Pe_factory(pars, x_d, x, n)

    # Instantaneous deterministic fitness gradient
    w = s*(a*P_s*P_r*P_e + P_s)

    return w


def det_fit(pars):
    # Deterministic fitness (assuming separation of timescales), Equation 8 14-12-31
    #
    # You can do stuff like this with it
    #   w(x_d,x) = det_fit(pars)
    #
    # You will find that the fitness = 1 (approx) whenever x'=x, e.g. 
    #   d = float(w.subs({x_d:x}).subs({x:-2}))

    # First, get the expression for n*(x), we'll call it n_s 
    w = instant_det_grwth_rate(pars, x_d = var('x'), x = var('x'))
    n_s = solve(w==1,n)[0].rhs()

    # Now get the deterministic fitness with n = n*, so we have w(x',x)
    w = instant_det_grwth_rate(pars, x_d = var('x_d'), x = var('x'), n = n_s)

    return w


def det_fit_grad(pars, w=None):
    # Deterministic fitness gradient (assuming separation of timescales), Equation 9 14-12-31
    #
    # You can use it to evaluate the fitness gradient at a particular point like this
    #   w_x = det_fit_grad(pars)
    #   c = w_x.subs({x:-1}) 
    #
    # You can also use it to find the deterministic evolutionarily singular strategy like this
    #   w_x = det_fit_grad(pars)
    #   x_s = find_root(w_x==0,-1,-4)
    #
    # Or another way to get the deterministic evolutionarily singular strategy
    #   w_x(x) = det_fit_grad(pars)
    #   sp.optimize.newton(lambda x: float(w_x(x)),-3)

    if w is None:
        # Returns w(x',x)
        w = det_fit(pars)

    # Fitness gradient
    w_x = diff(w,x_d)
    w_x = w_x.subs({x_d: var('x')}) # Evaluated at x'=x

    return w_x

def det_evol_sing_strategy(pars):
    # A solver for the deterministic evolutionarily singular strategy, solves Equation 9 14-12-31

    # Assign needed parameter values
    u_s = pars['u_s']
    u_r = pars['u_r']

    # Get the fitness gradient as a function
    w_x(x) = det_fit_grad(pars)

    # A reasonable start value is midway between the P_r and P_s curves
    # TODO: Can we get a better guess for this?
    # x_0 = (u_s + u_r)/2 - 1
    x_0 = u_s

    # Using the scipy secant method
    x_s = newton(lambda x: float(w_x(x)),x_0)

    # Note that the alternative, in which I provide the second derivative, 
    # is substantially slower
    #   x_s = newton(lambda x: float(w_x(x)), -3, lambda x: float(w_xx(x)))

    return x_s

# Stochastic growth rates, fitness, and fitness gradients
# --------------------------------------------------------

def instant_stoch_grwth_rate(pars, x_d=var('x_d'), x=var('x'), n=var('n'), epsilon=var('epsilon'), rho=var('rho')):
    # Instantaneous stochastic growth rate, Equation 10 14-12-31
    #
    # It can be used like this
    # w = instant_stoch_grwth_rate(pars)

    # Assign needed parameter values
    s = pars['s']
    a = pars['a']

    # Get each fitness component, full specification
    P_r = Pr_factory(pars, x_d=x_d, rho=rho)
    P_s = Ps_factory(pars, x_d=x_d, epsilon=epsilon)
    P_e = Pe_factory(pars, x_d, x, n)

    # Instantaneous fitness gradient
    w = s*(a*P_s*P_r*P_e + P_s)

    return w

def stoch_fit_grad(pars):
    # Stochastic fitness gradient (assuming separation of timescales), Equation 12 14-12-31
    #
    # You can evaluate the fitness gradient at a point like this
    #   g(x,n,epsilon,rho) = stoch_fit_grad(pars)
    #   c = float(g(-1, 382.73, 0, 0))

    # Instantaneous stochastic growth rate, Equation 10 14-12-31
    w = instant_stoch_grwth_rate(pars)

    # The stochastic evolutionarily singular strategy maximises the geometric mean 
    # fitness or the expected value of the logarithm of the growth rate, which we call
    # here f. See Equation 11 14-12-31
    f = log(w)

    # The fitness gradient, Equation 12 14-12-31
    g = diff(f,x_d).subs({x_d:x})

    return g

def stoch_res_grwth_rate(pars, w=None):
    # Resident population growth rate, Equation 25 14-12-31

    if w is None:
        # Assign needed parameter values
        s = pars['s']
        a = pars['a']

        # Get each fitness component, full specification
        P_r = Pr_factory(pars, x_d=x, rho=rho)
        P_s = Ps_factory(pars, x_d=x, epsilon=epsilon)
        P_e = Pe_factory(pars, x, x, n)

        F = n*s*(a*P_s*P_r*P_e + P_s)
    else:
        F = n*w

    return F

def stoch_approximations_float(pars, stoch_pars, state_vars, stoch_vars, det_subs):
    # Returns an approximation of the expected deviation from the deterministic steady-state E_v (Equation 24 14-12-31 and Section 2.3 'Methods'), and an approximation of the variance-covariance matrix between the state and stochastic variables V_z (Section 4.3.2 14-12-31 and Equations 24-26 'Methods')
    # TODO: Haven't checked yet if it handles more general cases than our question

    # Create matrix Gamma
    Gamma = diagonal_matrix(list(stoch_pars['G_'+str(i)] for i in stoch_vars))

    # Create matrix Sigma
    Sigma = matrix(len(stoch_vars),list(stoch_pars['C_omega_' + str(stoch_vars[min(i,j)]) + '_omega_' + str(stoch_vars[max(i,j)])] for i in range(0,len(stoch_vars)) for j in range(0,len(stoch_vars))))

    # Variance-covariance matrix of stochastic variables
    V_stoch = (identity_matrix(len(stoch_vars)^2) - Gamma.tensor_product(Gamma,subdivide=False)) \vector(Sigma.T.list())


    # TODO: should we make it so this can handle a vector?
    F = stoch_res_grwth_rate(pars) # For our case, we only have one F, F is not a vector

    # First derivatives of F 
    F_stoch = matrix(list(float(diff(F,i).subs(det_subs)) for i in stoch_vars)).T
    F_state = matrix(list(float(diff(F,i).subs(det_subs)) for i in state_vars)).T

    # Covariances between state and stochastic variables
    C_state_stoch = (identity_matrix(len(stoch_vars)) - Gamma.tensor_product(F_state,subdivide=False)) \ (Gamma.tensor_product(F_stoch, subdivide=False).T*V_stoch)

    # Variance-covariance matrix of state variables. There's only one element, but let's make it a bit more general than needed so it can be extended more easily later.
    V_state = (identity_matrix(len(state_vars)^2) - F_state.tensor_product(F_state,subdivide=False)) \ (F_state.tensor_product(F_stoch,subdivide=False).T*C_state_stoch + \
        F_stoch.tensor_product(F_state,subdivide=False).T*C_state_stoch + \
        F_stoch.tensor_product(F_stoch,subdivide=False).T*V_stoch )

    # Second derivatives of F
    F_state_state = matrix(list(float(diff(F,i,j).subs(det_subs)) for i in state_vars for j in state_vars)).T
    F_state_stoch = matrix(list(float(diff(F,i,j).subs(det_subs)) for i in state_vars for j in stoch_vars)).T
    F_stoch_stoch = matrix(len(stoch_vars),list(float(diff(F,i,j).subs(det_subs)) for i in stoch_vars for j in stoch_vars)).T

    # Expected deviation from n_s
    E_v = (identity_matrix(len(state_vars)) - F_state).inverse() * (
        (1/2) * F_state_state.T * V_state + \
        F_state_stoch.T * C_state_stoch + \
        (1/2) * matrix(F_stoch_stoch.T.list()) * V_stoch
        )

    V_z = block_matrix([
        [ matrix(V_state), matrix(C_state_stoch) ], \
        [ matrix(list(C_state_stoch),nrows=2), matrix(list(V_stoch),nrows=2) ] \
                    ],subdivide=False)

    # Put into named tuple for return
    stoch_approxes = namedtuple('Stochastic_approximations', ['E_v','V_z'])

    return stoch_approxes(E_v,V_z)

def stoch_approximations_symbolic(pars, stoch_pars, state_vars, stoch_vars, det_subs):
    # Returns an approximation of the expected deviation from the deterministic steady-state E_v (Equation 24 14-12-31 and Section 2.3 'Methods'), and an approximation of the variance-covariance matrix between the state and stochastic variables V_z (Section 4.3.2 14-12-31 and Equations 24-26 'Methods')
    # TODO: Haven't checked yet if it handles more general cases than our question

    # Create matrix Gamma
    Gamma = diagonal_matrix(SR, list(stoch_pars['G_'+str(i)] for i in stoch_vars))

    # Create matrix Sigma
    Sigma = matrix(SR, len(stoch_vars), list(stoch_pars['C_omega_' + str(stoch_vars[min(i,j)]) + '_omega_' + str(stoch_vars[max(i,j)])] for i in range(0,len(stoch_vars)) for j in range(0,len(stoch_vars))))

    # Variance-covariance matrix of stochastic variables
    V_stoch = (identity_matrix(SR, len(stoch_vars)^2) - Gamma.tensor_product(Gamma,subdivide=False)) \vector(Sigma.T.list())

    # TODO: should we make it so this can handle a vector?
    F = stoch_res_grwth_rate(pars) # For our case, we only have one F, F is not a vector

    # First derivatives of F 
    F_stoch = matrix(SR, list(diff(F,i).subs(det_subs) for i in stoch_vars)).T
    F_state = matrix(SR, list(diff(F,i).subs(det_subs) for i in state_vars)).T

    # Covariances between state and stochastic variables
    C_state_stoch = (identity_matrix(SR, len(stoch_vars)) - Gamma.tensor_product(F_state,subdivide=False)) \ (Gamma.tensor_product(F_stoch, subdivide=False).T*V_stoch)

    # Variance-covariance matrix of state variables. There's only one element, but let's make it a bit more general than needed so it can be extended more easily later.
    V_state = (identity_matrix(SR, len(state_vars)^2) - F_state.tensor_product(F_state,subdivide=False)) \ (F_state.tensor_product(F_stoch,subdivide=False).T*C_state_stoch + \
        F_stoch.tensor_product(F_state,subdivide=False).T*C_state_stoch + \
        F_stoch.tensor_product(F_stoch,subdivide=False).T*V_stoch )

    # Second derivatives of F
    F_state_state = matrix(SR, list(diff(F,i,j).subs(det_subs) for i in state_vars for j in state_vars)).T
    F_state_stoch = matrix(SR, list(diff(F,i,j).subs(det_subs) for i in state_vars for j in stoch_vars)).T
    F_stoch_stoch = matrix(SR, len(stoch_vars),list(diff(F,i,j).subs(det_subs) for i in stoch_vars for j in stoch_vars)).T

    # Expected deviation from n_s
    E_v = (identity_matrix(SR, len(state_vars)) - F_state).inverse() * (
        (1/2) * F_state_state.T * V_state + \
        F_state_stoch.T * C_state_stoch + \
        (1/2) * matrix(F_stoch_stoch.T.list()) * V_stoch
        )

    V_z = block_matrix([
        [ matrix(V_state), matrix(C_state_stoch) ], \
        [ matrix(list(C_state_stoch),nrows=2), matrix(list(V_stoch),nrows=2) ] \
                    ],subdivide=False)

    # Put into named tuple for return
    stoch_approxes = namedtuple('Stochastic_approximations', ['E_v','V_z'])

    return stoch_approxes(E_v,V_z)


def expected_stoch_fit_grad(pars, stoch_pars, state_vars, stoch_vars, sa = None):

    if sa == None:

        # Get the linearised approximations of the stochastic aspects E_v, V_z
        det_subs = {x: var('x_val'), n: var('n_s'), epsilon: 0, rho: 0} 
        sa = stoch_approximations_symbolic(pars, stoch_pars, state_vars, stoch_vars, det_subs)

    # Find the expected value of the stochastic fitness
    all_vars = state_vars + stoch_vars

    # The second derivatives of g to prepare for estimating the expected fitness value
    stoch_subs = {n: var('n_s'), epsilon: 0, rho: 0}
    g = stoch_fit_grad(pars)
    g_derivs = matrix(list((diff(g,i,j).subs(stoch_subs)) for i in all_vars for j in all_vars),nrows=len(all_vars))

    E_g = (g.subs(stoch_subs) + diff(g,n).subs(stoch_subs)*sa.E_v[0] + (1/2)*sum(sum(g_derivs.elementwise_product(sa.V_z)))).subs({x: var('x_val')}) 

    return E_g

def expected_stoch_fit_grad_float(x_val, pars, stoch_pars, state_vars, stoch_vars):
    # Expected value of the fitness gradient

    # Solve deterministic population dynamics
    w = instant_det_grwth_rate(pars, x_d = var('x'), x = var('x'))
    n_s = float(solve(w==1,n)[0].rhs().subs({x: x_val}))
    
    # Create a dictionary of substitutions, for where the
    # approximation is taken around the deterministic steady state
    det_subs = {x: x_val, n: n_s, epsilon: 0, rho: 0}

    # Get the linearised approximations of the stochastic aspects E_v, V_z
    sa = stoch_approximations_float(pars, stoch_pars, state_vars, stoch_vars, det_subs)

    # Find the expected value of the stochastic fitness
    all_vars = state_vars + stoch_vars

    # The second derivatives of g to prepare for estimating the expected fitness value
    stoch_subs = {n: n_s, epsilon: 0, rho: 0}
    g = stoch_fit_grad(pars)
    # slow ***
    g_derivs = matrix(list((diff(g,i,j).subs(stoch_subs)) for i in all_vars for j in all_vars),nrows=len(all_vars))

    # slow ***
    E_g = float( (g.subs(stoch_subs) + diff(g,n).subs(stoch_subs)*sa.E_v[0] + (1/2)*sum(sum(g_derivs.elementwise_product(sa.V_z)))).subs({x:x_val}) )

    return E_g

# ------------- Main -------------

# Stuff below

from scipy.optimize import newton

# TODO: see spare.sage for the stuff you've deleted

'''
# def go():
if True:

    # Parameter values
    pars = {
        'm': 1/2,
        'K': 100,
        'u_s': -5,
        'h_s': -1/2,
        'u_r':  5,
        'h_r':  1/2,
        'a': 3,
        's': 65/100
    }
    stoch_pars = {
        'G_epsilon': 1/10,
        'G_rho': 2/10,
        'C_omega_epsilon_omega_epsilon': 5,
        'C_omega_rho_omega_rho':  4,
        'C_omega_epsilon_omega_rho': 1,
    } # TODO: May need to check that order matches stoch_vars

    state_vars = [var('n')]
    stoch_vars = [var('epsilon'), var('rho')]

    # Find deterministic evolutionarily singular strategy
    #x_s = det_evol_sing_strategy(pars)
    x_s = -3.547657271113175 # save some time, this is what should be returned above

    # Note to self: for these parameters and at x_s, E_g
    # should equal 0.05489465879588512

    E_g = expected_stoch_fit_grad(pars, stoch_pars, state_vars, stoch_vars)

    w = instant_det_grwth_rate(pars, x_d = var('x'), x = var('x'))

    x_stoch = newton(lambda x_ss: float(E_g(x_val = x_ss, n_s = float(solve(w==1,n)[0].rhs().subs({x: x_ss})))), x_s)

    print x_stoch

    This below was what was here before 150415
    #x_stoch = newton(lambda x: expected_stoch_fit_grad_float(x, pars, stoch_pars, state_vars, stoch_vars), x_s)

    #x_stoch = newton(lambda x: float(E_g(x)),x_s)

    #E_g = expected_stoch_fit_grad_float(x_s, pars, stoch_pars, state_vars, stoch_vars)
    '''
