function E_g = calc_Eg(p,x);

% -- E_g = calc_Eg(p,x);
%
% The purpose of this function is to calculate the 
% expected value of the fitness gradient from the Taylor
% expansion for the stochastic model.
%
%
% INPUTS
%
% p: A dictionary of parameter values. Use the file params.m
% for an example of how to specify one of these.
%
% x: The breeding date.
%
%
% OUTPUTS
%
% E_g: The expected value of the fitness gradient.


% Derivatives of g, e.g. g_e means $g_e'$
[g,g_e,g_l,g_n,g_ee,g_el,g_en,g_ll,g_ln,g_nn] = calc_g(p,x);
[V_e,V_l,C_el,C_ve,C_vl,V_v,E_v] = calc_moments(p,x);

E_g = g + g_n*E_v + ...
    (1/2)*( g_nn*V_v + g_ee*V_e + g_ll*V_l ) + ...
    g_en*C_ve + g_ln*C_vl + g_el*C_el ;
