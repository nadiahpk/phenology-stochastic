function [V_e,V_l,C_el,C_ve,C_vl,V_v,E_v] = calc_moments(p,x);

% -- [V_e,V_l,C_el,C_ve,C_vl,V_v,E_v] = calc_moments(p,x);
%
% The purpose of this function is to calculate the moments
% of the environmental variables and the deviation from the
% deterministic population size. These are then used to
% approximate the expected value of the fitness gradient.
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
% V_e: Variance of environmental stochastic variable \epsilon
%
% V_l: Variance of environmental stochastic variable \lambda
%
% C_el: Covariance of environmental stochastic variables
% \epsilon and \lambda
%
% C_ve: Covariance of deviation of population size from n^*
% with environmental stochastic variable \epsilon 
%
% C_vl: Covariance of deviation of population size from n^*
% with environmental stochastic variable \lambda 
%
% V_v: Variance of deviation of population size from n^*
%
% E_v: Expected value of deviation of population size from n^*
%

% Autocorrelations
G_e = p.G_e;
G_l = p.G_l;

% Variances of underlying stochastic processes
V_oe = p.V_oe;
V_ol = p.V_ol;
C_oeol = p.C_oeol;

% Calculate moments for environmental process
V_e = V_oe/(1-G_e^2);
V_l = V_ol/(1-G_l^2);
C_el = C_oeol/(1-G_e*G_l);

% For the rest we'll need the derivatives of F
[F,F_e,F_l,F_n,F_ee,F_el,F_en,F_ll,F_ln,F_nn] = calc_F(p,x);

C_ve = G_e*(F_e*V_e + F_l*C_el)/(1-G_e*F_n);
C_vl = G_l*(F_l*V_l + F_e*C_el)/(1-G_l*F_n);
V_v = (1/(1-F_n^2))*( ...
    F_n*(F_e*C_ve + F_l*C_vl) + ...
    F_e*(F_n*C_ve + F_e*V_e + F_l*C_el) + ...
    F_l*(F_n*C_vl + F_e*C_el + F_l*V_l) );

E_v = (1/(1-F_n))*( ...
    .5*(F_nn*V_v + F_ee*V_e + F_ll*V_l) + ...
    F_en*C_ve + F_ln*C_vl + F_el*C_el );
