function [g,g_e,g_l,g_n,g_ee,g_el,g_en,g_ll,g_ln,g_nn] = calc_g(p,x);

% -- [g,g_e,g_l,g_n,g_ee,g_el,g_en,g_ll,g_ln,g_nn] = calc_g(p,x);
%
% Finds various derivatives of the stochastic fitness
% gradient.
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
% g: Stochastic fitness gradient.
%
% g_e, etc.: Derivatives of g. 
% For example, g_e means $g_\eps' =
% dg(x,n,\eps,\lambda)/d\eps$ evaluated at $n=n^*$,
% $\eps=0$, $\lambda=0$.
%

s = p.s; a = p.a;

[Ps,Ps_x,Ps_xe,Ps_xee,Ps_e,Ps_ee] = calc_Ps(p,x);
[Pr,Pr_x,Pr_xl,Pr_xll,Pr_l,Pr_ll] = calc_Pr(p,x);
[Pe,Pe_x,Pe_xn,Pe_xnn,Pe_n,Pe_nn] = calc_Pe(p,x);

w = s*(a*Ps*Pe*Pr + Ps);

w_x = s*Ps_x + s*a*( ...
    Ps_x*Pe*Pr + ...
    Pe_x*Ps*Pr + ...
    Pr_x*Ps*Pe );

w_xe = s*Ps_xe + s*a*( ...
    Ps_xe*Pe*Pr + Ps_e*( ...
        Pe_x*Pr + Pr_x*Pe ));
w_xl = s*a*( ...
    Pr_l*( ...
        Ps_x*Pe + Pe_x*Ps ) +
    Pr_xl*Ps*Pe );
w_xn = s*a*( ...
    Pe_n*( ...
        Ps_x*Pr + Pr_x*Ps ) +
    Pe_xn*Ps*Pr );

w_xee = s*Ps_xee + s*a*( ...
    Ps_xee*Pe*Pr + Ps_ee*( ...
        Pe_x*Pr + Pr_x*Pe ));
w_xel = s*a*( ...
    Pr_l*( ...
        Ps_xe*Pe + Ps_e*Pe_x ) +
    Pr_xl*Pe*Ps_e );
w_xen = s*a*( ...
    Pe_n*( ...
        Ps_xe*Pr + Ps_e*Pr_x ) +
    Pe_xn*Ps_e*Pr );

w_xll = s*a*( ...
    Pr_ll*( ...
        Ps_x*Pe + Pe_x*Ps ) + 
    Pr_xll*Ps*Pe );
w_xln = s*a*( ...
    Pe_n*( ...
        Pr_l*Ps_x + Pr_xl*Ps ) + 
    Pe_xn*Pr_l*Ps );

w_xnn = s*a*( ...
    Pe_nn*( ...
        Ps_x*Pr + Pr_x*Ps ) + 
    Pe_xnn*Ps*Pr );

w_e = s*(a*Pe*Pr*Ps_e + Ps_e);
w_l = s*a*Pe*Pr_l*Ps;
w_n = s*a*Pr*Ps*Pe_n;

w_ee = s*(a*Pe*Pr*Ps_ee + Ps_ee);
w_el = s*a*Pe*Pr_l*Ps_e; 
w_en = s*a*Pe_n*Pr*Ps_e; 

w_ll = s*a*Pe*Pr_ll*Ps;
w_ln = s*a*Pe_n*Pr_l*Ps;

w_nn = s*a*Ps*Pr*Pe_nn;

w1 = 1/w; w2 = -1/(w^2); w3 = 2/(w^3);

w1_e = w2*w_e;
w1_l = w2*w_l;
w1_n = w2*w_n;

w2_e = w3*w_e;
w2_l = w3*w_l;
w2_n = w3*w_n;

g = w1*w_x;

g_e = w2*w_e*w_x + w1*w_xe;
g_l = w2*w_l*w_x + w1*w_xl;
g_n = w2*w_n*w_x + w1*w_xn;

g_ee = w2_e*w_e*w_x + w2*w_ee*w_x + w2*w_e*w_xe + ...
        w1_e*w_xe + w1*w_xee;
%g_el = w2_l*w_e*w_x + w2*w_el*w_x + w2*w_e*w_xl + ...
        %w1_l*w_xe + w1*w_xel;

%g_ab = w1_b*w_xa + w1*w_xab + w_xb*w2*w_a + ...
        %w_x*w2_b*w_a + w_x*w2*w_ab;
g_el = w1_l*w_xe + w1*w_xel + w_xl*w2*w_e + ...
        w_x*w2_l*w_e + w_x*w2*w_el;

%g_en = w2_n*w_e*w_x + w2*w_en*w_x + w2*w_e*w_xn + ...
        %w1_n*w_xe + w1*w_xen;
g_en = w1_n*w_xe + w1*w_xen + w_xn*w2*w_e + ...
        w_x*w2_n*w_e + w_x*w2*w_en;

g_ll = w2_l*w_l*w_x + w2*w_ll*w_x + w2*w_l*w_xl + ...
        w1_l*w_xl + w1*w_xll;
g_ln = w2_n*w_l*w_x + w2*w_ln*w_x + w2*w_l*w_xn + ...
        w1_n*w_xl + w1*w_xln;

g_nn = w2_n*w_n*w_x + w2*w_nn*w_x + w2*w_n*w_xn + ...
        w1_n*w_xn + w1*w_xnn;

