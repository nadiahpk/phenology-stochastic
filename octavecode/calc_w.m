function [w_x,w] = calc_w(p,x);
% Deterministic fitness gradient and fitness

s = p.s; a = p.a;

[Ps,Ps_x,Ps_xe,Ps_xee,Ps_e,Ps_ee] = calc_Ps(p,x);
[Pr,Pr_x,Pr_xl,Pr_xll,Pr_l,Pr_ll] = calc_Pr(p,x);
[Pe,Pe_x,Pe_xn,Pe_xnn,Pe_n,Pe_nn] = calc_Pe(p,x);

w = s*(a*Ps*Pe*Pr + Ps);

w_x = s*Ps_x + s*a*( ...
    Ps_x*Pe*Pr + ...
    Pe_x*Ps*Pr + ...
    Pr_x*Ps*Pe );

