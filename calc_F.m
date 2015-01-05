function [F,F_e,F_l,F_n,F_ee,F_el,F_en,F_ll,F_ln,F_nn] = calc_F(p,x,n)

s = p.s; a = p.a;

if nargin < 3
    n = calc_n(p,x);
end

[Ps,Ps_x,Ps_xe,Ps_xee,Ps_e,Ps_ee] = calc_Ps(p,x);
[Pr,Pr_x,Pr_xl,Pr_xll,Pr_l,Pr_ll] = calc_Pr(p,x);
[Pe,Pe_x,Pe_xn,Pe_xnn,Pe_n,Pe_nn] = calc_Pe(p,x,n);

F = n*s*(a*Ps*Pr*Pe + Ps);

F_e = n*s*(a*Pr*Pe*Ps_e + Ps_e);
F_l = n*s*a*Ps*Pe*Pr_l;
F_n = s*(a*Ps*Pr*Pe + Ps) + n*s*a*Ps*Pr*Pe_n;

F_ee = n*s*(a*Pr*Pe*Ps_ee + Ps_ee);
F_el = n*s*a*Pr_l*Pe*Ps_e;
F_en = s*(a*Ps_e*Pr*Pe + Ps_e) + n*s*a*Ps_e*Pr*Pe_n;
F_ll = n*s*a*Ps*Pe*Pr_ll;
F_ln = s*a*Ps*Pr_l*Pe + n*s*a*Ps*Pr_l*Pe_n;
F_nn = s*a*Ps*Pr*Pe_n + s*a*Ps*Pr*Pe_n + n*s*a*Ps*Pr*Pe_nn;


