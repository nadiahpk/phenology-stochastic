function n = calc_n(p,x);

Ps = calc_Ps(p,x);
Pr = calc_Pr(p,x);

n = p.s*p.a*p.K*Ps*Pr/(1-p.s*Ps);
