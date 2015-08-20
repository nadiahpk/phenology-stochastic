function [Pe,Pe_x,Pe_xn,Pe_xnn,Pe_n,Pe_nn] = calc_Pe(p,x,n);

% -- [Pe,Pe_x,Pe_xn,Pe_xnn,Pe_n,Pe_nn] = calc_Pe(p,x,n);
% -- [Pe,Pe_x,Pe_xn,Pe_xnn,Pe_n,Pe_nn] = calc_Pe(p,x);
%
% The purpose of this function is to calculate Pe and the
% various derivatives of it. Note that for Pe_x we use
% Pe(x',x,n) evaluated at x'=x, then all subsequent
% derivatives of Pe_x involve x (not x') only. However for Pe_n and
% its derivative, we use Pe(x,x,n) directly. 
%
% INPUTS
%
% p: A dictionary of parameter values. Use the file params.m
% for an example of how to specify one of these.
%
% x: The breeding date. 
%
% n: Deterministic steady-state population size. Optional.
%
%
% OUTPUTS
%
% Pe: Pe(x,x,n(x)) = K/n
%
% Pe_x: dPe(x',x,n(x))/dx' eval at x'=x
%
% Pe_xn: dPe_x/dn
%
% Pe_xnn: d^2Pe_x/dn^2
%
% Pe_n: dPe(x,x,n)/dn. Note that the Pe used here and below
% already have x'=x set, i.e. we are using Pe = K/n
%
% Pe_nn: d^2Pe/dn^2 using Pe = K/n

if nargin < 3
    n = calc_n(p,x);
end

m = p.m;
K = p.K;

Pe = K/n;
Pe_n = -K*n^(-2);
Pe_nn = 2*K*n^(-3);

% The following were found with Sage
Pe_x = -m.*e.^((x + log(-K./((K./n - 1).*n))./m).*m + m.*x)./(e.^((x + log(-K./((K./n - 1).*n))./m).*m) + e.^(m.*x)).^2;

Pe_xn = (K./n - 1).*(K./((K./n - 1).*n.^2) - K.^2./((K./n - 1).^2.*n.^3)).*m.*n.*e.^((x + log(-K./((K./n - 1).*n))./m).*m + m.*x)./((e.^((x + log(-K./((K./n - 1).*n))./m).*m) + e.^(m.*x)).^2.*K) - 2.*(K./n - 1).*(K./((K./n - 1).*n.^2) - K.^2./((K./n - 1).^2.*n.^3)).*m.*n.*e.^(2.*(x + log(-K./((K./n - 1).*n))./m).*m + m.*x)./((e.^((x + log(-K./((K./n - 1).*n))./m).*m) + e.^(m.*x)).^3.*K);

Pe_xnn = -(K./n - 1).^2.*(K./((K./n - 1).*n.^2) - K.^2./((K./n - 1).^2.*n.^3)).^2.*m.*n.^2.*e.^((x + log(-K./((K./n - 1).*n))./m).*m + m.*x)./((e.^((x + log(-K./((K./n - 1).*n))./m).*m) + e.^(m.*x)).^2.*K.^2) + 6.*(K./n - 1).^2.*(K./((K./n - 1).*n.^2) - K.^2./((K./n - 1).^2.*n.^3)).^2.*m.*n.^2.*e.^(2.*(x + log(-K./((K./n - 1).*n))./m).*m + m.*x)./((e.^((x + log(-K./((K./n - 1).*n))./m).*m) + e.^(m.*x)).^3.*K.^2) - 6.*(K./n - 1).^2.*(K./((K./n - 1).*n.^2) - K.^2./((K./n - 1).^2.*n.^3)).^2.*m.*n.^2.*e.^(3.*(x + log(-K./((K./n - 1).*n))./m).*m + m.*x)./((e.^((x + log(-K./((K./n - 1).*n))./m).*m) + e.^(m.*x)).^4.*K.^2) - 2.*(K./n - 1).*(K./((K./n - 1).*n.^3) - 2.*K.^2./((K./n - 1).^2.*n.^4) + K.^3./((K./n - 1).^3.*n.^5)).*m.*n.*e.^((x + log(-K./((K./n - 1).*n))./m).*m + m.*x)./((e.^((x + log(-K./((K./n - 1).*n))./m).*m) + e.^(m.*x)).^2.*K) + (K./n - 1).*(K./((K./n - 1).*n.^2) - K.^2./((K./n - 1).^2.*n.^3)).*m.*e.^((x + log(-K./((K./n - 1).*n))./m).*m + m.*x)./((e.^((x + log(-K./((K./n - 1).*n))./m).*m) + e.^(m.*x)).^2.*K) + 4.*(K./n - 1).*(K./((K./n - 1).*n.^3) - 2.*K.^2./((K./n - 1).^2.*n.^4) + K.^3./((K./n - 1).^3.*n.^5)).*m.*n.*e.^(2.*(x + log(-K./((K./n - 1).*n))./m).*m + m.*x)./((e.^((x + log(-K./((K./n - 1).*n))./m).*m) + e.^(m.*x)).^3.*K) - 2.*(K./n - 1).*(K./((K./n - 1).*n.^2) - K.^2./((K./n - 1).^2.*n.^3)).*m.*e.^(2.*(x + log(-K./((K./n - 1).*n))./m).*m + m.*x)./((e.^((x + log(-K./((K./n - 1).*n))./m).*m) + e.^(m.*x)).^3.*K) - (K./((K./n - 1).*n.^2) - K.^2./((K./n - 1).^2.*n.^3)).*m.*e.^((x + log(-K./((K./n - 1).*n))./m).*m + m.*x)./((e.^((x + log(-K./((K./n - 1).*n))./m).*m) + e.^(m.*x)).^2.*n) + 2.*(K./((K./n - 1).*n.^2) - K.^2./((K./n - 1).^2.*n.^3)).*m.*e.^(2.*(x + log(-K./((K./n - 1).*n))./m).*m + m.*x)./((e.^((x + log(-K./((K./n - 1).*n))./m).*m) + e.^(m.*x)).^3*n);

% For debugging
if 1 == 0
    plot(x_d,Pe,'k')
    hold on
    plot(x_d,Pe_x,'r')
    plot(x_d,Pe_xn,'g')
    plot(x_d,Pe_xnn,'b')
    hold off

end
if 1 == 0
    plot(n,Pe,'k')
    hold on
    plot(n,Pe_n,'r')
    plot(n,Pe_nn,'g')
    hold off
end
