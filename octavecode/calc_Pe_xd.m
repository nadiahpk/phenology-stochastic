function [Pe,Pe_x] = calc_Pe_xd(p,x,x_d,n);

% Calculates value and derivatives of Pe(x',x,n)
% Returns
% Pe = Pe(x',x,n(x))
% dPe(x',x,n(x))/dx' eval at x'=x

if nargin < 4;
    n = calc_n(p,x);
    if nargin < 3;
        x_d = x;
    end
end

m = p.m;
K = p.K;

t = (1/m)*log((K./n)./(1-(K./n))); % Half-probability point
Pe = exp(m*(x+t))./(exp(m*(x+t))+exp(m*x_d));

Pe_x = -m.*e.^((x + log(-K./((K./n - 1).*n))./m).*m + m.*x_d)./(e.^((x + log(-K./((K./n - 1).*n))./m).*m) + e.^(m.*x_d)).^2;


% For debugging
if 1 == 0
    plot(x_d,Pe,'k')
    hold on
    plot(x_d,Pe_x,'r')
    hold off

end
