function [Pr,Pr_x,Pr_xl,Pr_xll,Pr_l,Pr_ll] = calc_Pr(p,x);

% -- [Pr,Pr_x,Pr_xl,Pr_xll,Pr_l,Pr_ll] = calc_Pr(p,x);
%
% Finds the recruitment probability Pr and various derivatives
% of it.
%
% INPUTS
%
% p: A dictionary of parameter values. Use the file params.m
% for an example of how to specify one of these.
%
% x: The breeding date. x may be a vector for debugging
%
% OUTPUTS
%
% Pr: Survival probability.
%
% Pr_x: \partial P_r(x',\lambda)/\partial x' evaluated at x'=x
%
% Pr_xe etc.: \partial Pr_x/\partial \lambda

u_r = p.u_r;
h_r = p.h_r;

lambda = 0; % The evaluation point

Pr = 1./(exp(h_r.*(x-(u_r+lambda)))+1);

Pr_x = -h_r.*exp((lambda+x-u_r).*h_r)./(exp((lambda+x-u_r).*h_r)+1).^2;

Pr_xl = h_r.^2.*e.^(-(lambda + u_r - x).*h_r)./(e.^(-(lambda + u_r - x).*h_r) + 1).^2 - 2.*h_r.^2.*e.^(-2.*(lambda + u_r - x).*h_r)./(e.^(-(lambda + u_r - x).*h_r) + 1).^3;

Pr_xll = -h_r.^3.*e.^(-(lambda + u_r - x).*h_r)./(e.^(-(lambda + u_r - x).*h_r) + 1).^2 + 6.*h_r.^3.*e.^(-2.*(lambda + u_r - x).*h_r)./(e.^(-(lambda + u_r - x).*h_r) + 1).^3 - 6.*h_r.^3.*e.^(-3.*(lambda + u_r - x).*h_r)./(e.^(-(lambda + u_r - x).*h_r) + 1).^4;

Pr_l = -Pr_x;

Pr_ll = -Pr_xl;

% For debugging
if 0 == 1
    plot(x,Pr,'k');
    hold on
    plot(x,Pr_x,'r');
    plot(x,Pr_xl,'g');
    plot(x,Pr_xll,'b');
    plot(x,Pr_l,'m');
    hold off
end

