function [Ps,Ps_x,Ps_xe,Ps_xee,Ps_e,Ps_ee] = calc_Ps(p,x);

% -- [Ps,Ps_x,Ps_xe,Ps_xee,Ps_e,Ps_ee] = calc_Ps(p,x);
%
% Finds the survival probability Ps and various derivatives
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
% Ps: Survival probability.
%
% Ps_x: \partial P_s(x',\eps)/\partial x' evaluated at x'=x
%
% Ps_xe etc.: \partial Ps_x/\partial \eps

u_s = p.u_s;
h_s = p.h_s;

epsilon = 0; % The evaluation point

Ps = 1./(exp(h_s.*(x-(u_s+epsilon)))+1);

Ps_x = -h_s.*exp((epsilon+x-u_s).*h_s)./(exp((epsilon+x-u_s).*h_s)+1).^2;

Ps_xe = h_s.^2.*e.^(-(epsilon + u_s - x).*h_s)./(e.^(-(epsilon + u_s - x).*h_s) + 1).^2 - 2.*h_s.^2.*e.^(-2.*(epsilon + u_s - x).*h_s)./(e.^(-(epsilon + u_s - x).*h_s) + 1).^3;

Ps_xee = -h_s.^3.*e.^(-(epsilon + u_s - x).*h_s)./(e.^(-(epsilon + u_s - x).*h_s) + 1).^2 + 6.*h_s.^3.*e.^(-2.*(epsilon + u_s - x).*h_s)./(e.^(-(epsilon + u_s - x).*h_s) + 1).^3 - 6.*h_s.^3.*e.^(-3.*(epsilon + u_s - x).*h_s)./(e.^(-(epsilon + u_s - x).*h_s) + 1).^4;

Ps_e = -Ps_x;

Ps_ee = -Ps_xe;

% For debugging
if 0 == 1
    plot(x,Ps,'k');
    hold on
    plot(x,Ps_x,'r');
    plot(x,Ps_xe,'g');
    plot(x,Ps_xee,'b');
    plot(x,Ps_e,'m');
    hold off
end
