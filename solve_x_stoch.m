function xs = solve_x_stoch(p,x0)

% -- xs = solve_x_stoch(p,x0)
% -- xs = solve_x_stoch(p)
%
% Solves the expected value of the fitness gradient equal to
% zero in order to to find the evolutionarily singular
% strategy for the breeding date in a stochastic
% environment.
%
% INPUTS
%
% p: A dictionary of parameter values. Use the file params.m
% for an example of how to specify one of these.
%
% x0: An initial guess for the breeding date. Optional.
%
%
% INPUTS
%
% xs: Evolutionarily singular breeding date strategy. 

if nargin < 2
    % Halfway between the survival and recruitment curves
    x0 = (p.u_s+p.u_r)/2;
end

f = @(x) calc_Eg(p,x);
[xs,fval] = fsolve(f,x0); % Finds xs such that calc_Eg returns 0
