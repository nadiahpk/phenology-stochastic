% Survival and recruitment sigmoids
p.u_s = -5; % Midpoint of sigmoid
p.h_s = -.5; % Slope of sigmoid
p.u_r = 5; 
p.h_r = .5;

% Territory competition
p.m = 0.5; % Slope term
%p.m = 4; % Slope term

% Other fixed parameters
p.K = 100;
p.s = 0.65; % lit
p.a = 3; % lit

% Stochastic parameters
%  - early season
p.G_e = 0.1; % Autocorrelation early-season
p.V_oe = 5; % Variance of underlying white-noise
%  - late season
p.G_l = 0.2; 
p.V_ol = 4;
%  - correlation between early and late
p.C_oeol = 1;
% xs = -3.7775; n = 215.77;

