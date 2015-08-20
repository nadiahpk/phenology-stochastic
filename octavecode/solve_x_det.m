function xs = solve_x_det(p,x0);
% Solve for the deterministic singular strategy xs

if nargin < 2
    % *** choose a better one later
    x0 = -2;
end

f = @(x) calc_w(p,x);
% An evolutionarily singular strategy is when that
% gradient is zero
[xs,fval] = fsolve(f,x0);
