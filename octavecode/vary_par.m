function M = vary_par(p,parName,parRange,nopts);

% Example:
%
% clear all
% params;
% parName = 'V_ol';
% parRange=[0,5];
% nopts=20;
% M=vary_par(p,parName,parRange,nopts);


str = ['output',parName,'.dat'];
fid = fopen(str,'w');
fprintf(fid,'# Generated with vary_par.m\n');
fprintf(fid,'#\n');
fprintf(fid,'# Default parameter values were:\n');
for [val,key] = p;
    str = ['# p.',key,' = ',num2str(val),';\n'];
    fprintf(fid,str);
end
fprintf(fid,'#\n');

xs = solve_x_det(p)
n=calc_n(p,xs)

fprintf(fid,['# xs_det = ',num2str(xs),'\n']);
fprintf(fid,['# n = ',num2str(n),'\n']);
fprintf(fid,'#\n');


fprintf(fid,'# Column entries are:\n');
fprintf(fid,'# 1:parVal \t 2: xs_stoch \t 3: E_v \n');
fprintf(fid,'\n');

parV = linspace(parRange(1),parRange(2),nopts);
M = []; % parVal, xs_det, xs_stoch, n_det, Ev
for i = 1:nopts;
    
    % Update parameters dictionary
    eval(['parVal = ',num2str(parV(i)),';']);
    eval(['p.',parName,' = parVal;']);

    % Calculate each bit
    xs_stoch = solve_x_stoch(p,xs);
    [V_e,V_l,C_el,C_ve,C_vl,V_v,E_v] = calc_moments(p,xs_stoch);

    % Store
    M = [M; parVal xs_stoch E_v];
    fprintf(fid,'%0.4e \t %0.4e \t %0.4e \n',M(i,:));
end

fclose(fid);
