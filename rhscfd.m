function [bcv, rhsf, coefs ] = rhscfd( gridx, t )
%RHSCFD Summary of this function goes here
%   Detailed explanation goes here
%   Dirichlet conditions

n = length(gridx) - 1; % number of gridx
neq = n - 1;	% assume Dirichlet conditions
m = neq;	% simplicity has charm
hx = gridx(2:end) - gridx(1:end-1);

for i = 1:m
    px = gridx(i+1);
    [rhsf(i, 1), coefs(i, 1), coefs(i, 2), coefs(i, 3)] = PDEcoefs(px, t);
end

%modification of boundary conditions
ax = gridx(1); bx = gridx(end);
ua = truevd(ax, t); ub = truevd(bx, t); bcv = zeros(neq, 1);
bcv(1) = ((2*coefs(1, 3) - coefs(1, 2)*hx(2))/(hx(1)*(hx(1)+hx(2))))*ua;
bcv(end) = ((2*coefs(neq, 3) + coefs(neq, 2)*hx(end-1))/...
    (hx(end)*(hx(end-1)+hx(end))))*ub;

end