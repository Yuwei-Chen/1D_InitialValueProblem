function [rhsf, coefu, coefux, coefuxx ] = PDEcoefs( x, t )
%coefficients returns the values of the PDE coefficient functions
%             and RHS function f(x)
%   DEs u_t = p(x)*u_xx + q(x)*u_x + r(x)*u + f(x, t), Dirichlet BCs

coefu     = t+1+1./(1+x); % 1+1./(1+x);
coefux    = t^2+cos(x); % cos(x);
coefuxx   = exp(x)+exp(-t); % -exp(x);

% rhsf = rhsfunc(x); %self-defined RHS function f(x,t) of parabolic PDE

% rhsf used to test when knowing the true value of functions
[true1, true2, true3, truet] = truevd(x, t);
Lu = coefu*true1 + coefux*true2 + coefuxx*true3;
rhsf = truet - Lu;
end