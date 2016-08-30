%Test function set-up
format compact
global Uno Uname;
Uno = 209; Uname = '';

%Initial Setup
ax = 0; bx = 1;%boundary
T = 1; %time
theta = 1/2; % 0 for forward, 1 for backward, 1/2 for Crank-Nicolson
ntimes = 6; errg = zeros(1, ntimes);
Nt = 100; ht = T/Nt; % uniform timestep

for nn = 1:ntimes
    n = 2^(nn+2); ngrid = n+1; nint(nn) = n;
    neq = n-1;
    
    %set up of non-uniform grid
    hx_u = (bx-ax)/n;
    gridx_u = ax + hx_u*[0:n];
    gridx = gridx_u.^2;
    hx = gridx(2:end) - gridx(1:end-1);
    
    %initial value of u
    u = truevd(gridx(2:end-1), 0)';
    t = 0;
    
    %solving IVP
    [bcv, rhsf, coefs ] = rhscfd( gridx, 0 );
    A = cfdmat(gridx, coefs);
    
    while t < T
        
        %timestep
        ht = ht;  %temporary time step
        t = t + ht;
        
        %backward matrix fot this timestep
        [bcv2, rhsf2, coefs2] = rhscfd( gridx, t );
        A2 = cfdmat(gridx, coefs2);
        
        %setup of system and solve
        rhs = ( speye(size(A))/ht+(1-theta)*A ) * u ...
         + rhsf*(1-theta) + rhsf2*theta + theta*bcv2 + (1-theta)*bcv;
        AA = speye(size(A))/ht - theta*A2;
        u = AA\rhs;
        
        bcv = bcv2; rhsf = rhsf2; coefs = coefs2; A = A2;
    end
    
    %Inf-norm abs error (approximation and true value)
    errg = errorfd(ngrid, gridx, n, u, T, nn, errg);
    
end

%OUTPUT
[udummy] = truevd(ax, T);
disp(['U = ' Uname ' = {' num2str(Uno) '}'])
disp(['domain [' num2str(ax)  ', ' num2str(bx) '] at time ' num2str(T)])

nint
format short e
disp('error on grid points')
errg
format short
disp('order of convergence')
errg(:, :) = max(errg(:, :),  0.222044604925e-15);
LogNintRatio = log(nint(1, 2:ntimes)./nint(1, 1:ntimes-1));
LogNintRatioMat = repmat(LogNintRatio, size(errg, 1), 1);
if ntimes > 1
    convg = log(errg(:, 1:ntimes-1)./errg(:, 2:ntimes))./LogNintRatioMat
end