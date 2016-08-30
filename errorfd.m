% function errg = errorfd(ngrid, gridx, n, uvct, nn, errg, truef);
function errg = errorfd(ngrid, gridx, n, uvct, T, nn, errg, truef);

if (nargin < 8) truef = 'truevd'; end;
[t1 t2 t3] = feval(truef, gridx(2:end-1), T);
errg(1, nn) = max(abs(t1-uvct'));
