% function [true1, true2, true3] = truevd(x,t)
%
% returns the values of the solution function, 1st and 2nd
% derivatives on x at time 

function [true1, true2, true3, truet] = truevd(x, t)

global Uno Uname;

o = ones(size(x)); z = zeros(size(x));

switch Uno
    
    case {309}
        Uname ='exp(-t)*x.^4';
        true1 = exp(-t)*x.^4;   true2 = 4*exp(-t)*(x.^3);  
        true3 = 12*exp(-t)*(x.^2);
        truet = -exp(-t)*(x.^4);
    
    case {209}
        Uname ='exp(-t)*sin(x)';
        true1 = exp(-t)*sin(x);   true2 = exp(-t)*cos(x);  
        true3 = -exp(-t)*sin(x);
        truet = -exp(-t)*sin(x);

    case {109}
        Uname ='exp(x)*t';
        true1 = exp(x).*t;     true2 = exp(x).*t;       true3 = exp(x).*t;
        truet = exp(x);
        
    case {101}
        Uname ='x*t';
        true1 = x.*t;       true2 = t.*o;       true3 = z;     truet = x;
        
    case {111}
        Uname ='x+t';
        true1 = x+t;   true2 = o;   true3 = z;
        truet = o;
        
    case {21}
        Uname ='x.^(2)*t';
        true1 = (x.^(2)).*t;   true2 = (2.*x).*t;   true3 = 2*t*o;
        truet = (x.^(2));
        
    case {22}
        Uname ='x.^2*t^2';
        true1 = x.^2*t^2;       true2 = 2.*x*t^2;       true3 = 2.*o*t^2;
        truet = 2*x.^2*t;
        
        
    case {11}
        Uname ='x.*x';
        true1 = x.*x;       true2 = 2*x;          true3 = 2*o;
        truet = z;
     
    case {10}
        Uname ='t';
        true1 = t.*o;          true2 = z;            true3 = z;
        truet = o;
        
    case {1}
        Uname ='x';
        true1 = x;          true2 = o;            true3 = z;
        truet = z;
        
    case {0}
        Uname ='1';
        true1 = o;          true2 = z;            true3 = z; 
        truet = z;
        
    otherwise
        error(['truevd: no such function ' num2str(Uno)])
end
