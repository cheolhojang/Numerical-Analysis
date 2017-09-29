function [rts] = quartic3032219997(C)
%   Cheolho Jang
%   Math128a
%   Discussion 101
%   MULLERS METHOD WITH DEFLATION

% Input of given values
a = C(1);
b = C(2);
c = C(3);


% Inserted into the function x^4 + a x^3 + b x^2 + c x -1 = 0
quarticfun = [1 a b c -1];



% Parameters
TOL = 10^(-8);
iterations = 10^3;



% Special exception of a=b=c=0
if (abs(a) == 0)
    if (abs(b) == 0)
        if (abs(c) == 0)
            rts(1) = 1;
            rts(2) = (-1);
            return;
        end
    end
end

% Selection of random numbers
randomnumbers = randn(1,3);

% While statement of evaluation at these random numbers to see if we need to continue
while (polyval(quarticfun, randomnumbers(2)) - polyval(quarticfun, randomnumbers(1)) == 0|| polyval(quarticfun, randomnumbers(3)) - polyval(quarticfun, randomnumbers(2)) == 0|| polyval(quarticfun, randomnumbers(3)) - polyval(quarticfun, randomnumbers(1)) == 0)
    randomnumbers = 2*randomnumbers;
end


% Assignment of variables
p0 = randomnumbers(1);
p1 = randomnumbers(2);
p2 = randomnumbers(3);
h1 = p1 - p0; 
h2 = p2 - p1;
delta1 = (polyval(quarticfun, p1) - polyval(quarticfun, p0))/h1; 
delta2 = (polyval(quarticfun, p2) - polyval(quarticfun, p1))/h2; 
d = (delta2 - delta1)/(h2 + h1);


% MULLERS METHOD for the first root
i = 3;
while (i <= iterations)
    b = delta2 + h2 * d;
    D = (b^2 -4*polyval(quarticfun, p2)*d)^(1/2);
    if (abs(b - D) < abs(b + D))
        E = b + D;
    else
        E = b - D;
    end
    h = -2*polyval(quarticfun, p2)/E;
	p = p2 + h;
    if (abs(h)/norm(p) < TOL)
		rts(1) = p;
        break;
    end
    p0 = p1; 	
	p1 = p2; 
    p2 = p;
	h1 = p1 - p0;
	h2 = p2 - p1;
	delta1 = (polyval(quarticfun, p1) - polyval(quarticfun, p0))/h1;
	delta2 = (polyval(quarticfun, p2) - polyval(quarticfun, p1))/h2;
	d = (delta1 - delta2)/(h2 + h1);
	i = i + 1;
end

% We have the first root, now we need to look for the others
% by dividing the quartic function and 'deflating' it




% Divide the original eq by x - rts(1) and replace the quartic function with a cubic one
divider = [1 -rts(1)];
[q, r] = deconv(quarticfun,divider);
cuba = q(1);
cubb = q(2);
cubc = q(3);
cubd = q(4);

%Special Cases for cubic equation
if (cuba == 0)
    if (cubb == 0)
        if (cubc == 0)
            rts(2) = nan;
            return;
        end
        rts(2) = -d/c;
        return;
    end
    cuba = q(2);
    cubb = q(3);
    cubc = q(4);
    rts(2) = (-cubb + sqrt(cubb^2 - 4*cuba*cubc))/(2*cuba);
    rts(3) = (-cubb - sqrt(cubb^2 - 4*cuba*cubc))/(2*cuba);
    return;
end



i = 3;

%Restart the iteration process for the new cubic equation
while (i <= N0)
    b = delta2 + h2 * d;
    D = (b^2 -4*polyval(q, p2)*d)^(1/2);
    if (abs(b - D) < abs(b + D))
        E = b + D;
    else
        E = b - D;
    end
    h = -2*polyval(q, p2)/E;
    p = p2 + h;
    if (abs(h)/norm(p) < TOL)
        rts(2) = p;
        break;
    end
    p0 = p1;    
    p1 = p2; 
    p2 = p;
    h1 = p1 - p0;
    h2 = p2 - p1;
    delta1 = (polyval(q, p1) - polyval(q, p0))/h1;
    delta2 = (polyval(q, p2) - polyval(q, p1))/h2;
    d = (delta1 - delta2)/(h2 + h1);
    i = i + 1;
end

b2 = [1 -rts(2)];
[q2, r2] = deconv(q2,b2);
quada = q2(1);
quadb = q2(2);
quadc = q2(3);
rts(3) = (-quadb + sqrt(quadb^2 - 4*quada*quadc))/(2*quada);
rts(4) = (-quadb - sqrt(quadb^2 - 4*quada*quadc))/(2*quada);
return;
end

