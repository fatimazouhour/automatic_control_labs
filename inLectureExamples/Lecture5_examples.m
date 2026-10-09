clear all
close all
clc
format default
A = [0 1 0 0; -1 0 0 0; 3 0 0 0; 0 0 2 0];
B=[1;0;0;0];
C=[1 0 0 0];

s= tf('s');
U= 1/(s^2+1);
y = zpk(minreal(C*inv(s*eye(4)-A)*B*U,1e-3))

[num_y,den_y]= tfdata(y,"v") % y3ne 3teene hene as a vector
[r,p] = residue(num_y,den_y)

% hl2 mnhot hol b aleb polynomial
theta  = angle(r)
deg = theta*180/pi
% now we write those as an equation manually and solve 

%% now this section is for me trying it alone nothing more
clear all
close all
clc

format compact
s = tf('s');
A = [0 1 0 0;-1 0 0 0;3 0 0 0;0 0 2 0]; B=[1 0 0 0 ]'; C=[1 0 0 0]; x0=[0 0 0 0]';
U = sqrt(2)/(s^2+2);
% y =??  note eno fine tale3a using H l2n 3nde zero initial cond
y = zpk(minreal(C*inv(s*eye(4)-A)*(x0+B*U),1e-3))
[numY,denY] = tfdata(y,'v');
[r,p] = residue(numY,denY)

% now ha 8ayyer to sin(tsqrt2))