clear all
clc
%define first the vectors

A = [0 1;-2 -3];
B = [1;0];

%introduce the variable s according to this statement
s= tf('s');
U = 2/s;
x0 = [2;2];
%the solution inlaplace domain can be found on slide 35 (L02)
%we can collect on the left (SI-A)^-1 [x(0) + BU(s)] can be computed using
%one line in matlab

X = zpk(minreal(inv(s*eye(2)-A)*(B*U+x0), 1e-3))
%return  htyneha rj3na 3mlnela comment bas tanshuf lshe8el till now

%1 is the first component and 2 is the second

[num_X1, den_X1] = tfdata(X(1),'v') % hon we want to compute the numerator (num) and denominator (den) so we can compute the coeffiients
[r1,p1] = residue(num_X1, den_X1)

% in P we have the roots in R the coefficients
%the first coefficient in r1 corresponds to the first in the other andso on
% hl2 sar fina nektob X= -3/(s+2) +2/(s+1) +3/s

% when we use the tfdatafunction bytl3 3na 2 vectors nfs el size aw eli howe

[num_X2, den_X2] = tfdata(X(2),'v')
[r2,p2] = residue(num_X2, den_X2)

%% example 2
clear all
clc

A = [-3 2;-2 -3]
B = [1;0]
C = [0 1]
x0 = [1;1]
s = tf('s');
U = 1/s

Y = zpk(minreal(C*inv(s*eye(2)-A)*(B*U+x0),1e-3))


% it shows a second degree polynomial

[num_Y, den_Y] = tfdata(Y,'v')
[r,p] = residue(num_Y,den_Y)
% in P we have the roots in R the coefficients
%now we just put together 
% notice eno sar 3na + imaginary part

M = abs(r(1)), 2*M
phi = angle(r(1))

%hek tl3na kl shi to write the solution
