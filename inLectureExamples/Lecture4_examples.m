%%  slide 34 example
clear all
close all
clc
format default
A= [0 1 0 0;-4 -4 0 0;0 0 -0.2 0; 0 0 0 -0.3]
eig(A)
% since the eigs are all negative --> asymptotically stable

%% slide 35 example 2
clear all
close all
clc

format default
s = tf("s")
A = [0 0 0;0 0 1;0 0 -5];
eig(A)
% multiplicity for  eig=0 is 2 -->minimal polynomial HAS to be computed
% we compute it as (sI-A)^-1 =Q

Q = zpk(minreal(inv(s*eye(3)-A),1e-3))

%% slide 39 example
clear all
close all
clc

format default 

A=[0 1 0;-1 0 0;0 0 -1];
eig(A)
% internally stable no need for further action to know


%%
clc
s = tf('s');
A = [0 1;-2 -3]; B= [1 0]'; x0=[2 2]';
x = zpk(minreal(inv(s*eye(2)-A)*(x0),1e-3))
[num1,den1] = tfdata(x(1),'v');
[r,p]= residue(num1,den1)

[num2,den2] = tfdata(x(2),'v');
[r2,p2]= residue(num2,den2)
lambda = eig(A)

M = abs(r(1)) , 2*M
phi = angle(r(1))

M2 = abs(r2(1)) , 2*M2
phi2 = angle(r2(1))

