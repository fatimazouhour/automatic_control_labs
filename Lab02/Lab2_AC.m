clear all
close all
clc

format default
A =[-2 0 0;0 0 1;0 0 0];
s = tf('s');
Q = zpk(minreal(inv(s*eye(3)-A),1e-3))
eig(A)

%% ex2
clear all
close all
clc
format default

A=[-1 0 0;0 0 1;0 -1 0];
s=tf('s');
q=zpk(minreal(inv(s*eye(3)-A),1e-3))

%% ex3 
clear all
close all
clc


A=[1 0 0;0 -2 3;0 -3 -2];
B=[0 1 1]';
C= [0 1 2];
s=tf('s');

q=zpk(minreal(inv(s*eye(3)-A),1e-3))
eig(A)
H = zpk(minreal(C*inv(s*eye(3)-A)*B,1e-3))
[numh,denh]=tfdata(H,'v');
[r,p]=residue(numh,denh)

%% p4
clear all
close all
clc
% abl ma tshufe sol shufe el x0 ade heye 
format default
s=tf('s');
A= [0 1 0 0;-1 -1 0 0;0 0 -10 0;0 0 0 -5];
B=[1 0 0 1]';
C=[1 0 1 1];
D=0;
U = 1/s +1/s^2;
x0=[0 0 1 1]';

%1.response y(t):
% sys = ss(A,B,C,D);
% H=tf(sys)
% y = zpk(minreal(H*U,1e-3))
eig(A)
Y=zpk(minreal(C*inv(s*eye(4)-A)*B*U,1e-3))
[numy,deny]= tfdata(Y,'v');
[r,p]=residue(numy,deny)
M=abs(r) ,2*M
theta = angle(r)
%2. response with another x0
% ana bas ha 8ayer mn fo2 el x0 mish aktar 
x= zpk(minreal(inv(s*eye(4)-A)*x0,1e-3))


%3.modal analysis + time constant


