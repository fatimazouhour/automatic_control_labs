%% slide 8 example

clear all
close all
clc

format default

A = [-1 1;0 1];
B=[2 0]';
C=[0 -1];

Mr = ctrb(A,B);
rho_Mr = rank(Mr)

% answer is rho_Mr=1< n=2  (la n3rf n nshuf matrices e.g. A heye nxn hon
% msln heye 2x2)

%hence here the system is unreachable and unstable 

%% slide 11
clear all
close all
clc
format default

A = [1 3;4 2];
B =[-1 2]';
C = [1 1];

Mr = ctrb(A,B)
rho_Mr = rank(Mr)
% rank=n=2 --> ok

% we now want k 
% we are asked to find k such that Lambda1=-2 and lambda2=-3

lambda_des =[-2 -3];
k = acker(A, B, lambda_des)
% or use the place()

sys_cont_c = ss(A-B*k,B,C,0);

N=1/dcgain(sys_cont_c)

%% FROM NOW ITS FROM THE PDF ABOUT L8 EXAMPLES

%% example 1
clear all
close all
clc

format default

A =[-0.1 -1;1 0];
B=[1 0]';
C = [0 1];
D =0;
% assume x0=[0 0];
x0 = [0; 0];
Mr = ctrb(A,B);
rho_Mr = rank(Mr)
% rank = n = 2 --> ok
%  1. define the given and get zeta (check L8 slide 15)

s_hat = 0.1;
ts2 = 6.5;
zeta = abs(log(s_hat))/(sqrt(pi^2+log(s_hat))^2)
% 2. now get wn 

wn = log((2/100)^(-1))/(zeta*ts2)

% 3. assign the eigenvalues wehde pos wehde neg

lmbda1 = -zeta*wn+j*wn*sqrt(1-zeta^2)
lmbda2 = -zeta*wn-j*wn*sqrt(1-zeta^2)
lambda_des = [lmbda1  lmbda2]

% 4. now we can calulate k
k = place(A,B,lambda_des)

% 5. check matrices (new)
Ac = A-B*k;
Bc= B;
Cc=C;
Dc=0;
sysN = ss(Ac,Bc,Cc,Dc);
N=1/dcgain(sysN)

% now we have the controlled system

sys_cntrl =ss(A-B*k,B*N,C,D);

% simulation

tsim =  linspace(0,20,10000);
[y,t,x0] = step(sys_cntrl,tsim);
%we plot

% figure(1),
% subplot(211)
% plot(t,x(:1))
% grid on, zoom on, hold on, xlabel('t(s)')

%% solving problem 1 alone 
clear all
close all
clc 

A = [-0.1 1;1 0]; B=[1 0]';C = [0 1]; x0=[0 0]'; D = 0;
s = tf('s');
% design control law 

% awal shi ntale3 k
% we check reachability / contability of the system 
Mr = ctrb(A,B);
rho_Mr = rank(A)  % = 2 = n y3ne tmm
shat = 0.1;
ts2 = 6.5;
zeta = abs(log(shat))/sqrt(pi^2+(log(shat)^2))
wn = log((2/100)^(-1))/(ts2*zeta)

% since n= 2 y3ne ana bde 2 labdas bs 
lambda1 = -zeta*wn+j*wn*sqrt(1-zeta^2)
lambda2 = -zeta*wn-j*wn*sqrt(1-zeta^2)

lambda_des = [lambda1 lambda2];
k = place(A,B,lambda_des)
sys = ss(A-B*k,B,C,D);
N= 1/dcgain(sys)
sysN = ss(A-B*k,B*N,C,D);
% yayyy easyyy
% we should do the simulation of the control system
sysx= ss(A,B,eye(2),D);
tsim = 20;
sim('Lecture8_example1_sim')
