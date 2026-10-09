%% example slide 12
clear all
close all
clc

format default
A= [0 1;-1 2]
B=[0 1]'
C=[0 1];
D=0;

Mo = ctrb(A,B);
rho = rank(Mo)
lambda_des = [-10 -20];

L = place(A',C',lambda_des) % be awareeno A w C lezm ykono transposed

%% pdf L9 examples

%% example 1
clear all
close all
clc

format default
A=[-0.4 -1;1 0];
B=[1 0]';
C=[0 1];
D =0;

x0=[-0.1 0.5]';
xhat=[0 0]';
sys = ss(A,B,C,D);

Mo= obsv(A,C);
rho = rank(Mo) % n = rank=2 --> ok
lambda_obs = [-0.5 -0.5];
L=acker(A',C',lambda_obs)';

sys_obv = ss(A-L*C,[B L],eye(2),0);
tsim = linspace(0,50,10000);
usim = 0.1*ones(size(tsim));
[y,t,x] = lsim(sys,usim,tsim,x0);

[yhat,that]=lsim(sys_obv,[usim',y],tsim,xhat);

% now to plot the results

figure(1),
subplot(211)
plot(t,x(:,1),'b','linew',1.5)
grid on ,zoom on, hold on,
plot(t,xhat(:,1),':r','linew',1.5)

%% example 2

clear all
close all
clc

format default 
A=[-0.4 -1;1 0];
B=[1 0]';
C=[0 1];
D=0;
x0=[-0.5 0]';
x0hat =[0 0]';

Mr = ctrb(A,B);
rho_Mr =rank(Mr) % n=2=rank --> ok
Mo= obsv(A,C);
rho = rank(Mo) % n =rank =2 --> ok

zeta = 0.63;
wn = 1.84;

lmda1 = -zeta *wn +j*wn*sqrt(1-zeta^2)
lmda2 = -zeta *wn -j*wn*sqrt(1-zeta^2)
lambda_des = [lmda1,lmda2];

k = place(A,B,lambda_des)

Ac = A-B*k;
Bc= B;
Cc=C;
Dc=0;
sys_c=ss(Ac,Bc,Cc,Dc);
N=1/dcgain(sys_c)

%compute L
lambda_obs =[-zeta*wn -zeta*wn]*5 %why *5??
L=acker(A',C',lambda_obs)'

