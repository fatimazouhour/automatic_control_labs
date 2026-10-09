%% P1
clear all
close all
clc

%define kl shi 3nde
A=[-0.2 -1;1 0];
B= [0.5 0]';
C=[0 1];
D=0;
x0=[0 0]';
s = tf('s');
r=1/s;

s_hat = 0.06;
ts1 = 2;

% checking reachability 
Mr = ctrb(A,B);
rhoM= rank(Mr)
% now get zeta and wn requirements

zeta = abs(log(s_hat))/sqrt(pi^2+log(s_hat)^2)
wn = log(0.02)^(-1)/(ts1*zeta)

% eig
lmbda1 = -zeta*wn+j*wn*sqrt(1-zeta^2)
lmbda2 = -zeta*wn-j*wn*sqrt(1-zeta^2)
lmbdades = [lmbda1 lmbda2];

k = place(A,B,lmbdades)
Ac = A-B*k;
Bc=B;
Cc=C;
Dc=0;
sys_N= ss(Ac,Bc,Cc,Dc);
N= 1/dcgain(sys_N)

% final
sys = ss(Ac,B*N,C,D);



%% P2
clear all
close all
clc

A=[-1.8 -2.8 -2; 1 0 0;0 1 0];
B=[1 0 0]';
C=[0 0.25 1];
D=0;
x0=[0 0 0]';
s=tf('s');

Mo=ctrb(A,B);
rho = rank(Mo);
shat=0;
ts1=1.5;
zeta=1;
wn = log(100/2)/(ts1*zeta)
lambda1 = -zeta*wn+j*wn*sqrt(1-zeta^2);
lambda2 = -zeta*wn-j*wn*sqrt(1-zeta^2);
lambda3 = -10*wn*zeta;

lambdades=[lambda1 lambda2 lambda3]
k= acker(A,B,lambdades);
sys_N = ss(A-B*k,B,C,D);
N=1/dcgain(sys_N)
sys_cont= ss(A-B*k,B*N,C,D);
%now i want to plot 
tsim = linspace(0,20,10000);
[y,t,x]= step(sys_cont,tsim);
figure(1)
subplot(311)
plot(t,x(:,1),'b','Linew',1.5)
grid on, zoom on, hold on, xlabel('t (s)'), ylabel('x_1(t)')
subplot(312)
plot(t,x(:,2),'b','Linew',1.5)
grid on, zoom on, hold on, xlabel('t (s)'), ylabel('x_2(t)')
subplot(313)
plot(t,x(:,3),'b','Linew',1.5)
grid on, zoom on, hold on, xlabel('t (s)'), ylabel('x_3(t)')

%% P3
clear all
close all
clc

A=[-1.8 -2.8 -2; 1 0 0;0 1 0];
B=[1 0 0]';
C=[0 4 1];
D=0;
x0=[0 0 0]';
s=tf('s');

Mo=ctrb(A,B);
rho = rank(Mo);
shat=0;
ts1=1.5;
zeta=1;
wn = 6.64/ts1

lambda1 = -zeta*wn+j*wn*sqrt(1-zeta^2)
lambda2 = -zeta*wn-j*wn*sqrt(1-zeta^2)
lambda3 = -10*wn*zeta
lambdades = [lambda1 lambda2 lambda3]
%% RESOLVE 
%% P2 
clear all
close all
clc
format compact 
s = tf('s');
A = [-1.8 -2.8 -2;1 0 0;0 1 0]; B = [1 0 0]';C=[0 0.25 1];
s_hat = 0;
ts1 = 1.5;
x0 = [0 0 0]';
Mr = ctrb(A,B);
zeta = 1;
rho_Mr = rank(Mr)
wn = log((1/100)^(-1))/(zeta*ts1);
lambda1 = -zeta*wn +j*wn*sqrt(1-zeta^2);
lambda2 = -zeta*wn -j*wn*sqrt(1-zeta^2);
lambda3 = -15*zeta*wn;
lambda_des = [lambda1 lambda2 lambda3]
k = acker(A,B,lambda_des)
sys_N = ss(A-B*k,B,C,0);
N = 1/dcgain(sys_N)

tsim = linspace(0,20,10000);
sys =ss(A-B*k,B*N,C,0);
[y,t,x] = step(sys,tsim);

w = tf(sys);
u=1/s;
y = zpk(minreal(w*s,1e-3));
[numy,deny] =tfdata(y,'v');
[z,p,k] = residue(numy,deny)


r = stepinfo(y,t,'SettlingTimeThreshold',0.01)
% figure(1), subplot(311), plot(t,y,'b','linew',1.5)
% grid on , zoom on, hold on
% xlabel('t(s)'),ylabel("y(t)")
% 
% subplot(312), plot(t,x(:,1),'b','linew',1.5)
% grid on , zoom on, hold on
% xlabel('t(s)'),ylabel("x1(t)")
% 
% subplot(313), plot(t,x(:,2),'b','linew',1.5)
% grid on , zoom on, hold on
% xlabel('t(s)'),ylabel("x2(t)")


