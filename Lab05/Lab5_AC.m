%% P1
clear all
close all
clc

format compact
s = tf('s');
A = [-0.2 -1;1 0]; B=[0.5 0]'; C=[0 1]; D=0; x0 = [0 0]';
s_hat =0.06;
ts1 = 2;
Mr = ctrb(A,B);
rho_Mr = rank(Mr)
zeta = abs(log(s_hat))/sqrt(pi^2+(log(s_hat))^2)
wn = log((1/100)^(-1))/(zeta*ts1);

lambda1 = -wn*zeta +j*wn*sqrt(1-zeta^2);
lambda2 = -wn*zeta -j*wn*sqrt(1-zeta^2);
lambda3 = -5*zeta*wn;
lambda_des = [lambda1 lambda2 lambda3]

Aaug = [0 -C;zeros(2,1) A]; 
Baug = [0 ;B];
Br = [1 ;zeros(2,1)];
Caug = [0 C];
% now i can get kaug 
kaug = acker(Aaug,Baug,lambda_des)
kq = kaug(1)
kx=kaug(2:end)

A_tildeaug = [0 -C;-B*kq A-B*kx]; B_tildeaug= Br; C_tildeaug =[0 C];
sys = ss(A_tildeaug,B_tildeaug,C_tildeaug,0);
tsim = linspace(0,20,10000);
[y,t,x_aug]= step(sys,tsim);

figure(1)
subplot(211), plot(t,x_aug(:,2),'b')
hold on , zoom on, grid on, xlabel('t(s)'), ylabel('x1(t)')
subplot(212)
plot(t,x_aug(:,3),'b')
grid on , zoom on, hold on,  xlabel('t(s)'), ylabel('x2(t)')

figure(2)
plot(t,y,'b')
grid on , zoom on, hold on,  xlabel('t(s)'), ylabel('y(t)')


% settling time ok , overshoot ok 
%% P2
clear all
close all
clc

format compact 
s=tf('s');
x0 = [-0.5;0];A=[-0.2 -1;1 0];B=[0.5 0]';C=[0 1];D=0;
zeta = 0.67;
wn=3.45;
Mr = ctrb(A,B);
rho_Mr= rank(Mr);
lambda1 = -zeta*wn+j*wn*sqrt(1-zeta^2);
lambda2 = -zeta*wn-j*wn*sqrt(1-zeta^2);
%n=2 --> only 2 lambdas
lambdas = [lambda1 lambda2]
k = acker(A,B,lambdas)
sysN = ss(A-B*k,B,C,D);
N=1/dcgain(sysN)
%now do observer
Mo = obsv(A,C);
rho_Mo = rank(Mo)
lambda_obsv_des = [-zeta*wn -zeta*wn]*5
L = acker(A',C',lambda_obsv_des)';
Ai=[A (-B*k);(L*C) (A-B*k-L*C)];
Bi = [B; B]*N;Ci = [C zeros(1,2)];
sys =ss(Ai,Bi,Ci,0);

x0hat=[0;0];
% simulation to get the rest 
tsim = linspace(0,20,10000);
usim = ones(size(tsim));

[y,t,xi]=lsim(sys,usim,tsim,[x0;x0hat]);
%now plot
r = stepinfo(sys,"SettlingTimeThreshold"',0.02)
figure(1)
plot(t,y,'b'), grid on, zoom on, hold on,xlabel('t(s)'),ylabel('y(t)')

%%    check!!!!!!!!!!!
clear all 
close all
clc
format compact 
A = [-0.2 -1;1 0]; B =[0 1]';x0 = [1 0.5]';D = 0;
Mr = ctrb(A,B);
rho_Mr = rank(Mr)
sys = ss(A,B,B',D);
[y,t,x] = initial(sys,x0);
normx = sqrt(x(:,1).^2+x(:,2).^2);
figure(1),
plot(t,normx), grid on, zoom on, xlabel('t(s)'),ylabel('normx2')

tolreg = 1e-5;
tregdes=0.8;
lambda=log(tolreg/norm(x0))/tregdes
lambda1 = lambda
lambda2 = 100*lambda1
lambdades =[lambda1 lambda2];
k =place(A,B,lambdades)

%%
clear all
close all
clc

format compact
s= tf('s');
A= [-1 0;0 10]; B=[1 1]';C=[1 0];D=0;
H = tf(ss(A,B,C,D));
lambdaA= eig(A)

Mr = ctrb(A,B); rhoMr = rank(Mr)
Mo = obsv(A,C); rhoMo = rank(Mo)

%% 
clear all
close all
clc

format compact
s= tf('s');
A= [-1 0;0 10]; B=[1;1];C=[1 0];D=0;
k = [0 11]; N= -0.1;
sys = ss(A-B*k,B*N,C,D);
tsim = linspace(0,20,10000);
[y,t,x] = step(sys,tsim);

figure(1)
plot(t,y,'g','linew',1.5)
grid on, hold on, zoom on
xlabel('t(s)'),ylabel('y(t)')

shat = (max(y)-y(end))/y(end);
r = stepinfo(sys)

%now PART 2
x0 =[0 0]';
Ai = A-B*k;
Bi = B*N;
r = 1/(s+1);
Y= minreal(zpk(C*inv(s*eye(2)-Ai)*(x0+Bi*r)));
[numy,deny]=tfdata(Y,'v');
[r,p]=residue(numy,deny)
