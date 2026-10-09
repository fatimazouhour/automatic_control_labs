%% AC simulation trial 

%% Q2

clear all
close all
clc

s=tf('s');
H=0.1*(s+2)/(2*s^3+6*s^2+5*s-1);
p=pole(H) % not bibo stable yss does not exist

%% Q3
clear all
close all
clc

s=tf('s');

k=10^(40/20);
H=k* (1-s/10)/(1-s/2)^2;

bode(H)  %to check

%unstable yss does not exist

%% Q4
clear all
close all
clc


A= [0 1;1/(cos(0))^2 0];
e=eig(A)


x = -pi/2:0.001:pi/2;

%this way 5brne 3na canale
line1= tan(x);
line2= 5*x; %pay attention to the sign elli lezm nhotta

figure(1)
plot(x,line1,'B','LineWidth',1.2)
hold on, grid on, zoom on   % aham shi nhot HOLD ON lbe2e ma droore hon
plot(x,line2,'r','LineWidth',1.2)
ylim([-20 20]) % 3shen el tan goes to inf fa byseer el line 8reeb so we have to limit it 3shen el result ttl3 ahsn

% they intersect at 3 pt 

%%  NOT RECOMMENDED
% other option fina n3ml function
f = @(x) tan(x)-5*x;
%tal3e root one by one
% el 2 heye wen you think mmkn tkoon el zero

%not recommended for accuracy 

%% Q5
clear all
close all
clc

format long  % 7ota bl treg 3shen e5d kaza decimal digit 

s=tf('s');
A=[0 1;-0.2 -0.5];
B=[0;1];
C=[1 0];
D=0;

lambda_des= [-15 -15];

k= acker(A,B,lambda_des)


N=1;
sysc= ss(A-B*k,B*N,C,D);
x0=[-0.1;0];

tsim = linspace(0,5,500000); %b third one kl ma ktrt points kl  ma more accurate
[y,t,x]= initial(sysc,x0,tsim);
norm= sqrt(x(:,1).^2+x(:,2).^2);

figure(1)
plot(t,norm,'B','LineWidth',1.2)
hold on, zoom on,grid on
yline(1e-5,':r')
yline(-1e-5,':r')
%fina n2raha bas saret less accurate 

%fina n3ml hek    
tol= 10^(-5);
range= find(abs(norm)<=tol,1);
treg=t(range)
%kl ma zdt pts bl linspace kl ma ken nb accurat aktr