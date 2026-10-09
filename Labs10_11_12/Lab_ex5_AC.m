%% !!!!EX5!!!!
clear all
close all
clc

s=tf('s');

Gp=25/(s^3+3.3*s^2+2*s);
p_gp=pole(Gp)


Ga=0.38;
Gs=2;
Gr=1;
Gd=1;
Gf=1/8;
Da=5.5*1e-3;
ap=2*1e-2;
as=0.1;
kd=4;
kp=dcgain(Gp*s);

wp=0.02;
ws=40;

%s2
rho_r=1.5*1e-1;
kc_s2= kd^2*1/(kp*Ga*rho_r)

%s3
rho_da=5.8;
kc_s3= Da/(5.8*Ga*Gf*Gs)

%s4
rho_dp=3.6*1e-4;
S= rho_dp/ap;
Ms= 20*log10(S)
wL= wp*10^(-Ms/40)
wc_S = 2*wL 

%s5
rho_ds=1.25*1e-4;
T= rho_ds*Gs/as;
Mh= 20*log10(T)
wh=ws*10^(Mh/40)
wc_T=wh/2


%s7
ts5=5;

%s8
shat=0.09; % changed by prof to make the question harder
zeta= abs(log(shat))/sqrt(pi^2+(log(shat))^2)

Tpmax=1/(2*zeta*sqrt(1-zeta^2))
Spmax= 2*zeta*sqrt(2+4*zeta^2+2*sqrt(1+8*zeta^2))/(sqrt(1+8*zeta^2)+4*zeta^2-1)

%s6
tr=2.5;
rho_tr= 1/sqrt(1-zeta^2) * (pi-acos(zeta)) *sqrt(sqrt(1+4*zeta^4)-2*zeta^2);
wc_tr= rho_tr/tr

%s7
rho_ts = -log(0.05)/zeta * sqrt(sqrt(1+4*zeta^4)-2*zeta^2);
wc_ts= rho_ts/ts5


%
wc_min=max([wc_ts  wc_tr  wc_S])
wc_max= wc_T

%
kc_in =0.015;

Lin= kc_in/s * Gp*Ga*Gf*Gs;
figure(1)
bode(Lin)
[num,den]=tfdata(Lin,'v'),
nyquist1(num,den)


figure(2)
myngridst(Tpmax,Spmax)
nichols(Lin)

% figure(3)
% nyquist(Lin)

% let me choose wc =0.9 
%v=1 so i can add a zero network 

wc_des= 0.9;


z= wc_des/20;
Rz= (1+s/z);
zd=wc_des/3.4;
md=10;
Rd= (1+s/zd)/(1+s/(md*zd))
L_new= Lin*Rz*Rd;

figure(2)
hold on
nichols(L_new)
