clear all
close all
clc

format compact 
s=tf('s');

Gp= 10/(s*(s+125));
p_Gp=pole(Gp)


Gs=5;
Ga=3;
Gda=1;
Gdp=1;
kd=1;

Gf=kd/Gs;

%s1

%s2

%s3

as=5*1e-3;
ws=2500;
rho_s3= 1*1e-4;
Mth= 20*log10(rho_s3*Gs/as)

wh= ws*10^(Mth/40);
wc_s3= wh/2;


%req transient 
shat=0.08;
ts5= 0.1;
tr= 0.05;

%s6
zeta = abs(log(shat))/sqrt(pi^2+(log(shat))^2);
Tpmax= 1/(2*zeta*sqrt(1-zeta^2));
Spmax= 2*zeta*sqrt(2+4*zeta^2+2*sqrt(1+8*zeta^2))/(sqrt(1+8*zeta^2)+4*zeta^2-1);

%s5
rho_tr= 1/sqrt(1-zeta^2)*(pi-acos(zeta))*sqrt(sqrt(1+4*zeta^4)-2*zeta^2);
wc_tr=rho_tr/tr;

%s5
rho_ts5= -log(5)/zeta * sqrt(sqrt(1+4*zeta^4)-2*zeta^2);
wc_ts5= rho_ts5/ts5;


wc_min= max(wc_ts5,wc_tr)
wc_max= wc_s3
%
nu =1;
Lin = 1/(s^nu)*Gp*Ga*Gf*Gs;

figure(1),
bode(Lin)

figure(2),
[num,den]=  tfdata(Lin,'v');
nyquist1(num,den)

% figure(3)
% myngridst(Tpmax,Spmax)
% hold on,
% nichols(Lin)


% 1st trial
wc_des = 50;
wc_norm1= 100;

z=wc_des/wc_norm1

Rz= (1+s/z);
kc_new= 10^(42/20)
L_new= 1/s^nu*Gp*Gs*Ga*Gf*Rz*kc_new;

figure(4)
myngridst(Tpmax,Spmax)
hold on,
range=logspace(-3,4,500);
nichols(L_new,range)

%evaluation 
T=minreal(zpk(L_new/(1+L_new)));
S=minreal(zpk(1/(1+L_new)));

figure(5)
step(T*kd,1),
figure(6)
bodemag(T),



Gc= zpk(kc_new/s^nu*Rz)

