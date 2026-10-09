%lab problem 1

clear all
close all
clc

s=tf('s');

Gp = 25/(s^3+3.3*s^2+2*s);
[r,p,kp]=tfdata(Gp,'v')

%given 
Gs=1;
Ga=0.095;
Gr=1;
Dao=5.5*1e-3;
ap=2*1e-2;
as=1e-1;

wp=0.02;
ws=40;

%1
kd=1;
Gf=kd/Gs;

%2
rho_s2= 1.5*1e-1;
Ro=1;
kc_2= kd^2*Ro/(12.5*Ga*rho_s2)

%3
rho_s3= 1.5*1e-2;
Kc_3= Dao/(Ga*Gs*Gf*rho_s3)

%4
rho_s4= 5*1e-4;
rho_s5=rho_s4;

Ms_tf = 20*log10(rho_s5/ap)
wl= wp*10^(- Ms_tf/40);

wc_s = 2*wl

%5
Mt_tf= 20*log10(rho_s5*Gs/as)
wh= ws*10^(Ms_tf/40)

wc_t=wh/2

%transient req
shat = 10;
tr=3;
ts5=12;

%7

zeta=abs(log(shat))/(sqrt(pi^2+(log(shat))^2));
Tpmax= 1/(2*zeta*sqrt(1-zeta^2));
Spmax= 2*zeta*sqrt(2+4*zeta^2+2*sqrt(1+8*zeta^2))/(sqrt(1+8*zeta^2)+4*zeta^2-1);

%6
rho_tr= 1/sqrt(1-zeta^2)*(pi-acos(zeta))*sqrt(sqrt(1+4*zeta^4)-2*zeta^2);
wc_tr= rho_tr/tr

%7
rho_ts5= -log(5)/zeta*sqrt(sqrt(1+4*zeta^4)-2*zeta^2);
wc_ts5= rho_ts5/ts5

%
wc_min= max([wc_tr,wc_ts5,wc_s])
wc_max= wc_t

nu=0;
kcin=6;
Lin= kcin/s^nu *Gp*Gs*Gf*Ga;
[rl,pl]=tfdata(Lin,'v')

w_des= 1;
%new L
w_norm1= 1.8;
z=w_des/w_norm1
md=6;
Rl= (1+s/z)/(1+s/(md*z));
mi=10^(18.2/20);
Ri= 1/mi;
L= kcin/s^nu *Gp*Gs*Gf*Ga*Rl*Ri;

figure(1)
bode(Lin)

figure(2)
[num,den]= tfdata(Lin,'v')
nyquist1(num,den)

figure(3)
myngridst(Tpmax,Spmax)
hold on
nichols(Lin)
hold on
nichols(L)


%% now to check 
figure(6)
T= minreal(zpk(L/(L+1)));
S=minreal(zpk(1/(1+L)));
step(T*kd)
stepinfo(T*kd,"SettlingTimeThreshold",0.05)
figure(7)
bodemag(S)

figure(8)
bodemag(T)