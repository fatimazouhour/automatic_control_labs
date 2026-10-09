clear all
close all
clc

s=tf('s');

Gp= 100/(s^2+5.5*s+4.5);
[rp,pp]= tfdata(Gp,'v')

%given
Gs=1;
Ga=0.014;
Gr=1;
Gd=1;
Dao= 1.5*1e-3;

wp=0.03;
ws=60;

ap=16*1e-2;
as=2*1e-1;

kd=1;
Ro=1;

shat=12;
ts5=8;
tr=2;

Gf=1/Gs;
kp=dcgain(Gp);
%2
r1= 1.5*1e-1;
kc2=kd^2*Ro/(kp*Ga*r1)

%3
r2= 4.5*1e-3;


%4
r4=2*1e-3
Ms = 20*log10(r4/ap)

wl= 10^(-Ms/40)*wp;
wc_s= wl*2

%5
r5= 8*1e-4;
Mt= 20*log10(r5*Gs/as)

wh= 10^(Mt/40)*ws;
wc_t= wh/2
%
zeta= abs(log(shat))/sqrt(pi^2+(log(shat))^2),
Tpmax= 1/(2*zeta*sqrt(1-zeta^2))
Spmax= 2*zeta*sqrt(2+4*zeta^2+2*sqrt(1+8*zeta^2))/(sqrt(1+8*zeta^2)+4*zeta^2-1),

%tr cond
r_tr= 1/sqrt(1-zeta^2)*(pi-acos(zeta))*sqrt(sqrt(1+4*zeta^4)-2*zeta^2);
wc_tr= r_tr/tr

r_ts= -log(5)/zeta *sqrt(sqrt(1+4*zeta^4)-2*zeta^2);
wc_ts5= r_ts/ts5;

%
wc_min= max([wc_ts5,wc_tr,wc_s])
wc_max= wc_t

kc=22;
%initial contr
nu=1
Lin= kc/s^nu *Gp*Gs*Gf*Ga

figure(1)
bode(Lin),


wdes=1.2;
wnorm1= 0.75;
z=wdes/wnorm1,
Rz=(1+s/z);

wnorm2= 100;
pi= wdes/wnorm2,
mi= 10^(13/20),
Ri=(1+s/(mi*pi))/(1+s/pi);

L= kc/s^nu*Gp*Gs*Ga*Gf*Rz*Ri;

figure(2)
[num,den]= tfdata(Lin,'v');
nyquist1(num,den)

figure(3)

myngridst(Tpmax,Spmax)
hold on,
nichols(Lin)
hold on 
nichols(L)

T=minreal(zpk(L/(1+L)),1e-3);
S= minreal(zpk(1-T),1e-3);
kd=1;
t=linspace(0,20,50000);
figure(4)
step(kd*T)
r=stepinfo(kd*T,'SettlingTimeThreshold',0.05)

figure(5)
bodemag(T)
figure(6)
bodemag(S)