%% lab 3 mainly related with subjects discussed in L06 and L067
%% Problem 1

clear all
close all
clc

A= [0 1;0.1 -0.3];
B= [0 1]';
C= [7 7];
D=0;

z=tf('z',1);
x0=[0 0]';
u=z/(z-1);

x=zpk(minreal(inv(z*eye(2)-A)*(z*x0+B*u),1e-3));
[numx1,denx1]=tfdata(x(1),'v');
[rx1,px1,kx1] =residuez(numx1,denx1)

[numx2,denx2]=tfdata(x(2),'v');
[rx2,px2,kx2] =residuez(numx2,denx2)

%now for y
y=zpk(minreal(C*inv(z*eye(2)-A)*(z*x0+B*u),1e-3));
[numy,deny]=tfdata(y,'v');
[ry,py,ky]=residuez(numy,deny)

% for the modal analysis
lambda=eig(A)

sys=ss(A,B,C,D,1);
H=tf(sys)
[numh,denh]=tfdata(H,'v');
[r,p,k]=residuez(numh,denh)

%% p2
clear all
close all
clc

s=tf('s');
h=1/(s^3+2*s^2+5.25*s+4.25);
wo=0.1;
[m,f]= bode(h,wo)
y_hat=m*3
f/180*pi
ys2 = dcgain(h)

%part 3
%on paper

%% problem 2
clear all
close all
clc

s = tf('s');
h = 10/(s^2+1.6*s+4)
ts = linspace(0,20,10000);
[y,t,x]=step(h,ts);

figure(1)
plot(t,y,'b','linew',1.5)
grid on, zoom on, hold on, xlabel('t(s)'),
ylabel('y(t)')
%% problem 3

%on paper 
%% problem 4
clear all
close all
clc

%% problem 5 
clear all
close all
clc
s= tf('s');
A=[0 1 0 0;-1 0 0 0;3 0 0 0;0 0 2 0];
B=[1 0 0 0]';
C=[1 0 0 0];
D=0;
x0= [0 0 0 0]';
u = sqrt(2)/(s^2+2)

y = zpk(minreal((C*inv(s*eye(4)-A))*(x0+B*u),1e-3));
[numy,deny] = tfdata(y,'v');
[r,p]= residue(numy,deny)

M = abs(r), 2*M
theta = angle(r)