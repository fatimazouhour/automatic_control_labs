%% Problem 1
clear all
close all
clc

format default

A= [1 2;4 3];
B = [5 ;8];
C=[-1 3];
D=[8];

x0 = [2;2];
s = tf('s');
U=2/s;
x =zpk(minreal(inv((s*eye(2))-A)*(x0+B*U),1e-3));
y = zpk(minreal(C*x+D*U,1e-3));
% to get the equations for x

[num_x,den_x] = tfdata(x(1),'v');
[rx,px]= residue(num_x,den_x)

[num_x2,den_x2] = tfdata(x(2),'v');
[rx2,px2]= residue(num_x2,den_x2)

%now for y
[num_y, den_y] = tfdata(y,'v');
[ry, py] = residue(num_y,den_y)

%% Problem 2
clear all
close all
clc

A = [0 1; -1 -1];
B=[4 0;10 1];
C=[1 0];


s= tf('s');
U = [0 ;2];
x0 = [0 ;0];

x = zpk(minreal(inv(s*eye(2)-A)*(x0+B*U),1e-3));
y = zpk(minreal(C*inv(s*eye(2)-A)*(x0+B*U),1e-3));

[numy,deny] = tfdata(y(1),'v');
[ry,py] = residue(numy,deny)

M= abs(ry(1)) , 2*M
ang = angle(ry(1))

%% Problem 3

clear all
close all
clc

A = [0 6;-1 -5];
B = [0 1]';

s= tf('s');
X= zpk(minreal(inv(s*eye(2)-A),1e-3));

[n1,d1] = tfdata(X(1,1),'v');
[r1,p1] = residue(n1,d1)
[n2,d2] = tfdata(X(1,2),'v');
[r2,p2] = residue(n2,d2)
[n3,d3] = tfdata(X(2,1),'v');
[r3,p3] = residue(n3,d3)
[n4,d4] = tfdata(X(2,2),'v');
[r4,p4] = residue(n4,d4)
