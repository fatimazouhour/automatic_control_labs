%%
clear all
close all
clc
format compact
s =tf('s');
u = 1/s;
H=(1+s/2)/(s*(1+s/10)^2)
Y = zpk(minreal(H*u,1e-3))
[numy,deny] = tfdata(Y,'v');
[r,p,k]=residue(numy,deny)
bode(H)
%%
clear all
close all
clc
format compact 
s = tf('s');
H = 10*(1+s/5)/((s+1)*(1+s/20)^2);
U = 1/s + 300/(s^2+300^2);
[m,f] = bode(H,300)
bode(H)
%%
clear all
close all
clc
format compact
s = tf('s');
H = 10*(1-s/5)/((s+1)*(1-s/20));
bode(H)
