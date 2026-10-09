%% 
clear all
close all
clc
format compact 
s= tf('s');
%part 1 on paper

%2
A=[0 1;-3.5 -0.75]
eig(A)
Aprime=[0 1;-3.5 0] 
eig(Aprime)
%%
clear all
close all
clc

format compact
m = 0.1;l=0.5;beta=0.1;g=9.81;
u =m*g*l*sin((2*pi)/3)
x1 = 2*pi/3;x2 =0;
A= [0 1;(-g/l*cos(x1))  -beta/(m*l^2)];
B = [0;1/(m*l^2)]