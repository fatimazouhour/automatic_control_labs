%% example slide 70
clear all
close all 
clc

s = tf('s');
H = 1/(s^2+s+1);
H.inputdelay = 2;
SYSX = zpk(pade(H,1))