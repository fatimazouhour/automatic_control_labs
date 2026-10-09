%%
clear all
close all
clc

format default
s = tf("s");
H = 1/((s+2)*(s+10));
k = dcgain(H)

% el dcgain bta2teene el lim (s-->0) of s^r * H(s)
