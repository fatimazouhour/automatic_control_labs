% lecture 10 example pdf  
clear all
close all
clc

format default
% define what we have given
m=0.1;
l=0.5;
b =0.1;
g=9.81;
u_eq=0;
k=0;
x_eq=[k*pi;0]; %from i/o equations
x0=x_eq +[5/180*pi;0];

% tsim=10;
% out=sim('sim_pendulum_20');

% ask about the pendulum file