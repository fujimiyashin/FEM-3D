clc, clear;
close all;

G = 80;
nu = 0.3;

dx = 1;
dy = 1;
dz = 1;

len_x = 2;
len_y = 10;
len_z = 2;

mx = len_x*dx;
my = len_y*dy;
mz = len_z*dz;

node_sum = (len_x+1)*(len_y+1)*(len_z+1);

unit_stiffness_p = zeros(24,1);
unit_stiffness_k = zeros(24,24);
stiffness_k = zeros(3*node_sum, 3*node_sum);
stiffness_p = zeros(3*node_sum, 1);

node_lis = [];
cubic_lis = [];
color = "#a0e0f0";
magnitude = 1;