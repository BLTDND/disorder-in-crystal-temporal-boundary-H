format long
syms k

params.C_before =4.4;
params.C_after  = -7;
params.N_t=100;
params.kappa=2.3; 
params.k0=pi/3;%the central carrier wave frequency
params.N=200;
params.disorder_level=12;
params.repeat_times=30;

[R]=cal_disorder(params);






















