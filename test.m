clear all; close all;
addpath(genpath('../QIC-IDADE'));
N = 3; % number of ions in the trap
w_tar = ones(1,N); % Target eigenfrequencies
B_init = ones(N); % Initial guess for the eigenvectors
A_conv = ones(N)./2; % Matrix for the trap w/o tweezers
B_ens = cat(3, B_init, B_init-2, B_init+1,B_init+4); % create the ensemble of B matrices

[out] = full_IDA(w_tar,B_ens,A_conv,0.01,100);

DEparams = [0.5, 0.1, 0.3];
A_DE = DE(out.A_bars,w_tar,A_conv,DEparams);


