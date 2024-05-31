%%% This script describes how the spokes ptx pulses can be designed 
%%% and an older version of it was used for creating the 7T liver imaging results
%%% reported in X. Wu et al. Quantitative Imaging in Medicine and Surgery
%%% 2014.
% Description: 
% Author: Xiaoping Wu <xpwu@cmrr.umn.edu> 
% Maintainer: 
clear
addpath('./src/');
addpath('./regtools/')

%
load fieldmaps.mat
% pulse design specs
Specs.mask = mask; % region of interest where the pulses
Specs.mask0= mask0;

Specs.b1maps = b1maps;
Specs.b0map= b0map;
Specs.fox= 1e-3*[450 450];% field of excitation in m used for field mapping
Specs.dt = 10e-6; % pulse dwell time in sec
Specs.nominalfox = 8:2:18; % a range of fox in cm over which the kspace placement is optimized.
Specs.subrfduration=1; % sub pulse duration in ms
Specs.numofspokes=2;
%-------------------------------

rfDesigner= pTX3DPulseDesigner(Specs.b1maps,Specs.mask,1e2*Specs.fox);
rfDesigner.DwellTime=1e6* Specs.dt;
rfDesigner.NominalFOX= Specs.nominalfox;
rfDesigner.SubRFDuration= Specs.subrfduration;
rfDesigner.B0Map= Specs.b0map;

rfDesigner.NumOfSpokes= Specs.numofspokes;

[rfs,grad,kp,myrho]= rfDesigner.design;


% run Bloch simulation to examine in-plane excitation fidelity
mxypat = run_bloch_sim (rfs,grad(1:2,:),Specs.b1maps,Specs.mask0,...
    Specs.fox, Specs.b0map,0,[],Specs.dt);
figure, myimagesc((abs(mxypat)),Specs.mask), colormap jet,axis square

% run Bloch simulation with high resolution field mapping.
load FieldmapsHR.mat
MxyPat = run_bloch_sim (rfs,grad(1:2,:),b1Maps,Mask0,...
    Specs.fox, b0Map,0,[],Specs.dt);
figure, myimagesc((abs(MxyPat))), colormap jet,axis square

%%% now you are ready to write your pulses into a text file for
%%% your sequence to load in and apply them. 

