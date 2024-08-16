clear all; clc; close all;
%--------------------------------------------------------------------------
% Load example file
%--------------------------------------------------------------------------
load examples\example1;

%--------------------------------------------------------------------------
% Set up parameters
%--------------------------------------------------------------------------
Setting.h_min          = 30;
Setting.h_max          = 2e2;
Setting.h0             = 1*Setting.h_min;
Setting.K              = 20;
Setting.g              = 0.5;
Setting.Smoothing      = true;
Setting.SmoothingParam = 1e-2;

%--------------------------------------------------------------------------
% Start mesh generation
%--------------------------------------------------------------------------
fixedPoints = [];
Mesh1D = MeshGeneration1D(Points,fixedPoints,Setting);

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
figure; plot(Points(:,1),Points(:,2),'k');
hold on; plot(Mesh1D.X,Mesh1D.Y,'-*b');
daspect([1 1 1]);

