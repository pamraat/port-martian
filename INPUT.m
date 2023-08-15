%% Input conditions for the problem to be entered here
%% Problem ID
shot = 'R625'; %Enter the exp run numner
solmethod = 'NEQ'; % 'EQ' or 'NEQ'
viscflag = 'aTRUE'; % accounts for boundary layer by changing theta
cleanflag = 'TRUE'; % cleans unwanted CEA generated files 
%% Defining the Nozzle
% The nozzle definition assumes a conical nozzle. The function defining
% the area variation is defined in "constants.m".

% VSSC nozzle
rstar = 12.5e-3; %throat radius, m
rexit = 500e-3; %exit radius, m 
Theta = 10; %nozzle halfangle, deg
R =15e-3; %Nose radius of heat flux probe, m
T_w = 300 ; %Wall temp of hat flux probe,K
% T5 nozzle
%   rstar = 15.0e-3; %throat radius, m
%   rexit = 150e-3; %exit radius, m 
%   theta = 7; %nozzle halfangle, deg

%% Initial conditions (in shock tube) 
u1 = 3192; %incident shock veocity, m/s
p1 = .06; %driven tube fill pressure, bar
T1 = 300; %%driven tube  temperature, K

p5e = 52.86; % measured p5, bar
P02 = 0.0656; % measured pitot pressure, bar
Q = 230; % stagnation point heat flux, W/cm^2






