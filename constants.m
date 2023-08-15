%% Definitions of contants/parameters for use in the problem
%% Definition of nozzle area profile
Astar = pi*rstar*rstar;        % Throat area (m^2)
l = (rexit-rstar)/tan(Theta*pi/180);   % Nozzle length in meters
if strcmp(viscflag,'TRUE')
    A = @(x) pi*((tan(theta*pi/180)*x*l + rstar).^2)/Astar;    
else
    A = @(x) pi*((tan(Theta*pi/180)*x*l + rstar).^2)/Astar;    
end
 

%% Initialising Constants

R0 = 8.314; % Universal Gas Constant (J/mole/K)

% Assuming Martian air as the test gas.
mf1 = [1.93, 0.0, 0.0, 0.0, 0.0, 1.89, 0.146, 0.0, 0.0557, ...
    0.0, 95.97, 0.0]*1e-2; % test gas mass fraction
% mf1 = [1.28, 0.0, 0.0, 0.0, 0.0, 75.47, 23.2, 0.0, 0.0, ...
%     0.0, 0, 0.0]*1e-2; % test gas mass fraction

% 17 species chemical kinetic model of Martian air
% Source: 
speciesVec = {'AR',	'C', 'N', 'O', 'C2', 'N2', 'O2', 'CN', 'CO', 'NO', 'CO2', 'NCO', 'C+', 'O+', 'CO+', 'NO+', 'E-'}; % For Martian Atmosphere
% speciesvec = {'O2', 'N2', 'O', 'N', 'NO','Ar'}; %format for cea
speciesvec = {'Ar',	'C', 'N', 'O', 'C2', 'N2', 'O2', 'CN', 'CO', 'NO', 'CO2', 'NCO'}; %format for cea % For Martian Atmosphere
% molecularWeightVec = [15.9994*2, 14.0067*2, 15.9994, 14.0067, ...
% 14.0067+15.99994 39.948]*1e-3;  %kg/mole  
molecularWeightVec = [39.948, 12.01, 14.0067, 15.9994, 12.01*2, 14.0067*2, ...
15.9994*2, 12.01+14.0067, 12.01+15.9994, 14.0067+15.9994, 12.01+15.9994*2, ...
14.0067+12.01+15.9994, 12.01, 15.9994, 12.01+15.9994, 14.0067+15.9994, 5.4858e-4]*1e-3;  %kg/mole 

%% Reading species properties
%
% Species are numbered like this.
% $1~-~Ar$, $2~-~C$, $3~-~N$, $4~-~O$  $5~-~C_2 $6~-~N_2$ $7~-~O_2$ $8~-~CN$
% $9~-~CO$ $10~-~NO$ $11~-~CO_2$ $12~-~NCO$ $13~-~C+$ $14~-~O+$ $15~-~CO+$
% $16~-~NO+$ $17~-~e-$

coefVec = readSandia(speciesVec,'martian_nasa.dat'); %to be used in getCp_T_mod.m

% Reading Chemical kinetic data (Gupta: NASA RP 1232, 1990)
% Species are numbered like this.
% $1~-~Ar$, $2~-~C$, $3~-~N$, $4~-~O$  $5~-~C_2 $6~-~N_2$ $7~-~O_2$ $8~-~CN$
% $9~-~CO$ $10~-~NO$ $11~-~CO_2$ $12~-~NCO$ $13~-~C+$ $14~-~O+$ $15~-~CO+$
% $16~-~NO+$ $12~-~e-$
tmp = csvread('martian_kinetic_data.csv', 1, 1, [1 1 29 70]);
alpha = tmp(:,1:24);
beta = tmp(:,25:48);
z = tmp(1:8,49:64);
Kf = tmp(:,65:67); % A_{f,r}, B_{f,r} and T_{D_{f,r}}
Kb = tmp(:,68:70);

AREA_INC_STEP = 5.0e-5*4;