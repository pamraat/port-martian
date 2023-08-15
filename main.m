    
%% Non-equilibrium (chemical) nozzle flow calculation
%
% Species under consideration are  $Ar$, $C$, $N$, $O$  $C_2 $N_2$ $O_2$ $CN$
% $CO$ $NO$ $CO_2$ $NCO$ $C+$ $O+$ $CO+$ $NO+$ $E-$
%
% Dependencies of the equilibrium solver:
% Executables: 1. cea. thermo.lib 3. trans.lib CEA executable and its
% libraries
%
% Dependencies of the non-equilibrium solver:
% Executables: 1. cea. thermo.lib 3. trans.lib 
% Data Files:  1. martian_kinetic_data.csv 2. martian_nasa.dat - kinetic database (Gupta, Park 1994)
% and NASA polynomial database
%
% Usage:
% Enter valid inpt data in INPUT.m
% Ensure all dependent files are in the current folder
% Results are recorded in the file results.dat
% ================================
% © K Satheesh and Devendra Ghate (Edited by Pamraat Parmar for Martian Atmosphere)
% ================================

clear all;
global R0 mdotstardash rhostardash Tstardash m0 rstar Theta theta l coefVec speciesVec speciesvec viscflag;
global alpha beta Kf Kb z;

INPUT

if strcmp(viscflag,'TRUE')
    dp02 = 1;
    theta = Theta;
    while dp02>1e-2
        constants
       % Equilibrium solver
        if strcmp(solmethod,'EQ')
            equilibriumSolver;
        elseif strcmp(solmethod,'NEQ')  
            nonequilibriumSolver;
        end
%% Free stream pitot pressure
        [rho2,u2,h0,p02,T02,rho_02,mu_02,k_02,cp_02,mf_02] = pitot(uinf,pinf,tinf,rhoinf,speciesvec,minf);%h0-Jkg^-1K^-1,p02-bar
        dp02 = abs((p02-P02)/P02);
        theta = theta-0.025;
        if p02>P02
            error('Switch off the ''viscflag'' option')
        end
    end  
else
    constants
    if strcmp(solmethod,'EQ')
        equilibriumSolver;
    elseif strcmp(solmethod,'NEQ')  
        nonequilibriumSolver;
    end
    
end

results
if strcmp(cleanflag,'TRUE')
    [cmd,op] = system('del *.inp *.plt *.out');
end

