%% Equilibrium solver solves equilibrium equations using NASA CEA code.
% The executable for the program "cea_run" should be present in the
% directory along with the requisite library files.
format short g
ar_sup = logspace (0.04,log10(A(1)),9);
format
if strcmp(viscflag,'TRUE')
    x = (sqrt(ar_sup*Astar/pi)-rstar)./tan(theta*pi/180);
else
    x = (sqrt(ar_sup*Astar/pi)-rstar)./tan(Theta*pi/180);
end
[p5,T5,h5,rho5,gamma5,mf5,rhoinf, tinf, uinf, minf,pinf] =  cea_run(u1,p1,T1,speciesvec,mf1, ar_sup,p5e); %cea eq solution

