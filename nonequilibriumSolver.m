
% Non-equilibrium solver solves 12 species (martian atmosphere) chemical non-equilibrium
% equations. The solver is initialised using the equilibrium solution
% from the NASA CEA program (cea_run) and the solution process starts at
% a location downstream of the throat.
% Ref: Lordi,Mates,Moselle, NASA CR-472, 1966

du = 2e-4; % check for supersonic solution in nozzle
ar_sup = 1.0001; %supersonic area ratio

counter = 0;
while du>1e-4 || id == 1
    id = 0;
    if strcmp(viscflag,'TRUE')
        xi = (sqrt(ar_sup)-1)*rstar/(tan(theta*pi/180))/l; %non-eq simulation starts from this x
    else
         xi = (sqrt(ar_sup)-1)*rstar/(tan(Theta*pi/180))/l;
    end
    [p5,T5,h5,rho5,gamma5,mf5,rhostardash, Tstardash, ustardash, massfrac,pstardash] = cea_run(u1,p1,T1,speciesvec,mf1, ar_sup,p5e); %run cea for initial conditions
   
    if abs(sum(massfrac)-1)>1e-4
        error('component mass fractions do not sum to 1.0');  
    end
%     massfrac(17)=0;
    gamma0 = massfrac./molecularWeightVec;  % massfrac/molewt, mole/kg 
    m0 = 1/sum(gamma0); % Mixture Molecular weight,kg/mole
    pstardash = rhostardash*R0*Tstardash/m0; %pressure, Pa
    p0dash = p5e*1e5;%pa
    % H0 = (hstardash+ustardash*ustardash/2)/R0/Tstardash*m0;
    mdotstardash = rhostardash*Astar*ustardash; %mass flow rate kg/sec

    %% Solving the EOM

    y0 = [1, 1, gamma0]; 
    options = odeset('Mass',@mass,'RelTol',1e-6 ,'stats','off');
    [x, y] = ode23tb(@rhs, [xi 1], y0, options);
    %[LASTMSG, LASTID] = lastwarn;
    %if strcmp(LASTID,'MATLAB:ode23tb:IntegrationTolNotMet')

    if max(x)<1 %when simulation stops before nozzle exit
        ar_sup = ar_sup+ AREA_INC_STEP; 
        id=1;
    end
    rhodash = y(:,1)*rhostardash;
    udash = mdotstardash./(rhodash.*feval(A,x)*Astar);
    du = ustardash-udash(end) ;
    ar_sup = ar_sup+ AREA_INC_STEP;
    counter = counter + 1;
end
%% Results
ar = A(x); %nozzle area ratio
molewt = 1./sum(y(:,3:2+length(speciesVec)),2); %kg/mole
Tdash = y(:,2)*Tstardash;
pdash = rhodash.*R0.*Tdash./molewt; %pa
pinf = min(pdash); %pa
uinf = max(udash); %m/sec
tinf = min(Tdash); %K
rhoinf = min(rhodash); %kg m^-3
minf = y(size(y,1),3:2+length(speciesVec)).*molecularWeightVec; % composition in mas frac