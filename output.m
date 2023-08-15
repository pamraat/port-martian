% This script post-processes the solution and calculates the
% relavent output conditions.
% A log-log plot of "P/P0" Vs "A/Astar"
%
% © K Satheesh & Devendra Ghate

if strcmp(solmethod, 'EQ')
    pinf = pinf*1e5; %Pa
    udash = shock(2:size(shock,1),3).*shock(2:size(shock,1),4);%M*a
    rhodash = shock(2:size(shock,1),1);
    P = shock(2:size(shock,1),11);
    Tdash = shock(2:size(shock,1),2);
    %loglog([1 ar_sup],shock(2:12,2)/T5,'r.'); %plot temp
    system('mv shock.plt eq_nozzle.plt');
elseif strcmp(solmethod, 'NEQ')
    molewt = 1./sum(y(:,3:8),2); %kg/mole
    Tdash = y(:,2)*Tstardash;
    P = rhodash.*R0.*Tdash./molewt; %pa
    atom_O = y(:,3)*2+y(:,5)+y(:,7);
    atom_N = y(:,4)*2+y(:,6)+y(:,7);
    [Cp, h] = getCph_T_mod(Tstardash, R0, coefVec, speciesVec);
    Hstardash = sum(y(1,3:8).*h)*R0*Tstardash;
       %% testing energy conservation
    for i=1:size(Tdash,1)
        [Cp, h] = getCph_T_mod(Tdash(i), R0, coefVec, speciesVec);
        gamma = y(i,3:8);
        hdash = h*R0*y(i,2)*Tstardash; 
        Hdash = sum(gamma.*hdash);
        H0dash(i) = Hdash + udash(i)*udash(i)/2; 
    end
    err = (max(H0dash)-min(H0dash))/min(H0dash);
      %% free stream properties
    pinf = min(P); %pa
    uinf = max(udash); %m/sec
    tinf = min(Tdash); %K
    rhoinf = min(rhodash); %kg m^-3
    minf = y(size(y,1),3:8).*molecularWeightVec; % composition in mas frac
       
    %plot(A(x),Tdash/T5,'r')
end

%% Free stream pitot pressure
[rho2,u2,h0,p02,T02,rho_02,mu_02,k_02,cp_02,mf_02] = pitot(uinf,pinf,tinf,rhoinf,speciesvec,minf); %Jkg^-1K^-1,bar

% [Cp, h] = getCph_T_mod(Tstardash, R0, coefVec, speciesVec);
% hdash = h*R0*Tstardash; 
% Hdash = sum(gamma0.*hdash);
% H0dash = Hdash + ustardash*ustardash/2
% mdot = rhodash.*udash.*A(x);
%plot(t,gamma_v(:,1)*R0*Td/1e8,'*');
%figure;
%plot(t,gamma_v(:,3)*R0*Td/10^8,'o');

%molefr=gamma_v*R0*Td/Pd %output in mole fractions

%m = 1/sum(gamma);
% Stagnation heat flux
[q] = fay_riddel(R,T_w,p02,rho_02,T02,h0,rhoinf,rho2,u2,cp_02,k_02,mf_02); 

