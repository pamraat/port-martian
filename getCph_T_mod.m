
%% Computes Cp and h using the NASA polynomial coefficients obtained fron readsandia.m
% From aerosfrosh code by Mattheww Campbell, Stanford
function [Cp, h, h_0] = getCph_T_mod(T, Ru, coefVec, speciesVec)
% This finds the specific heat and enthalpy at the specified temperature
% loop over all species to compute for each species first
    HoverRT = zeros(1,size(coefVec,3));
    H_0 = zeros(1,size(coefVec,3)); %h at 0K
    CpoverR = zeros(1,size(coefVec,3));
    for i=1:size(coefVec,3) % loop over 3rd dimension of coefVec 
        % grab temperature limits
            Tlow = coefVec(3,1,i);
            Tcom = coefVec(3,2,i);
            Thigh = coefVec(3,3,i);
        % check temperature limits
            if T < Tlow
                %fprintf('%3s: T (%1.2f K) < Tlow (%1.0f K)\n',speciesVec{i},T,Tlow);
            elseif T > Thigh
                %fprintf('%3s: T (%1.2f K) > Thigh (%1.0f K)\n',speciesVec{i},T,Thigh);
            end
        % decide which set of data to use
            if T < Tcom % low temp set
                j = 2;
            else % high temp set for T >= Tcom
                j = 1;
            end
        % grab temperature coefficients
            a1 = coefVec(j,1,i);
            a2 = coefVec(j,2,i);
            a3 = coefVec(j,3,i);
            a4 = coefVec(j,4,i);
            a5 = coefVec(j,5,i);
            a6 = coefVec(j,6,i);
            a7 = coefVec(j,7,i);
        % compute Cp/R
            CpoverR(i) = a1 + a2*T + a3*(T^2) + a4*(T^3) + a5*(T^4); %Cp[j/mole/K]/Ro
        % compute h/RT
            HoverRT(i) = a1 + a2*T/2 + a3*(T^2)/3 + a4*(T^3)/4 + a5*(T^4)/5 + a6/T; %h[j/mole]/(Ro*T)
           H_0(i) = coefVec(j,6,i)*Ru;
    end
% now compute mixture properties
    Cp = CpoverR;
    h = HoverRT;
    h_0 = H_0; %h@ T=0K
    %Cp = (X*CpoverR')*Ru*(1000/gasMolWt); % (no units)*(J/mol*K)*(1000 g/kg)*(mol/g) = J/kg*K
    %h = (X*HoverRT')*T*Ru*(1000/gasMolWt); % J/kg
    
    
