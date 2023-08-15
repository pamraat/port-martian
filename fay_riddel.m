%% Calculates Fay & Riddel stagnation point heat flux 
%Ref: Olivier, Shock waves, 3, 1999
function [q] = fay_riddel(R,T_w,p02,rho_02,T02,h0,rhoinf,rho2,u2,cp_02,k_02,mf_02)  %Olivier, Shock Waves (1993) 3:129-139
    global R0  coefVec speciesVec
    a = 1.458e-6; b = 110.4; %Sutehrland constants
    [Cp, h, h_0] = getCph_T_mod(298.15, R0, coefVec, speciesVec);
     h_f298 = h*R0*298.15;% h_T = h_f298 + (h_T-h_298)
     %estimation of h_f0 from h_f298 Ref: Gupta
     h_f0(1) = 0;
     h_f0(2) = 0;
     h_f0(3) = h_f298(:,3)-(h_f298(:,3)-h_0(:,3))+.5*(h_f298(:,1)-h_0(:,1)); %O
     h_f0(4) = h_f298(:,4)-(h_f298(:,4)-h_0(:,4))+.5*(h_f298(:,2)-h_0(:,2)); %N
     h_f0(5) = h_f298(:,5)-(h_f298(:,5)-h_0(:,5))+.5*(h_f298(:,1)-h_0(:,1)+h_f298(:,2)-h_0(:,2)); %NO
     h_D = mf_02(1)*h_f0(3)+mf_02(2)*h_f0(4)+mf_02(3)*h_f0(5);

     %Assuming Le = 1
     rho_w = p02*1e5/(287*T_w);
     mu = @(T) a*T^1.5/(T+b); %Sutherland formula
     mu_w = mu(T_w);
     mu_s = mu(T02);
     rho_s = rho_02;
     h_s = h0;
     h_w = 1.004e3*T_w; %J/kg
     delta = R*rhoinf/rho2;
     b = delta+R;
     dvdy = (u2/R)*(1+(2+(b/R)^3)/2/((b/R)^3-1));
     %dvdy = (2*p02*1e5/rho_02)^.5/R;
     q = .94*(rho_w*mu_w)^.1*(rho_s*mu_s)^.4*(h_s-h_w)*(dvdy)^.5; 
     %Accounting for mass diffusion
     %Diffussion co-efficent, Gupta NASA RP1232 (assuming N2-N system instead of air)
     %T02  = 300;p02=1;cp_02=1e3;k_02=.02;rho_02 = 1.2
     a_d = 0;b_d = 0.0195;c_d = 1.488;d_d = -10.365; %curve fit co-efficients
     D = exp(d_d)*T02^(a_d*(log(T02))^2+b_d*log(T02)+c_d); %cm^2 atm sec^-1
     D = D*1.013*1e-4/p02; % m^2 sec^-1
     L_02 = D*rho_02*cp_02/k_02; %Le at 02

    %  T02  = 300;cp_02=1e3;k_02=.01;rho_02 = p02/R0/T02; %note k_02 os not correct
    % 
    %  D = exp(d_d)*T02^(a_d*(log(T02))^2+b_d*log(T02)+c_d); %cm^2 atm sec^-1
    %  D = D*1.013*1e-4/p02; % m^2 sec^-1
    %  L_w = D*rho_02*cp_02/k_02; %Le at wall
    %  Le = 0.5*(L_w+L_02);
    %  q = q*(1+  (L_02^0.52-1)*h_D/h_s);   