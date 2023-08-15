%% calculates and prints the  relavent output data.
% Properties calculated here:
% 1.Pitot pressure
% 2.Stagnation point heat flux
% Formated output is generated as results.dat
% © K Satheesh & Devendra Ghate
[rho2,u2,h0,p02,T02,rho_02,mu_02,k_02,cp_02,mf_02] = pitot(uinf,pinf,tinf,rhoinf,speciesvec,minf);
[q] = fay_riddel(R,T_w,p02,rho_02,T02,h0,rhoinf,rho2,u2,cp_02,k_02,mf_02); 
if strcmp(solmethod, 'EQ')
    disp('Type of Solver: Equilibrium');
else
    disp('Type of Solver: Non-equilibrium');
end
disp(' ');
disp('Chemical kinetic model: 12 species Martian air (With ionisation)');
disp(' ');
disp(['Area ratio of the nozzle: ',num2str(A(1))]);
disp(' ');
disp('Experimental data provided by the user:');
disp(['    Incident shock velocity:          ',num2str(u1),' m/s']);
disp(['    Driven tube fill pressure:        ',num2str(p1),' bar']);
disp(['    Driven tube temperature:          ',num2str(T1),' K']);
disp(['    Measured reservoir pressure (p5): ',num2str(p5e),' bar']);
disp(' ');
if strcmp(solmethod, 'EQ')
    load eq_nozzle.plt 
    load eq_nozzle1.plt;
    rhodash = eq_nozzle(2:size(eq_nozzle,1),1);
    Tdash = eq_nozzle(2:size(eq_nozzle,1),2);
    udash = eq_nozzle(2:size(eq_nozzle,1),3).*eq_nozzle(2:size(eq_nozzle,1),4);  
    pdash = eq_nozzle(2:size(eq_nozzle,1),17);
    massf = [eq_nozzle(2:size(eq_nozzle,1),5:16) eq_nozzle(2:size(eq_nozzle1,1),:)];
    ar = [1 ar_sup]'; %nozzle area ratio
    loglog([1 ar_sup],eq_nozzle(2:size(eq_nozzle,1),17)/p5e,'g'); %plot pressure
    xlabel('Area Ratio'); ylabel('P/P_0');
    title('Static pressure in nozzle');
elseif strcmp(solmethod, 'NEQ')
    fprintf('Noneq simulation started %d%% downstream of throat.\n',xi*100);
    massf = y(:,3:19)*diag(molecularWeightVec); %mass fractions
    %% testing energy conservation 
    atom_O = y(:,9)*2+y(:,6)+y(:,11)+y(:,12)+y(:,13)*2+y(:,14)+y(:,16)+y(:,17)+y(:,18);
    atom_N = y(:,8)*2+y(:,5)+y(:,10)+y(:,12)+y(:,14)+y(:,18);
    [Cp, h] = getCph_T_mod(Tstardash, R0, coefVec, speciesVec);
    Hstardash = sum(y(1,3:19).*h)*R0*Tstardash;
    for i=1:size(Tdash,1)
        [Cp, h] = getCph_T_mod(Tdash(i), R0, coefVec, speciesVec);
        gamma = y(i,3:19);
        hdash = h*R0*y(i,2)*Tstardash; 
        Hdash = sum(gamma.*hdash);
        H0dash(i) = Hdash + udash(i)*udash(i)/2; 
    end
    err = (max(H0dash)-min(H0dash))/min(H0dash);
    disp(['Residual for the energy equation ',num2str(err)]);

    loglog((A(x)),pdash./p0dash,'r.');
    xlabel('Area Ratio'); ylabel('P/P_0');
    title('Static pressure in nozzle');
    %plot(A(x),Tdash/T5,'r')
end
if strcmp(viscflag,'TRUE')
    fprintf('Nozzle half angle changed to %d degree.\n',theta+dtheta);
end
disp(['Free-stream static pressure (pinf): ',num2str(pinf),' Pa']);
disp(['Free-stream static temperature (tinf): ',num2str(tinf),' K']);
disp(['Free-stream density  (rhoinf): ',num2str(rhoinf),' kg m^-3']);
disp(['Free-stream pitot pressure (p02): ',num2str(p02),' bar']);
disp(['Free-stream velocity (u): ',num2str(uinf),' msec^-1']);

disp(['Free-stream mass fractions' 10.....
    9 'Ar:' 9 ,num2str(minf(1)), 10 ...
    9 'C:' 9 ,num2str(minf(2)),10 ...
    9 'N:' 9 ,num2str(minf(3)), 10 ...
    9 'O:' 9 ,num2str(minf(4)),10 ...
    9 'C2:' 9 ,num2str(minf(5)), 10 ...
    9 'N2:' 9 ,num2str(minf(6)), 10 ...
    9 'O2:' 9 ,num2str(minf(7)), 10 ...
    9 'CN:' 9 ,num2str(minf(8)), 10 ...
    9 'CO:' 9 ,num2str(minf(9)), 10 ...
    9 'NO:' 9 ,num2str(minf(10)), 10 ...
    9 'CO2:' 9 ,num2str(minf(11)), 10 ...
    9 'NCO:' 9 ,num2str(minf(12)), 10 ...
    9 'C+:' 9 ,num2str(minf(13)), 10 ...
    9 'O+:' 9 ,num2str(minf(14)), 10 ...
    9 'CO+:' 9 ,num2str(minf(15)), 10 ...
    9 'NO+' 9 ,num2str(minf(16)), 10 ...
    9 'E-:'  9,num2str(minf(17))]);
disp(['Stagnation Enthalpy (h0): ',num2str(h0),' J/kg/K']);
disp(['Stagnation heat flux (q): ',num2str(q),' W/m^2']);
disp(' ');
%% Saving the data in results.dat
%u          save('solution.mat','x','udash','rhodash','P','Tdash','uinf','pinf','tinf','rhoinf','minf','p02','h0');
fid = fopen('results.dat','w');
fprintf(fid,'Shot number %s\n',shot);
if strcmp(solmethod, 'EQ')
    fprintf(fid,'%%Type of Solver: Equilibrium\n');
else
    fprintf(fid,'%%Type of Solver: Non-equilibrium\n');
end
fprintf(fid,'%%Chemical kinetic model: 12 species Martian air (With ionisation)\n');
fprintf(fid,'%%Area ratio of the nozzle: %d\n',A(1));
fprintf(fid,'%%Experimental data provided by the user\n');
fprintf(fid,'%% \t Incident shock velocity: %d ms^-1       \n',u1);
fprintf(fid,'%%\t Driven tube fill pressure:%d bar        \n',p1);
fprintf(fid,'%%\t Driven tube temperature:%d K          \n',T1);
fprintf(fid,'%%\t Reservoir pressure (p5):%d bar \n',p5e);
fprintf(fid,'%%\t Free stream pitot pressure (p02):%d bar \n',P02);
fprintf(fid,'%%\t Free stream stagnationpoint heat flux (q):%d Wm^-2 \n',Q);
fprintf(fid,'%%Nozzle inlet conditions\n');
fprintf(fid,'%% \t Pressure(P5): %d bar       \n',p5);
fprintf(fid,'%% \t Temperature(T5): %d K       \n',T5);
fprintf(fid,'%% \t Enthalpy(h5): %d kJ/kgK       \n',h5);
fprintf(fid,'%% \t Density(rho5): %d kg/m^3       \n',rho5);
massfr = num2str(mf5);
fprintf(fid,'%% \t Species Mass fraction [Ar C N O C2 N2 O2 CN CO NO CO2 NCO C+ O+ CO+ NO+ e-]: %s        \n',massfr);
fprintf(fid,'%% \t Gamma(gamma5): %d        \n',gamma5);
if strcmp(solmethod, 'NEQ')
    fprintf(fid,'%%Noneq simulation started %d%% downstream of throat.\n',xi*100);
end
if strcmp(viscflag,'TRUE')
    fprintf(fid,'%%Nozzle half angle changed to %d degree.\n',theta+.05);
end
fprintf(fid,'%%Concentrations listed in mass fractions\n');
OUT = [ar rhodash Tdash pdash udash massf];
fprintf(fid,'%12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n','%A/A*', 'rho(kg/m3)','T(K)','P(pa)','u(m/s)', 'Ar', 'C', 'N', 'O', 'C2', 'N2', 'O2', 'CN', 'CO', 'NO', 'CO2', 'NCO', 'C+', 'O+', 'CO+', 'NO+', 'e-');
fprintf(fid,'%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n',OUT');
fclose(fid);
disp('Output written to file results.dat');
disp('© K Satheesh & Devendra Ghate');
disp('Edited for martian atmosphere by Pamraat Parmar');
