%% Calculates stanation properties behind the normal shock at free stream conditions(02-pitot)
% Uses NASA CEA to isentropicaly decelerate the gas behind normal shock
%function [p5,T5,h5,rhostardash, Tstardash, ustardash, mf5, massfrac] = pitot(u1,p1,T1,tgas,mf1,rhodash,udash)
function [rho2,u2,h0,p02,T02,rho_02,mu_02,k_02,cp_02,mf_02] = pitot(u1,p1,T1,rho1,tgas,mf1)
    %% Create .inp file for cea shock tube run to compute properties behind normal shock at free stream conditions
    fileID = fopen('shock.inp','w');
    fprintf(fileID, '\n reac');
    for j=1:12
        fprintf(fileID,' name=%s  weight=%d t(k) %d\n',tgas{j},mf1(j),T1);
    end
    fprintf(fileID, 'prob ions case=3  p,bar=%d shock incd eq   u1=%d  \n output massf trace=1.e-12 \n only Ar C N O C2 N2 O2 CN CO NO CO2 NCO C+ O+ CO+ NO+ e- \n end',p1/1e5,u1);
    fclose(fileID);
    %% Run ./cea with shock.inp as input file
    [cmd,op] = system('cea'); 
    %% Read shock tube conditions from cea output 
    content = fileread( 'shock.out' ) ;
%     strexp = regexp(content,'U2.+\*C2.+\*N', 'once'); %check if C2 is present among the products

%     if isempty(strexp)
%         a  = regexp(content,...
%         'U2\,[^0-9]+([^a-zA-Z]+).+P\,[^0-9]+([^a-zA-Z]+).+T\,[^0-9]+([^a-zA-Z]+).+RHO\,[^0-9]+([^a-zA-Z]+).+H\,[^0-9]+([^a-zA-Z]+).+S\,[^0-9]+([^a-zA-Z]+).+GAMMAs[^0-9]+([\S]+).+\*e-[^0-9]+([\S]+).+\*Ar[^0-9]+([\S]+).+\*C[^0-9]+([\S]+).+\*CN[^0-9]+([\S]+).+\*CO[^0-9]+([\S]+).+\*CO+[^0-9]+([\S]+).+\*CO2[^0-9]+([\S]+).+\*N[^0-9]+([\S]+).+\NCO[^0-9]+([\S]+).+\*NO[^0-9]+([\S]+).+\*NO+[^0-9]+([\S]+).+\*N2[^0-9]+([\S]+).+\*O[^0-9]+([\S]+).+\*O+[^0-9]+([\S]+).+\*O2[^0-9]+([\S]+)'...
%         ,'tokens');%looking for string 'U2,', skips all char till 'P,', skips all non-digit, token () is all non alphabet char ..
%         a = a{:};
%         %a1 = {a{1:15} '0' a{16:23}};
%         a1 = {a{1:10} '0' a{10:14} '0' a{15:22}};
%     else
%         a  = regexp(content,...
%         'U2\,[^0-9]+([^a-zA-Z]+).+P\,[^0-9]+([^a-zA-Z]+).+T\,[^0-9]+([^a-zA-Z]+).+RHO\,[^0-9]+([^a-zA-Z]+).+H\,[^0-9]+([^a-zA-Z]+).+S\,[^0-9]+([^a-zA-Z]+).+GAMMAs[^0-9]+([\S]+).+\*e-[^0-9]+([\S]+).+\*Ar[^0-9]+([\S]+).+\*C[^0-9]+([\S]+).+\*C+[^0-9]+([\S]+).+\*CN[^0-9]+([\S]+).+\*CO[^0-9]+([\S]+).+\*CO+[^0-9]+([\S]+).+\*CO2[^0-9]+([\S]+).+\*C2[^0-9]+([\S]+).+\*N[^0-9]+([\S]+).+\NCO[^0-9]+([\S]+).+\*NO[^0-9]+([\S]+).+\*NO+[^0-9]+([\S]+).+\*N2[^0-9]+([\S]+).+\*O[^0-9]+([\S]+).+\*O+[^0-9]+([\S]+).+\*O2[^0-9]+([\S]+)'...
%         ,'tokens');%looking for string 'U5,', skips all char till 'P,', skips all non-digit, token () is all non alphabet char ..
%         a1 = a{:};
%     end
        a  = regexp(content,...
        'U2\,[^0-9]+([^a-zA-Z]+).+P\,[^0-9]+([^a-zA-Z]+).+T\,[^0-9]+([^a-zA-Z]+).+RHO\,[^0-9]+([^a-zA-Z]+).+H\,[^0-9]+([^a-zA-Z]+).+S\,[^0-9]+([^a-zA-Z]+).+GAMMAs[^0-9]+([\S]+).+\*Ar[^0-9]+([\S]+).+\*CO[^0-9]+([\S]+).+\*CO2[^0-9]+([\S]+).+\*N[^0-9]+([\S]+).+\NCO[^0-9]+([\S]+).+\*NO[^0-9]+([\S]+).+\*NO+[^0-9]+([\S]+).+\*N2[^0-9]+([\S]+).+\*O[^0-9]+([\S]+).+\*O2[^0-9]+([\S]+)'...
        ,'tokens');%looking for string 'U2,', skips all char till 'P,', skips all non-digit, token () is all non alphabet char ..
        a = a{:};
        %a1 = {a{1:15} '0' a{16:23}};
        a1 = {a{1:7} '0' a{8} '0' '0' '0' a{9} '0' a{10} '0'  a{11:16} '0' a{17}};
%     strexp = regexp(content,'U2.+\*N.+\*NO', 'once'); %check if N is present among the products
%     if isempty(strexp)
%         a  = regexp(content,...
%         'U2\,[^0-9]+([^a-zA-Z]+).+P\,[^0-9]+([^a-zA-Z]+).+T\,[^0-9]+([^a-zA-Z]+).+RHO\,[^0-9]+([^a-zA-Z]+).+H\,[^0-9]+([^a-zA-Z]+).+S\,[^0-9]+([^a-zA-Z]+).+GAMMAs[^0-9]+([\S]+).+\*Ar[^0-9]+([\S]+).+\*NO[^0-9]+([\S]+).+\*N2[^0-9]+([\S]+).+\*O[^0-9]+([\S]+).+\*O2[^0-9]+([\S]+)'...
%         ,'tokens');%looking for string 'U2,', skips all char till 'P,', skips all non-digit, token () is all non alphabet char ..
%         a = a{:};
%         a1 = {a{1:6} '0' a{7:10}};
%     else
%         a  = regexp(content,...
%         'U2\,[^0-9]+([^a-zA-Z]+).+P\,[^0-9]+([^a-zA-Z]+).+T\,[^0-9]+([^a-zA-Z]+).+RHO\,[^0-9]+([^a-zA-Z]+).+H\,[^0-9]+([^a-zA-Z]+).+S\,[^0-9]+([^a-zA-Z]+).+GAMMAs[^0-9]+([\S]+).+\*Ar[^0-9]+([\S]+).+\*N[^0-9]+([\S]+).+\*NO[^0-9]+([\S]+).+\*N2[^0-9]+([\S]+).+\*O[^0-9]+([\S]+).+\*O2[^0-9]+([\S]+)'...
%         ,'tokens');%looking for string 'U5,', skips all char till 'P,', skips all non-digit, token () is all non alphabet char ..
%         a1 = a{:}; 
%     end

    for i=1:size(a1,2)
        a2(i) = cellstr(regexprep(char(a1(i)),'([0-9])-','$1e-')); %replacing '-' by 'e-'
        a2(i) = cellstr(regexprep(char(a2(i)),'([0-9])\s([0-9])','$1e+$2')); %replacing '\s ' by 'e+'
    end
    %A=zeros(size(a2,1),size(a2,2));
    A = str2double(a2); %thermodynamic properties
    [u2, p2, T2, rho2, h2,S2, cp_cv_2, mf2(1:17)] = deal(A(1,1), A(1,2), A(1,3) ,A(1,4), A(1,5),A(1,6),A(1,7), A(8:24)); %shock tube stagnation conditions (bar,K,kg/m3,kJ/kgK) 
    mf2 = [ mf2(2) mf2(3) mf2(10) mf2(15) mf2(9) mf2(14) mf2(17) mf2(5) mf2(6) mf2(12) mf2(8)...
        mf2(11) mf2(4) mf2(16) mf2(7) mf2(13) mf2(1)];
    system('move shock.out normalshock_free_stream.out');
    %% Estimating p02 with cea sp run
    h0 = h2*1e3+0.5*u2*u2;
    %tgas = {'Ar',	'C', 'N', 'O', 'C2', 'N2', 'O2', 'CN', 'CO', 'NO', 'CO2', 'NCO', 'C+', 'O+', 'CO+', 'NO+', 'e-'};
    p = linspace(p2,1.5*rho1*u1*u1,6)/1e5; %bar
    allOneString = sprintf('%d,' , p);
    allOneString = allOneString(1:end-1);% strip final comma
    fileID = fopen('shock.inp','w');
    fprintf(fileID, '\n reac');
    for j=1:11
        fprintf(fileID,' name=%s  weight=%d \n',tgas{j},mf2(j));
    end
    fprintf(fileID, 'prob ions case=3  sp  s/r=%d   ',S2/8.314);
    fprintf(fileID, ' p,bar=%s  \n ',allOneString);
    fprintf(fileID, 'output massf trace=1.e-12   plot p h s \n ');
    fprintf(fileID, 'only Ar C N O C2 N2 O2 CN CO NO CO2 NCO C+ O+ CO+ NO+ e- \n end');
    fclose(fileID);
    
    [cmd,op] = system('cea'); 

    fid = fopen('shock.plt');
    content = fread(fid) ;
    fclose(fid);
    X = char(content.') ;
    Y = strrep(X, '#', '%') ;% replace string # with string %
    fid2 = fopen('shock.plt','wt') ;
    fwrite(fid2,Y) ;
    fclose (fid2) ;
    load shock.plt
    p02 = interp1(shock(:,2),shock(:,1),h0/1e3,'spline'); %interpolate h-p to get p02 corresponding to h0
    system('move shock.plt isentropic_h_p.plt');
    
    %% Get transport properties at stagnation (02) condition for use in Fay-Riddel

    fileID = fopen('shock.inp','w');
    fprintf(fileID, '\n reac');
    for j=1:11
        fprintf(fileID,' name=%s  weight=%d \n',tgas{j},mf2(j));
    end
    fprintf(fileID, 'prob ions case=3  sp  s/r=%d   ',S2/8.314);
    fprintf(fileID, ' p,bar=%d   \n',p02);
    fprintf(fileID, 'output tran massf trace=1.e-12   plot t rho s visc conductivity cp C N O C2 CN CO NO CO2 NCO C+ O+ CO+ NO+ e- \n ');
    fprintf(fileID, 'only Ar C N O C2 N2 O2 CN CO NO CO2 NCO C+ O+ CO+ NO+ e- \nend');
    fclose(fileID);
    [cmd,op] = system('cea'); 

    fid = fopen('shock.plt');
    content = fread(fid) ;
    fclose(fid);
    X = char(content.') ;
    Y = strrep(X, '#', '%') ;% replace string # with string %
    fid2 = fopen('shock.plt','wt') ;
    fwrite(fid2,Y) ;
    fclose (fid2) ;
    load shock.plt
    T02 = shock(:,1);
    rho_02 = shock(:,2);
    mu_02 = shock(:,4)*1e-4;%Ns/m^2
    k_02 = shock(:,5)*0.1;%W/m/K
    cp_02 = shock(:,6)*1e3;%J/kg/K
    mf_02 = shock(:,7:20); %mass frac of O N NO 
