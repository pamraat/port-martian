%% Prepares the run of NASA CEA  and formats the output data from CEA for further use. Resuts from this are used for
% 1. Calculating shock tube end wall conditions, includig isentropic
% relaxation to measured P5
% 2. Computing the equilibrium nozzle flow
% 3. Initialisation of the Non equilibrium simulation
% Note: the original CEA code has been modified to enable execution without
%  a filename input.
% Ref: https://www.grc.nasa.gov/www/CEAWeb/
function [p5,T5,h5,rho5,gamma5,mf5,rhostardash, Tstardash, ustardash,massfrac,pstardash] = cea_run(u1,p1,T1,tgas,mf1, ar_sup,p5e)
    
    global speciesvec;
    %% Create .inp file for cea shock tube run, to find nozzle inlet conditions
    fileID = fopen('shock.inp','w');
    fprintf(fileID, '\n reac');
    for j=1:length(speciesvec)
        fprintf(fileID,' name=%s  weight=%d t(k) %d\n',tgas{j},mf1(j),T1);
    end
    fprintf(fileID, 'prob case=3  p,bar=%d shock incd eq ref eq ions u1=%d  \n output massf trace=1.e-12 \n only Ar C N O C2 N2 O2 CN CO NO CO2 NCO C+ O+ CO+ NO+ e-\n end',p1,u1);
    fclose(fileID);

    %% Run 'cea' with shock.inp as input file
    [cmd,op] = system('cea'); 
    
    %% Read shock tube conditions from cea output 
             content = fileread( 'shock.out' ) ;
    strexp = regexp(content,'U5.+\*C[+].+\*CN', 'once'); %check if C+ is present among the products

    if isempty(strexp)
        a  = regexp(content,...
        'U5\,.+P\,[^0-9]+([^a-zA-Z]+).+T\,[^0-9]+([^a-zA-Z]+).+RHO\,[^0-9]+([^a-zA-Z]+).+H\,[^0-9]+([^a-zA-Z]+).+GAMMAs[^0-9]+([\S]+).+\*e-[^0-9]+([\S]+).+\*Ar[^0-9]+([\S]+).+\*C[^0-9]+([\S]+).+\*CN[^0-9]+([\S]+).+\*CO[^0-9]+([\S]+).+\*CO+[^0-9]+([\S]+).+\*CO2[^0-9]+([\S]+).+\*C2[^0-9]+([\S]+).+\*N[^0-9]+([\S]+).+\NCO[^0-9]+([\S]+).+\*NO[^0-9]+([\S]+).+\*NO+[^0-9]+([\S]+).+\*N2[^0-9]+([\S]+).+\*O[^0-9]+([\S]+).+\*O+[^0-9]+([\S]+).+\*O2[^0-9]+([\S]+)'...
        ,'tokens');%looking for string 'U5,', skips all char till 'P,', skips all non-digit, token () is all non alphabet char ..
        a = a{:};
        a1 = {a{1:3} '0' a{4:21}};
    else
        a  = regexp(content,...
        'U5\,.+P\,[^0-9]+([^a-zA-Z]+).+T\,[^0-9]+([^a-zA-Z]+).+RHO\,[^0-9]+([^a-zA-Z]+).+H\,[^0-9]+([^a-zA-Z]+).+GAMMAs[^0-9]+([\S]+).+\*e-[^0-9]+([\S]+).+\*Ar[^0-9]+([\S]+).+\*C[^0-9]+([\S]+).+\*C+[^0-9]+([\S]+).+\*CN[^0-9]+([\S]+).+\*CO[^0-9]+([\S]+).+\*CO+[^0-9]+([\S]+).+\*CO2[^0-9]+([\S]+).+\*C2[^0-9]+([\S]+).+\*N[^0-9]+([\S]+).+\NCO[^0-9]+([\S]+).+\*NO[^0-9]+([\S]+).+\*NO+[^0-9]+([\S]+).+\*N2[^0-9]+([\S]+).+\*O[^0-9]+([\S]+).+\*O+[^0-9]+([\S]+).+\*O2[^0-9]+([\S]+)'...
        ,'tokens');%looking for string 'U5,', skips all char till 'P,', skips all non-digit, token () is all non alphabet char ..
        a1 = a{:};
    end
     
    for i=1:size(a1,2)
        a2(i) = cellstr(regexprep(char(a1(i)),'([0-9])-','$1e-')); %replacing '-' by 'e-'
        a2(i) = cellstr(regexprep(char(a2(i)),'([0-9])\s([0-9])','$1e+$2')); %replacing '\s ' by 'e+'
    end
    %A=zeros(size(a2,1),size(a2,2));
    A = str2double(a2); %thermodynamic properties
    [p5, T5, rho5, h5, cp_cv_5, mf5] = deal(A(1,1), A(1,2), A(1,3) ,A(1,4), A(1,5), A(1,6:22)); %shock tube stagnation conditions (bar,K,kg/m3,kJ/kgK)
    mf5temp = [ mf5(2) mf5(3) mf5(10) mf5(15) mf5(9) mf5(14) mf5(17) mf5(5) mf5(6) mf5(12) mf5(8)...
        mf5(11) mf5(4) mf5(16) mf5(7) mf5(13) mf5(1)];
    clearvars mf5;
    mf5 = mf5temp;
    clearvars mf5temp;
    [cmd,op] = system('move shock.out shock_tube.out');
    
    %% Isentropic relaxation to measured pressure p5e
    if p5>p5e
        p5e = p5/p5e;
        %  Create .inp file for cea eq-nozzle relaxation run
        %Fuel = {'Ar',	'C', 'N', 'O', 'C2', 'N2', 'O2', 'CN', 'CO', 'NO', 'CO2', 'NCO', 'C+', 'O+', 'CO+', 'NO+', 'e-'};
        fileID = fopen('shock.inp','w');
        fprintf(fileID, 'prob ions case=4 ro equilibrium\n o/f 0.000\n p,bar %d \n  pip %d \n only Ar C N O C2 N2 O2 CN CO NO CO2 NCO C+ O+ CO+ NO+ e-\n',p5,p5e);
        fprintf(fileID, '\n reac\n');

 
        for j=1:12
            fprintf(fileID,'fuel %s  wt%%=%d t,k=%d\n',speciesvec{j},mf1(j),T5);
        end
        fprintf(fileID, 'oxid  Air wt%%=100 t,k=298.15 \n output massf trace=1.e-12 \n plot\n p t mach son Ar C N O C2 N2 O2 CN CO NO CO2 NCO h rho gam\n end');
%         CEA can print only 20 data sets in one go, but we have 17
%         species, and 7 conditions, total 24 data sets
%         so we will re run the study to extract 5 data sets later, we will
%         extract 19 first.
        fclose(fileID);
        % Run cea nozzle flow calculation
        [cmd,op] = system('cea');
        [cmd,op] = system('move shock.out shocktube_relax.out');
        % relaxed nozzle inlet conditions
        fid = fopen('shock.plt');
        content = fread(fid);
        fclose(fid);
        X = char(content.');
        Y = strrep(X, '#', '%') ;% replace string # with string %
        fid2 = fopen('shock.plt','wt') ;
        fwrite(fid2,Y);
        fclose (fid2);
        load shock.plt
        p5 = shock(3,1);
        T5 = shock(3,2);
        ustardash = shock(3,3)*shock(3,4); 
        mf5(1:12) = shock(3,5:16);
        h5 = shock(3,17);
        rho5 = shock(3,18);
        gamma5 = shock(3,19); % The 19 data sets extracted
        
        
        
        fileID = fopen('shock.inp','w');
        fprintf(fileID, 'prob ions case=4 ro equilibrium\n o/f 0.0001\n p,bar %d \n  pip %d \nonly Ar C N O C2 N2 O2 CN CO NO CO2 NCO C+ O+ CO+ NO+ e-\n',p5,p5e);
        fprintf(fileID, '\n reac\n');   
        speciesvec1=speciesvec(1:12);
        for j=1:12
            fprintf(fileID,'fuel %s  wt%%=%d t,k=%d\n',speciesvec1{j},mf1(j),T5);
        end
        fprintf(fileID, 'oxid  Air wt%%=100 t,k=298.15 \n output massf trace=1.e-12\n plot\n C+ O+ CO+ NO+ e- \n end');
        fclose(fileID);
        % Run cea nozzle flow calculation
        [cmd,op] = system('cea');
        [cmd,op] = system('move shock.out shocktube_relax.out');
        % relaxed nozzle inlet conditions
        fid = fopen('shock.plt');
        content = fread(fid);
        fclose(fid);
        X = char(content.');
        Y = strrep(X, '#', '%') ;% replace string # with string %
        fid2 = fopen('shock.plt','wt') ;
        fwrite(fid2,Y);
        fclose (fid2);
        load shock.plt
        mf5(13:17) = shock(3,:);
        [cmd,op] = system('move shock.plt shocktube_relax.plt');        
    end
    %%  Create .inp file for cea eq-nozzle  run. Generates supersonic inlet condition for the Non Eq run or perform complete Eq nozzle simulation
    %Fuel = {'Ar',	'C', 'N', 'O', 'C2', 'N2', 'O2', 'CN', 'CO', 'NO', 'CO2', 'NCO', 'C+' 'O+' 'CO+' NO+' 'e-'};
    fileID = fopen('shock.inp','w');
    fprintf(fileID, 'prob ions case=4 ro equilibrium\n o/f 0.0001\n p,bar %d \n  supar',p5);
    fprintf(fileID, ' %f',ar_sup);
    fprintf(fileID, '\n reac\n');
    for j=1:12
        fprintf(fileID,'fuel %s  wt%%=%d t,k=%d\n',speciesvec{j},mf1(j),T5);
    end
    fprintf(fileID, 'oxid  Air wt%%=100 t,k=298.15 \n output massf trace=1.e-12\n plot\n  rho t mach son Ar C N O C2 N2 O2 CN CO NO CO2 NCO p ae\n end');
    fclose(fileID);
    % Run cea nozzle flow calculation
    %[cmd,op] = system('move shock.out shock_tube.out');
    [cmd,op] = system('cea');
    [cmd,op] = system('move shock.out nozzle_throat.out');
    % Exit conditions
    fid = fopen('shock.plt');
    content = fread(fid) ;
    fclose(fid);
    X = char(content.') ;
    Y = strrep(X, '#', '%') ;% replace string '#' with string '%'
    fid2 = fopen('shock.plt','wt') ;
    fwrite(fid2,Y) ;
    fclose (fid2) ;
    load shock.plt
    rhostardash = shock(size(shock,1),1);
    Tstardash = shock(size(shock,1),2);
    ustardash = shock(size(shock,1),3)*shock(size(shock,1),4); 
    massfrac = shock(size(shock,1),5:16);
    pstardash =  shock(size(shock,1),17)*1e5;%Pa
    system('move shock.plt eq_nozzle.plt'); %file with output data
    
    
    fileID = fopen('shock.inp','w');
    fprintf(fileID, 'prob ions case=4 ro equilibrium\n o/f 0.0001\n p,bar %d \n  supar',p5);
    fprintf(fileID, ' %f',ar_sup);
    fprintf(fileID, '\n reac\n');
    for j=1:12
        fprintf(fileID,'fuel %s  wt%%=%d t,k=%d\n',speciesvec{j},mf1(j),T5);
    end
    fprintf(fileID, 'oxid  Air wt%%=100 t,k=298.15 \n output massf trace=1.e-12\n plot\n C+ O+ CO+ NO+ e- \n end');
    fclose(fileID);
    % Run cea nozzle flow calculation
    [cmd,op] = system('cea');
    [cmd,op] = system('move shock.out nozzle_throat.out');
    % Exit conditions
    fid = fopen('shock.plt');
    content = fread(fid) ;
    fclose(fid);
    X = char(content.') ;
    Y = strrep(X, '#', '%') ;% replace string '#' with string '%'
    fid2 = fopen('shock.plt','wt') ;
    fwrite(fid2,Y) ;
    fclose (fid2) ;
    load shock.plt
    massfrac(13:17) = shock(size(shock,1),:);
    system('move shock.plt eq_nozzle1.plt'); %file with output data
