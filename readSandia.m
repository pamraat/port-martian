%% Retrieves the NASA polynomial coefficents for the Cp vs T curve fit
% Data is retrived from the file NASA.DAT 
% there are two fits for data: a high and a low temperature fit
% Each fit has 7 coefficients; we need to reassemble these. 
% Note: we are not reading in the date, and also there is a value in
% the 5th position in line 4 which is the enthalpy of formation divided
% by the gas constant R; this can be computed using the equations from
% the coefficients already read into the file.
% for more information, see the following resources:
%    http://garfield.chem.elte.hu/burcat/THERM.DAT
%    http://www.galcit.caltech.edu/EDL/public/formats/chemkin.html
%    http://www.galcit.caltech.edu/EDL/public/formats/nasaold.html
%    http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19780009781_1978009781.pdf
% From aerosfrosh code by Mattheww Campbell, Stanford Univerity
function [coefVec, nameEle, numEle] = readSandia(speciesVec, nasaFile)
% This function reads the thermodynamic data file



% make a checklist vector to ensure we got everything
    checkList = zeros(1,length(speciesVec));
% open the file
    fileID = fopen(nasaFile);
% read the first line
    firstLine = fgets(fileID);
% is this the right starter file?
    tempIndex = strfind(firstLine,'THERMO'); % index where 'THERMO' begins
    if isempty(tempIndex) % if 'THERMO' does not occur in first line
        error('%1s is empty',nasaFile)
    end
% read the next line
    secondLine = fgets(fileID);
% find the general temperature limits
    genTlow = str2num(secondLine(1:10)); % general low temp limit (K)
    genTcom = str2num(secondLine(11:20)); % general common temp (K)
    genThigh = str2num(secondLine(21:30)); % general high temp limit (K)
% read the file
    closeNow = false;
    while (closeNow ~= true)
        % reset variables
            elements = {''};
            numAtoms = [];
            Tlow = 0;
            Tcom = 0;
            Thigh = 0;
            TlowWarn = false;
            TcomWarn = false;
            ThighWarn = false;
        % get the first line of the species
            L1 = fgets(fileID);
        % check for ending
            endText = upper(strtrim(L1(1:3)));
            if strcmp(endText,'END')
                closeNow = true;
        % have we found everything?
            elseif sum(checkList) == length(checkList) % all 1's so this sums to the length if full! 
                closeNow = true; 
        % if not at the end, keep reading the file
            else
                % grab the species name
                    speciesName = upper(strtrim(L1(1:18))); % remove white space, make UPPERCASE
                % check the name; sometimes there is other text beside it
                    refIndex = strfind(speciesName,'REF'); % find where 'REF' occurs
                    if ~isempty(refIndex) % if 'REF' occured, then... 
                        speciesName = strtrim(speciesName(1:refIndex-1)); % take a substring
                    end
                % grab the elements and number of atoms of each
                    elpos = [25, 30, 35, 40, 74];
                    for i=1:5
                        % generate indicies
                            tmpIdx1 = elpos(i)+0; %element char 1
                            tmpIdx2 = elpos(i)+1; %element char 2
                            tmpIdx3 = elpos(i)+2; %number of elem digit 1
                            tmpIdx4 = elpos(i)+4; %number of elem digit 3
                        % grab the element name and the number of atoms
                            tmpEl = upper(strtrim(L1(tmpIdx1:tmpIdx2))); % UPPERCASE
                            tmpNum = str2num(L1(tmpIdx3:tmpIdx4));
                        % make sure the element is there 
                        % not all compounds have 5 elements
                        % if it's there, store it.  if not, throw it away. 
                        check = strcmp(tmpEl,''); % true if nothing, false if something
                        % some files have a zero in column 74
                        % this throws off the algorithm so check for it 
                        check2 = strcmp(tmpEl,'0'); 
                        if (check == false) && (check2 == false)
                            elements{1,i} = tmpEl;
                            numAtoms(1,i) = tmpNum;
                        end
                    end
                % grab the phase
                    phase = upper(L1(45)); % UPPERCASE letter
                % grab the low, middle, and high temperatures (K)
                    Tlow = str2num(L1(46:55));
                    Thigh = str2num(L1(56:65));
                    Tcom = str2num(L1(66:73));
                % check the inputs
                    if sum(size(Tlow))>2 || sum(size(Tcom))>2 || sum(size(Thigh))>2
                        error('SANDIA.DAT file format issue for %7s; Check for tabs instead of spaces.',speciesName)
                    end
                % if we did not see Tlow etc, use the general value
                    if isempty(Tlow); 
                        Tlow = genTlow;
                        TlowWarn = true; % flag that temperature was replaced
                    end
                    if isempty(Tcom);
                        Tcom = genTcom;
                        TcomWarn = true; % flag that temperature was replaced
                    end
                    if isempty(Thigh);
                        Thigh = genThigh;
                        ThighWarn = true; % flag that temperature was replaced
                    end
                % get line two
                    L2 = fgets(fileID);
                % grab the temperature coefficients here
                    Ha1 = str2num(L2(1:15)); % high temperature coefficient a1
                    Ha2 = str2num(L2(16:30));
                    Ha3 = str2num(L2(31:45));
                    Ha4 = str2num(L2(46:60));
                    Ha5 = str2num(L2(61:75));
                % get line three
                    L3 = fgets(fileID);
                % grab the temperature coefficients here
                    Ha6 = str2num(L3(1:15));
                    Ha7 = str2num(L3(16:30));
                    La1 = str2num(L3(31:45)); % low temperature coefficient a1
                    La2 = str2num(L3(46:60));
                    La3 = str2num(L3(61:75));            
                % get line four
                    L4 = fgets(fileID);
                % grab the temperature coefficients here
                    La4 = str2num(L4(1:15));
                    La5 = str2num(L4(16:30));
                    La6 = str2num(L4(31:45));
                    La7 = str2num(L4(46:60));              
                % if we want this species, store it, otherwise disregard
                    for i=1:length(speciesVec);
                        % if the same name as one of our wanted species, 
                        % and if it is a gas
                        % and if we have not grabbed it yet
                        if strcmp(speciesName,speciesVec{i}) && strcmp(phase,'G') && checkList(i)==0; 
                            % coefficeint and temperature matrix
                            coefVec(1,:,i) = [Ha1 Ha2 Ha3 Ha4 Ha5 Ha6 Ha7];
                            coefVec(2,:,i) = [La1 La2 La3 La4 La5 La6 La7];
                            coefVec(3,:,i) = [Tlow Tcom Thigh 0 0 0 0];
                            % element name and number vectors
                            for j=1:length(elements)
                                nameEle(i,j) = elements(1,j);  % use parentheses not curly brackets
                                numEle(i,j) = numAtoms(1,j);
                            end
                            % display warnings about changing temperatures
                            if TlowWarn
                                fprintf('WARNING: Defaulting Tlow to %4.2f K for %7s\n',genTlow,speciesName);
                            end
                            if TcomWarn
                                fprintf('WARNING: Defaulting Tcom to %4.2f K for %7s\n',genTcom,speciesName);
                            end
                            if ThighWarn
                                fprintf('WARNING: Defaulting Thigh to %4.2f K for %7s\n',genThigh,speciesName);
                            end
                            % check off that we got it
                            checkList(i) = 1; 
                        end
                    end
            end
    end
% close the file
    fclose(fileID);
% did we get everything?
    if sum(checkList) ~= length(checkList)
        for i=1:length(checkList)
            if checkList(i) == 0
                error('Species not in %5s: %5s',nasaFile,speciesVec{i})
            end
        end
    end
    
end