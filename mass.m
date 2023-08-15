%% Function outputs the state dependent mass matrix for the ODE solver.
% See notation.pdf for details
function mm = mass(x, y)
    global R0 mdotstardash rhostardash Tstardash m0 rstar Theta theta l coefVec speciesVec viscflag...
    ;
    Astar = pi*rstar*rstar;        % Throat area (m^2)
    if strcmp(viscflag,'TRUE')
        A = @(x) pi*((tan(theta*pi/180)*x*l + rstar).^2)/Astar; 
    else
        A = @(x) pi*((tan(Theta*pi/180)*x*l + rstar).^2)/Astar; 
    end
    rhodash = y(1)*rhostardash;
    Tdash = y(2)*Tstardash;
    udash = mdotstardash/rhodash/A(x)/Astar;  
    u = udash/sqrt(R0*Tstardash/m0); 
    [Cp, hn, h_0] = getCph_T_mod(Tdash, R0, coefVec, speciesVec);
    h = hn*y(2);          
    
%     H = m0*sum(gamma.*h');
%     u = sqrt(2*(H0-H));
    m = 1/sum(y(3:2+length(speciesVec)));

    mm = zeros(2+length(speciesVec),2+length(speciesVec));
    mm(1,1) = -u*u/m0/y(1);
    mm(1,2) = sum(y(3:2+length(speciesVec)).*Cp');
    mm(1,3:2+length(speciesVec)) = h;
    mm(2,1) = (1/m - u*u/m0/y(2))/y(1);
    mm(2,2) = 1/m/y(2);
    mm(2,3:2+length(speciesVec)) = 1;
    mm(3:2+length(speciesVec),3:2+length(speciesVec)) = eye(length(speciesVec),length(speciesVec));

