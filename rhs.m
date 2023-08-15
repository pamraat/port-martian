%% Defines the RHS of governing equation
% See notation.pdf for details
function dy = rhs(x,y)
% This functions takes in y = [ rho T \gamma_i] and
% gives output 

    global R0 mdotstardash rhostardash Tstardash m0 rstar theta Theta l coefVec speciesVec viscflag;
    global alpha beta Kf Kb z;
    Af = Kf(:,1);
    Bf = Kf(:,2);
    Tdf = Kf(:,3);
    Astar = pi*rstar*rstar;        % Throat area (m^2)
    if strcmp(viscflag,'TRUE')
        A = @(x) pi*((tan(theta*pi/180)*x*l + rstar).^2)/Astar;     
        dA = @(x) 2*pi*tan(theta*pi/180)*l*(tan(theta*pi/180)*x*l + rstar)/Astar; 
    else
        A = @(x) pi*((tan(Theta*pi/180)*x*l + rstar).^2)/Astar;     
        dA = @(x) 2*pi*tan(Theta*pi/180)*l*(tan(Theta*pi/180)*x*l + rstar)/Astar; 
    end
    rhodash = y(1)*rhostardash;
    Tdash = y(2)*Tstardash;
    udash = mdotstardash/rhodash/A(x)/Astar; 
    u = udash/sqrt(R0*Tstardash/m0); 
    
%     gamma = y(3:8);
%     [Cp, h] = getCph_T_mod(Tdash, R0, coefVec, speciesVec);
%     h = h*y(2); 
%     H = m0*sum(gamma.*h');
%     u = sqrt(2*(H0-H));

    Kc=eqcons(Tdash,coefVec);
    Fr = sum(alpha')'-1;
    Br = sum(beta')'-1;
    Kfr = Tdash.^Bf.*Af.*exp(-Tdf./Tdash).*(1e6).^(-Fr);
%     Kbr = Tdash.^Bb.*Ab.*exp(-Tdb./Tdash).*(1e6).^(-Br);
    for i=1:29
        Kbr(i) = Kfr(i)/Kc(i);
    end

    Gamma = zeros(24,1); Gamma(1:16) = y(4:19);
    for i=1:8
        Gamma(i+16) = sum(z(8,:).*y(4:19)');
    end

    Rfr = zeros(29,1); Rbr = zeros(29,1);
    for i=1:29
        Rfr(i) = Kfr(i)*prod((Gamma*rhodash).^(alpha(i,:)'));
        Rbr(i) = Kbr(i)*prod((Gamma*rhodash).^(beta(i,:)'));
    end

    dy = zeros(19,1);
    % These are the gammas.
    dy(1) = u*u/m0/A(x)*dA(x);
    dy(2) = dy(1)/y(2);

    for i=1:16
        dy(i+3) = l/rhodash/udash*sum((beta(:,i) - alpha(:,i)).*(Rfr-Rbr));
    end
    dy(3) = 0; % Ar  assumed inert
