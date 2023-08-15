function Kc=eqcons(T,coefVec)
Ru=8.314;
% HoverRT = zeros(1,size(coefVec,3));
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
%             HoverRT(i) = a1 + a2*T/2 + a3*(T^2)/3 + a4*(T^3)/4 + a5*(T^4)/5 + a6/T; %h[j/mole]/(Ro*T)
        H(i) = Ru*T*(a1 + a2*T/2 + a3*(T^2)/3 + a4*(T^3)/4 + a5*(T^4)/5 + a6/T);
        S(i)= Ru*(a1*log(T)+a2*T+a3*T^2/2+a4*T^3/3+a5*T^4/4+a7);
        G(i)=H(i)-T*S(i);
        H_0(i) = a6*Ru;
    end
    
    delG(1)=2*G(2)-G(5);
    delG(2)=2*G(3)-G(6);
    delG(3)=2*G(4)-G(7);
    delG(4)=G(2)+G(3)-G(8);
    delG(5)=G(2)+G(4)-G(9);
    delG(6)=G(3)+G(4)-G(10);
    delG(7)=G(9)+G(4)-G(11);
    delG(8)=G(3)+G(9)-G(12);
    delG(9)=G(7)+G(3)-G(4)-G(10);
    delG(10)=G(3)+G(10)-G(6)-G(4);
    delG(11)=G(7)+G(2)-G(9)-G(4);
    delG(12)=G(5)+G(4)-G(9)-G(2);
    delG(13)=G(4)+G(8)-G(9)-G(3);
    delG(14)=G(8)+G(3)-G(6)-G(2);
    delG(15)=G(10)+G(2)-G(8)-G(4);
    delG(16)=G(5)+G(3)-G(8)-G(2);
    delG(17)=G(7)+G(9)-G(11)-G(4);
    delG(18)=G(4)+G(12)-G(8)-G(7);
    delG(19)=G(9)+G(12)-G(11)-G(8);
    delG(20)=G(3)+G(12)-G(8)-G(10);
    delG(21)=G(12)+G(4)-G(10)-G(9);
    delG(22)=G(2)+G(12)-G(9)-G(8);
    delG(23)=G(16)+G(17)-G(4)-G(3);
    delG(24)=G(15)+G(17)-G(4)-G(2);
    delG(25)=G(13)+G(10)-G(2)-G(16);
    delG(26)=G(14)+G(6)-G(3)-G(16);
    delG(27)=G(2)+G(15)-G(13)-G(9);
    delG(28)=G(13)+2*G(17)-G(13);
    delG(29)=G(14)+2*G(17)-G(4);
        
    Kc=exp(-delG/(Ru*T));
% now compute mixture properties
%     Cp = CpoverR;
%     h = H/(Ru*T);
%     h_0 = H_0; %h@ T=0K