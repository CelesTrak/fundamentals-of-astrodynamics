% ----------------------------------------------------------------------------
%
%                             function GravAccelPines
%
%   this function finds the acceleration for the gravity field using the normalized
%   Pines approach. the arrays are indexed from 0 to coincide with the usual nomenclature
%   (eq 8-21 in my text). Fortran and MATLAB implementations will have indices of 1 greater
%   as they start at 1. note that this formulation is able to handle degree and order terms
%   larger then 170 due to the formulation. 
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs        description                                   range  /  units
%    recef       position ecef vector of satellite             km
%    gravarr     gravity coefficients normalized
%    degree      degree of gravity field                       1..2159+
%    order       order of gravity field                        1..2159+
% 
% outputs
%    accel       ecef frame acceeration (km / s^2
%    LegGottN    normalized alfs   (not need to pass back)
%
%  References
%    Eckman, Brown, Adamo 2016 NASA report
%
%  [accel] = GravAccelPines(recef, gravarr, degree, order);
% ----------------------------------------------------------------------------

function [accel] = GravAccelPines(recef, gravarr, degree, order)
    constastro;

    %recef = rnp*recef; %RAE
    magr = norm(recef);
    S = recef(1) / magr;
    T = recef(2) / magr;
    U = recef(3) / magr;
    LegPinesN = zeros(degree+3, degree+3); %RAE
    %ANMx = zeros(NMAX+3, NMAX+3); %dav
    
    LegPinesN(1, 1) = sqrt(2.0); %norm
    for m = 0:degree+2 %RAE
        mi = m + 1;
        if (m ~= 0) % DIAGONAL RECURSION
            LegPinesN(mi, mi) = sqrt(1.0+(1.0 / (2.0*m))) * LegPinesN(mi-1, mi-1); %norm
        end
        if (m ~= degree+2) % FIRST OFF-DIAGONAL RECURSION %RAE
            LegPinesN(mi+1, mi) = sqrt(2.0*m+3.0) * U*LegPinesN(mi, mi); %norm
        end
        if (m < degree+1) % COLUMN RECURSION %RAE
            for L = m+2:degree+2 %RAE
                Li = L + 1;
                ALPHA_NUM = (2.0*L + 1.0)*(2.0*L - 1.0);
                ALPHA_DEN = (L - m)*(L + m);
                ALPHA = sqrt(ALPHA_NUM / ALPHA_DEN);
                BETA_NUM = (2.0*L + 1.0)*(L - m - 1.0)*(L + m - 1.0);
                BETA_DEN = (2.0*L - 3.0)*(L + m)*(L - m);
                BETA = sqrt(BETA_NUM / BETA_DEN);
                LegPinesN(Li, mi) = ALPHA*U*LegPinesN(Li - 1, mi) - BETA*LegPinesN(Li - 2, mi); %norm
            end
        end
    end
    for L = 0:degree+2 %RAE
        Li = L + 1;
        LegPinesN(Li, 1) = LegPinesN(Li, 1)*sqrt(0.5); %norm
    end
    
    %ANMx = ANM;  % dav
    
    RM = zeros(order+2, 1); %RAE
    IM = zeros(order+2, 1); %RAE
    RM(1) = 0.0;
    IM(1) = 0.0;
    RM(2) = 1.0; %RAE
    IM(2) = 0.0; %RAE
    for m = 1:order %RAE
        mp2 = m + 2; %RAE
        RM(mp2) = S*RM(mp2 - 1) - T*IM(mp2 - 1);
        IM(mp2) = S*IM(mp2 - 1) + T*RM(mp2 - 1);
    end
    RHO = mu / (re*magr);
    reor = re / (magr);
    G1 = 0.0;
    G2 = 0.0;
    G3 = 0.0;
    G4 = 0.0;
    for L = 0:degree
        Li = L + 1;
        G1TEMP = 0.0;
        G2TEMP = 0.0;
        G3TEMP = 0.0;
        G4TEMP = 0.0;
        SM = 0.5;
        if (L>order) %RAE
            nmodel = order; %RAE
        else %RAE
            nmodel = L; %RAE
        end %RAE
        
        for m = 0:nmodel %RAE
            mi = m + 1;
            mp2 = m + 2; %RAE
            DNM = gravarr.cNor(Li, mi)*RM(mp2) + gravarr.sNor(Li, mi)*IM(mp2);
            ENM = gravarr.cNor(Li, mi)*RM(mp2-1) + gravarr.sNor(Li, mi)*IM(mp2-1);
            FNM = gravarr.sNor(Li, mi)*RM(mp2-1) - gravarr.cNor(Li, mi)*IM(mp2-1);
            ALPHA = sqrt(SM*(L - m)*(L + m + 1)); %norm
            G1TEMP = G1TEMP + LegPinesN(Li, mi)*(m)*ENM;
            G2TEMP = G2TEMP + LegPinesN(Li, mi)*(m)*FNM;
            G3TEMP = G3TEMP + ALPHA*LegPinesN(Li, mi + 1)*DNM; %norm
            G4TEMP = G4TEMP + ((L+m+1)*LegPinesN(Li, mi) + ALPHA*U*LegPinesN(Li, mi + 1))*DNM;  %norm
            %ANMx(Li, mi) = (L+m+1)*ANMx(Li, mi);  % dav
            %ANMx(Li, mi+1) = ALPHA*U*ANMx(Li, mi+1);   % dav
            if(m == 0)
                SM = 1.0;
            end %norm
        end

        RHO = reor*RHO;
        G1 = G1 + RHO*G1TEMP;
        G2 = G2 + RHO*G2TEMP;
        G3 = G3 + RHO*G3TEMP;
        G4 = G4 + RHO*G4TEMP;
    end

    accel(1) = G1 - G4 * S;
    accel(2) = G2 - G4 * T;
    accel(3) = G3 - G4 * U;
    %accel=rnp'*G_F; %RAE

    return
    

