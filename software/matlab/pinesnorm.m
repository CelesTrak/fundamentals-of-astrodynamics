% inputs
%    mu       gravitaional paramater
%    req      equatorial radius
%    x        position cvector of satellite km
%    cnm, snm gravity coefficients normalized
%    nnax     degree
%    mmax     order
%    rnp      3x3 identity matrix
% 
% outputs
%    accel    acceeration (km/s^2
%    ANM      normalized alfs
%
function [ANM, accel] = pinesnorm(MU, REQ, recef, CNM, SNM, NMAX, MMAX)
    %recef = rnp*recef; %RAE
    RMAG = norm(recef);
    S = recef(1)/RMAG;
    T = recef(2)/RMAG;
    U = recef(3)/RMAG;
    ANM = zeros(NMAX+3, NMAX+3); %RAE
    %ANMx = zeros(NMAX+3, NMAX+3); %dav
    
    ANM(1, 1) = sqrt(2.0); %norm
    for M = 0:NMAX+2 %RAE
        mi = M + 1;
        if (M ~= 0) % DIAGONAL RECURSION
            ANM(mi, mi) = sqrt(1+(1/(2*M))) * ANM(mi-1, mi-1); %norm
        end
        if (M ~= NMAX+2) % FIRST OFF-DIAGONAL RECURSION %RAE
            ANM(mi+1, mi) = sqrt(2*M+3) * U*ANM(mi, mi); %norm
        end
        if (M < NMAX+1) % COLUMN RECURSION %RAE
            for N = M+2:NMAX+2 %RAE
                Li = N + 1;
                ALPHA_NUM = (2*N+1)*(2*N-1);
                ALPHA_DEN = (N-M)*(N+M);
                ALPHA = sqrt(ALPHA_NUM/ALPHA_DEN);
                BETA_NUM = (2*N+1)*(N-M-1)*(N+M-1);
                BETA_DEN = (2*N-3)*(N+M)*(N-M);
                BETA = sqrt(BETA_NUM/BETA_DEN);
                ANM(Li, mi) = ALPHA*U*ANM(Li-1, mi) - BETA*ANM(Li-2, mi); %norm
            end
        end
    end
    for N = 0:NMAX+2 %RAE
        Li = N + 1;
        ANM(Li, 1) = ANM(Li, 1)*sqrt(0.5); %norm
    end
    
    %ANMx = ANM;  % dav
    
    RM = zeros(MMAX+2, 1); %RAE
    IM = zeros(MMAX+2, 1); %RAE
    RM(1) = 0.0D0;
    IM(1) = 0.0D0;
    RM(2) = 1.0D0; %RAE
    IM(2) = 0.0D0; %RAE
    for M = 1:MMAX %RAE
        M_RI = M + 2; %RAE
        RM(M_RI) = S*RM(M_RI-1) - T*IM(M_RI-1);
        IM(M_RI) = S*IM(M_RI-1) + T*RM(M_RI-1);
    end
    RHO = (MU)/(REQ*RMAG);
    RHOP = (REQ)/(RMAG);
    G1 = 0.0D0;
    G2 = 0.0D0;
    G3 = 0.0D0;
    G4 = 0.0D0;
    for N = 0:NMAX
        Li = N + 1;
        G1TEMP = 0.0D0;
        G2TEMP = 0.0D0;
        G3TEMP = 0.0D0;
        G4TEMP = 0.0D0;
        SM = 0.5;
        if (N>MMAX) %RAE
            nmodel=MMAX; %RAE
        else %RAE
            nmodel=N; %RAE
        end %RAE
        for M = 0:nmodel %RAE
            mi = M + 1;
            M_RI = M + 2; %RAE
            DNM = CNM(Li, mi)*RM(M_RI) + SNM(Li, mi)*IM(M_RI);
            ENM = CNM(Li, mi)*RM(M_RI-1) + SNM(Li, mi)*IM(M_RI-1);
            FNM = SNM(Li, mi)*RM(M_RI-1) - CNM(Li, mi)*IM(M_RI-1);
            ALPHA = sqrt(SM*(N-M)*(N+M+1)); %norm
            G1TEMP = G1TEMP + ANM(Li, mi)*(M)*ENM;
            G2TEMP = G2TEMP + ANM(Li, mi)*(M)*FNM;
            G3TEMP = G3TEMP + ALPHA*ANM(Li, mi+1)*DNM; %norm
            G4TEMP = G4TEMP + ((N+M+1)*ANM(Li, mi) + ALPHA*U*ANM(Li, mi+1))*DNM;  %norm
            %ANMx(Li, mi) = (N+M+1)*ANMx(Li, mi);  % dav
            %ANMx(Li, mi+1) = ALPHA*U*ANMx(Li, mi+1);   % dav
            if(M == 0)
                SM = 1.0;
            end %norm
        end
        RHO = RHOP*RHO;
        G1 = G1 + RHO*G1TEMP;
        G2 = G2 + RHO*G2TEMP;
        G3 = G3 + RHO*G3TEMP;
        G4 = G4 + RHO*G4TEMP;
    end
    accel = [G1 - G4*S;G2 - G4*T;G3 - G4*U];
    %accel=rnp'*G_F; %RAE
    
    return
    

