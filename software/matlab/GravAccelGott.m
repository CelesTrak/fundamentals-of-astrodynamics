% ----------------------------------------------------------------------------
%
%                             function GravAccelGott
%
%   this function finds the acceleration for the gravity field using the normalized
%   Gottlieb approach. the arrays are indexed from 0 to coincide with the usual nomenclature
%   (eq 8-21 in my text). Fortran and MATLAB implementations will have indices of 1 greater
%   as they start at 1. note that this formulation is able to handle degree and order terms
%   larger then 170 due to the formulation. includes two-body contribution.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs        description                                   range / units
%    recef       position ecef vector of satellite             km
%    gravarr     gravity coefficients normalized
%    degree      degree of gravity field                       1..2159+
%    order       order of gravity field                        1..2159+
% 
% outputs
%    accel    eci frame acceeration (km/s^2
%    LegGottN  normalized alfs  (not need to pass back)
%
%  References
%    Eckman, Brown, Adamo 2016 NASA report
%
%  [accel] = GravAccelGott(recef, gravarr, degree, order, norm1, norm2, norm11, normn10, norm1m, norm2m, normn1 );
% ----------------------------------------------------------------------------

function [accel] = GravAccelGott(recef, gravarr, degree, order, norm1, norm2, norm11, normn10, norm1m, norm2m, normn1 )
    constastro;
    
    % calculate partial of acceleration wrt state
    partials = 'n';

    % this can be done ahead of time
    % these are 0 based arrays since they start at 2
    % for L = 2:degree+1 %RAE
    %     norm1(L) = sqrt((2*L+1.0) / (2*L-1.0)); % eq 3-1 RAE
    %     norm2(L) = sqrt((2*L+1.0) / (2*L-3.0)); % eq 3-2 RAE
    %     norm11(L) = sqrt((2*L+1.0) / (2*L))/(2*L-1.0); % eq 3-3 RAE
    %     normn10(L) = sqrt((L+1.0)*L * 0.5); % RAE
    % 
    %     for m = 1:L %RAE
    %         norm1m(L,m) = sqrt((L-m)*(2*L+1.0) / ((L+m)*(2*L-1.0))); % eq 3-4 RAE
    %         norm2m(L,m) = sqrt((L-m)*(L-m-1.0)*(2*L+1.0) / ((L+m)*(L+m-1.0)*(2*L-3.0))); % eq 3-5 RAE
    %         normn1(L,m) = sqrt((L+m+1.0)*(L-m)); % part of eq 3-9 RAE
    %     end %RAE
    % end %RAE
    
    %x = rnp*xin; %RAE
    r = sqrt(recef(1)^2 + recef(2)^2 + recef(3)^2);
    ri = 1.0/r;
    xor = recef(1)*ri;
    yor = recef(2)*ri;
    zor = recef(3)*ri;
    sinlat = zor;
    reor = re*ri;
    reorn = reor;
    muor2 = mu*ri*ri;
    
    LegGottN(1,1) = 1.0; %RAE
    LegGottN(1,2) = 0.0; %RAE
    LegGottN(1,3) = 0.0; %RAE
    LegGottN(2,2) = sqrt(3.0); %RAE    * coslat if comparing legpoly
    LegGottN(2,3) = 0.0; %RAE
    LegGottN(2,4) = 0.0; %RAE

    % --------------- sectoral 
    for L = 2:degree %RAE
        Li = L+1; %RAE
        LegGottN(Li,Li) = norm11(L)*LegGottN(L,L)*(2*L-1.0); %eq 3-15 * coslat if comparing legpoly
        LegGottN(Li,Li+1) = 0.0; %RAE
        LegGottN(Li,Li+2) = 0.0; %RAE
    end

    ctil(1) = 1.0; %RAE
    stil(1) = 0.0; %RAE
    ctil(2) = xor; %RAE
    stil(2) = yor; %RAE
    sumh = 0.0;
    sumgm = 1.0;  
    sumj = 0.0;
    sumk = 0.0;
    LegGottN(2,1) = sqrt(3)*sinlat; %RAE %norm 

    for L=2:degree
        Li = L + 1; %RAE
        reorn = reorn*reor;
        n2m1 = L + L - 1;
        nm1 = L - 1;
        np1 = L + 1;

        % --------------- tesserals(ni, m=ni-1) initial value
        LegGottN(Li,L) = normn1(L,L-1)*sinlat*LegGottN(Li,Li); % eq 3-16 RAE norm, but eq 3-9 is normn1 doesn't match eq 3-12???

        % --------------- zonals (ni, m=1)
        LegGottN(Li,1) = (n2m1*sinlat*norm1(L)*LegGottN(L,1) - nm1*norm2(L)*LegGottN(nm1,1))/L; % eq 3-17 RAE %norm
%fprintf(1,'%11.7f  %11.7f  %11.7f  %11.7f  ',norm1(n),LegGottN(n,1),norm2(n),LegGottN(nm1,1))

        % --------------- tesseral (n,m=2) initial value
        LegGottN(Li,2) = (n2m1*sinlat*norm1m(L,1)*LegGottN(L,2) - L*norm2m(L,1)*LegGottN(nm1,2))/(nm1); %RAE %norm

        sumhn = normn10(L)*LegGottN(Li,2)*gravarr.cNor(Li,1); %norm %RAE
        sumgmn = LegGottN(Li,1)*gravarr.cNor(Li,1)*np1; %RAE

        if partials == 'y'
            sumvn = pn(0) * gravarr.cNor(Li,1);
            summn = pn(2) * gravarr.cNor(Li,1)*upsn(O);
            sumpn = sumhn * np1;
            sumln = sumgamn * (np1 + 1.0);
        end

        if (order>0)
            for m = 2:L-2 % eckman wrong?, should be n-1?
                mi = m+1; %RAE
                % --------------- tesseral (n,m)
                LegGottN(Li,mi) = (n2m1*sinlat*norm1m(L,m)*LegGottN(L,mi) - ...
                    (nm1+m)*norm2m(L,m)*LegGottN(nm1,mi)) / (L-m); % eq 3-18 RAE %norm
            end  % for

            sumjn = 0.0;
            sumkn = 0.0;
            if partials == 'y'
                sumnn = 0.0;
                sumon = 0.0;
                sumqn = 0.0;
                sumrn = 0.0;
                sumsn = 0.0;
                sumtn = 0.0;    
            end
            ctil(Li) = ctil(2)*ctil(Li-1) - stil(2)*stil(Li-1); %RAE
            stil(Li) = stil(2)*ctil(Li-1) + ctil(2)*stil(Li-1); %RAE
            if L < order
                lim = L;
            else
                lim = order;
            end

            for m = 1:lim
                mi = m+1; %RAE
                mm1 = mi-1; %RAE
                mp1 = mi+1; %RAE
                mxpnm = m * LegGottN(Li,mi); %RAE
                bnmtil = gravarr.cNor(Li,mi) * ctil(mi) + gravarr.sNor(Li,mi) * stil(mi); %RAE
                sumhn = sumhn + normn1(L,m) * LegGottN(Li,mp1) * bnmtil; %RAE %norm
                sumgmn = sumgmn + (L+m+1) * LegGottN(Li,mi) * bnmtil; %RAE
                bnmtm1 = gravarr.cNor(Li,mi) * ctil(mm1) + gravarr.sNor(Li,mi) * stil(mm1); %RAE
                anmtm1 = gravarr.cNor(Li,mi) * stil(mm1) - gravarr.sNor(Li,mi) * ctil(mm1); %RAE
                sumjn = sumjn + mxpnm * bnmtm1;
                sumkn = sumkn - mxpnm * anmtm1;

                if partials == 'y'
                    npmpl = Li + mp1;
                    pnm = pn(mi);
                    pnmp1 = pn(mp1);
                    cnm = c(mi, 1);
                    snm = s(mi, 1);
                    ctil = ctil(mm1);
                    stil = stil(mm1);

                    bnmtml = cnm * ctil + snm * stil;
                    pnmbnm = pnm * bnmtll;
                    sumvn = sumvn + pnmbnm;
                    if m < L
                        z_pnmpl = zn(m) * pn(mpl);
                        sumhn = sumhn + z_pnmpl * bnmtil;
                        sumpn = sumpn + npmpl * z_pnmpl * bnmtil;
                        sumqn = sumqn + M * z_pnmp1 * bnmtml;
                        sumrn = sumrn- M * z_pnmpl * anmtml;
                    end
                    sumln = sumln + npmpl * (mpl + npl) *pnmbnm;
                    summn = summn + pn(mp2) * bnmtil*upsn(m);
                    sumsn = sumsn + npmpl * mxpnm * bnmtm1;
                    sumtn = sumtn - npmp1 * mxpnm * anmtml;
                    if m >= 2
                        mm2 = m - 2;
                        sumnn = sumnn + mml * mxpnm * (cnm * ctil(mm2) + snm*stil(mm2));
                        sumon = sumon + mml * mxpnm * (cnm * stil(mm2) - snm*ctil(mm2));
                    end
                end

            end  % for m

            sumj = sumj + reorn*sumjn;
            sumk = sumk + reorn*sumkn;
            if partials == 'y'
                sumn = sumn + reorn* sumnn;
                sumo = sumo + reorn* sumon;
                sumq = sumq + reorn* sumqn;
                sumr = sumr + reorn* sumrn;
                sums = sums + reorn* sumsn;
                sumt = sumt + reorn* sumtn; 
            end
        end  % if

        % these will have values when m == 0
        sumh = sumh + reorn*sumhn;
        sumgm = sumgm + reorn*sumgmn;
        if partials == 'y'
            sumv = sumv + reorn * sumvn;
            suml = suml + reorn * sumln;
            summ = summ + reorn * summn;
            sump = sump + reorn * sumpn;
        end
    end  % for

    lambda = sumgm + sinlat*sumh;
    % muor2 = -1.0;
    accel(1) = -muor2*(lambda*xor - sumj);
    accel(2) = -muor2*(lambda*yor - sumk);
    accel(3) = -muor2*(lambda*zor - sumh);
    %accel = rnp'*g; %RAE

    if partials == 'y'
        pot = muor * sumv;
        % Need to construct second partial matr1x_3x3
        gg = -(summ * sinlat + sump + sumh);
        ff = suml + lambda + sinlat * (sump + sumh - gg);
        d1 = sinlat * sumq + sums;
        d2 = sinlat * sumr + sumt;
        dgdx (1, 1) = muor3 * ((ff* xovr- 2.0 * d1) * xovr - lambda + sumn);
        dgdx (2, 2) = muor3 * ((ff * yovr - 2.0 * d2) * yovr - lambda - sumn);
        dgdx (3, 3) = muor3 * ((ff * zovr + 2.0 * gg) * zovr - lambda + summ);
        temp = muor3 * ((ff * yovr - d2) * xovr - d1 * yovr - sumo);
        dgdx (1, 2) = temp;
        dgdx (2, 1) = temp;
        temp = muor3 * ((ff* xovr - d1) * zovr + gg * xovr + sumq);
        dgdx (1, 3) = temp;
        dgdx (3, 1) = temp;
        temp = muor3 * ((ff * yovr - d2) * zovr + gg * yovr + sumr);
        dgdx (2, 3) = temp;
        dgdx (3, 2) = temp;
    end

end  % GravAccelGott
