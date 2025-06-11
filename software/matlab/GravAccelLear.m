% ----------------------------------------------------------------------------
%
%                             function GravAccelLear
%
%   this function finds the acceleration for the gravity field using the normalized
%   Lear approach. the arrays are indexed from 0 to coincide with the usual nomenclature
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
%    accel       eci frame acceeration                         km/s^2
%    pnm, ppnm   normalized alfs   (not need to pass back)
%
%  References
%    Eckman, Brown, Adamo 2016 NASA report
%
%  [accel] = GravAccelLear(recef, gravarr, degree, order, norm1, norm2, norm11, norm1m, norm2m);
% ----------------------------------------------------------------------------

function [accel] = GravAccelLear(recef, gravarr, degree, order, norm1, norm2, norm11, norm1m, norm2m)
    constastro;

    % this can be done ahead of time
    % these are 0 based arrays since they start at 2
    % for L  =  2:degree %RAE
    %     norm1(L) = sqrt((2.0*L + 1.0) / (2.0*L - 1.0)); %RAE
    %     norm2(L) = sqrt((2.0*L + 1.0) / (2.0*L - 3.0)); %RAE
    %     norm11(L) = sqrt((2.0*L + 1.0) / (2.0*L)) / (2.0*L - 1.0); %RAE
    %     for m  =  1:L %RAE
    %         norm1m(L,m) = sqrt((L - m)*(2.0*L + 1.0) / ((L + m)*(2.0*L - 1.0))); %RAE
    %         norm2m(L,m) = sqrt((L - m)*(L - m - 1.0)*(2.0*L + 1.0) / ((L + m)*(L + m - 1.0)*(2.0*L - 3.0))); %RAE
    %     end %RAE
    % end %RAE

    for L = 2:degree %RAE
        pnm(L - 1,L) = 0;
    end

    %recef = rnp*recef; %RAE
    e1 = recef(1)^2 + recef(2)^2;
    magr2 = e1 + recef(3)^2;
    magr = sqrt(magr2);
    r1 = sqrt(e1);
    sphi = recef(3) / magr;
    cphi = r1 / magr;
    sm(1) = 0.0;
    cm(1) = 1.0;
    if (r1 ~= 0.0)
        sm(1) = recef(2) / r1;
        cm(1) = recef(1) / r1;
    end
    reor(1) = re / magr;
    reor(2) = reor(1)^2;
    sm(2) = 2.0*cm(1)*sm(1);
    cm(2) = 2.0*cm(1)^2 - 1.0;
    root3 = sqrt(3.0); %RAE
    root5 = sqrt(5.0); %RAE
    pn(1) = root3*sphi; %norm
    pn(2) = root5*(3.0*sphi^2 - 1.0) * 0.5; %norm

    ppn(1) = root3; %norm
    ppn(2) = root5*3*sphi; %norm
    pnm(1,1) = root3; %norm
    pnm(2,2) = root5*root3*cphi * 0.5; %norm
    pnm(2,1) = root5*root3*sphi; %norm %RAE
    ppnm(1,1) =  -root3*sphi; %norm
    ppnm(2,2) =  -root3*root5*sphi*cphi; %norm
    ppnm(2,1) = root5*root3*(1.0 - 2.0*sphi^2); %norm

    if (degree >= 3) %RAE
        for L = 3:degree %RAE
            nm1 = L - 1;
            nm2 = L - 2;
            reor(L) = reor(nm1)*reor(1);
            sm(L) = 2.0*cm(1)*sm(nm1) - sm(nm2);
            cm(L) = 2.0*cm(1)*cm(nm1) - cm(nm2);
            e1 = 2.0*L - 1.0;
            pn(L) = (e1*sphi*norm1(L)*pn(nm1) - nm1*norm2(L)*pn(nm2)) / L; %norm
            ppn(L) = norm1(L)*(sphi*ppn(nm1) + L*pn(nm1)); %norm
            pnm(L,L) = e1*cphi*norm11(L)*pnm(nm1,nm1); %norm
            ppnm(L,L) =  -L*sphi*pnm(L,L);
        end
        for L = 3:degree %RAE
            nm1 = L - 1;
            e1 = (2.0*L - 1.0)*sphi;
            e2 =  - L*sphi;
            for m = 1:nm1
                e3 = norm1m(L,m)*pnm(nm1,m); %norm
                e4 = L + m;
                e5 = (e1*e3 - (e4 - 1.0)*norm2m(L,m)*pnm(L - 2,m)) / (L - m); %norm
                pnm(L,m) = e5;
                ppnm(L,m) = e2*e5 + e4*e3;
            end
        end
    end
    asph(1) = -1.0;
    asph(3) = 0.0;
    for L = 2:degree %RAE
        ni = L + 1; %RAE
        e1 = gravarr.cNor(ni,1)*reor(L); %RAE
        asph(1) = asph(1) - (L + 1)*e1*pn(L);
        asph(3) = asph(3) + e1*ppn(L);
    end
    asph(3) = cphi*asph(3);
    t1 = 0.0;
    t3 = 0.0;
    asph(2) = 0.0;
    for L = 2:degree %RAE
        ni = L + 1; %RAE
        e1 = 0.0;
        e2 = 0.0;
        e3 = 0.0;
        if (L>order) %RAE
            nmodel = order; %RAE
        else %RAE
            nmodel = L; %RAE
        end %RAE

        for m = 1:nmodel %RAE
            mi = m + 1; %RAE
            tsnm = gravarr.sNor(ni,mi); %RAE
            tcnm = gravarr.cNor(ni,mi); %RAE
            tsm = sm(m);
            tcm = cm(m);
            tpnm = pnm(L,m);
            e4 = tsnm*tsm + tcnm*tcm;
            e1 = e1 + e4*tpnm;
            e2 = e2 + m*(tsnm*tcm - tcnm*tsm)*tpnm;
            e3 = e3 + e4*ppnm(L,m);
        end
        t1 = t1 + (L + 1.0)*reor(L)*e1;
        asph(2) = asph(2) + reor(L)*e2;
        t3 = t3 + reor(L)*e3;
    end

    e4 = mu / magr2;
    asph(1) = e4*(asph(1) - cphi*t1);
    asph(2) = e4*asph(2);
    asph(3) = e4*(asph(3) + t3);

    e5 = asph(1)*cphi - asph(3)*sphi;
    accel(1,1) = e5*cm(1) - asph(2)*sm(1); %RAE
    accel(2,1) = e5*sm(1) + asph(2)*cm(1); %RAE
    accel(3,1) = asph(1)*sphi + asph(3)*cphi; %RAE
    %accel = rnp'*agr; %RAE
    return

end  %  GravAccelLear