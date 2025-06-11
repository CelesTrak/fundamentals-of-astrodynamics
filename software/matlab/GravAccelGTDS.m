% ----------------------------------------------------------------------------
%
%                           function GravAccelGTDS
%
%   this function finds the acceleration for the gravity field. the acceleration is
%   found in the body fixed frame. rotation back to inertial is done after this
%   routine. this is the gtds approach and it uses unnormalized values.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs        description                                   range / units
%    recef       - position vector ECEF                          km
%    gravarr     - gravitational coefficients
%    normArr     - normalization coefficients
%    degree      - degree of gravity field                       1..85
%    order       - order of gravity field                        1..85
%
%  outputs       :
%    apert       - ecef perturbation acceleration                km / s^2
%
%  locals :
%    L, m        - degree and order indices
%    trigArr     - array of trigonometric terms
%    LegArr      - array of Legendre polynomials
%
%  coupling      :
%   LegPoly      - find the unnormalized Legendre polynomials through recursion
%   TrigPoly     - find the trigonmetric terms through recursion
%
%  references :
%    vallado       2022, 600
%
%  [aPert] = GravAccelGTDS ( recef, gravarr, normArr, degree, order);
% ----------------------------------------------------------------------------

function [aPert] = GravAccelGTDS ( recef, gravarr, normArr, degree, order)
    constastro;

    sumM1 = 0.0;
    sumM2 = 0.0;
    sumM3 = 0.0;

    % --------------------find latgc and lon----------------------
    [latgc, latgd, lon, hellp] = ecef2ll(recef);

    % ---------------------Find Legendre polynomials --------------
    [legarrGU] = legpolyGTDS (latgc, degree+2, order+2);

    [trigArr] = trigpolyGTDS(latgc, lon, degree+2);

    magr = mag(recef);
    oor = 1.0 / magr;
    reor = re * oor;
    dRdr = 0.0;
    dRdlat = 0.0;
    dRdlon = 0.0;

    for L = 2: degree
        Li = L + 1;

        % will do the power as each L is found
        temp = reor * reor;
        sumM1 = 0.0;
        sumM2 = 0.0;
        sumM3 = 0.0;

        for m = 0 : L
            mi = m + 1;
            % take normalized coefficients and revert to unnormalized
            temparg = gravarr.cNor(Li, mi)*normArr(Li, mi) * trigArr(mi, 1+1) ...
                + gravarr.sNor(Li, mi)*normArr(Li, mi) * trigArr(mi, 0+1);
            sumM1 = sumM1 + legarrGU(Li, mi) * temparg;
            sumM2 = sumM2 + (legarrGU(Li, mi + 1) - trigArr(mi, 2+1) * legarrGU(Li, mi)) * temparg;  % yes m+1 in 1st legp
            sumM3 = sumM3 + m * legarrGU(Li, mi) * (gravarr.sNor(Li, mi)*normArr(Li, mi) * trigArr(mi, 1+1) ...
                - gravarr.cNor(Li, mi)*normArr(Li, mi) * trigArr(mi, 0+1));
        end % for m

        dRdr = dRdr + temp * (L + 1) * sumM1;  
        dRdlat = dRdlat + temp * sumM2;  
        dRdlon = dRdlon + temp * sumM3;  
    end % for L

    muor = mu * oor;
    dRdr = -muor * oor * dRdr;
    dRdlat = muor * dRdlat;
    dRdlon = muor * dRdlon;

    %fprintf(1,'----------Non - spherical perturbative acceleration ------------ \n');
    RDelta = sqrt(recef(1) * recef(1) + recef(2) * recef(2));
    oordeltasqrt = 1.0 / RDelta;
    temp1 = oor * dRdr - recef(3) * oor * oor * oordeltasqrt * dRdlat;

    aPert(1) = temp1 * recef(1) - oordeltasqrt * oordeltasqrt * dRdlon * recef(2);
    aPert(2) = temp1 * recef(2) + oordeltasqrt * oordeltasqrt * dRdlon * recef(1);
    aPert(3) = oor * dRdr * recef(3) + oor * oor * RDelta * dRdlat;

    % add in two-body
    % two-body term
    % tmp(1) = mu / (magr^3);
    % tmp(2) = mu / (magr^3);
    % tmp(3) = mu / (magr^3);
    % aPert(1) = aPert(1) + tmp(1)*recef(1);
    % aPert(2) = aPert(2) + tmp(2)*recef(2);
    % aPert(3) = aPert(3) + tmp(3)*recef(3);

end  % GravAccelGTDS
