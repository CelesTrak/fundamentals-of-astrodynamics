% ----------------------------------------------------------------------------
%
%                           function GravAccelGTDS
%
%   this function finds the acceleration for the gravity field.
%
%  author        : david vallado                    719-573-2600  10 oct 2019
%
%  inputs        description                                   range / units
%    recef       - position vector ECEF                          km
%    order       - size of gravity field                         1..360
%    normal      - normalized in file                            'y', 'n'
%    gravData.c  - gravitational coefficients (in gravData)
%    gravData.s  - gravitational coefficients (in gravData)
%
%  outputs       :
%    apert       - efc perturbation acceleration                  km / s^2
%
%  locals :
%    conv        - conversion to normalize
%    L, m        - degree and order indices
%    trigArr     - array of trigonometric terms
%    LegArr      - array of Legendre polynomials
%    VArr        - array of trig terms
%    WArr        - array of trig terms
%
%  coupling      :
%   LegPoly      - find the unnormalized Legendre polynomials through recursion
%   TrigPoly     - find the trigonmetric terms through recursion
%
%  references :
%    vallado       2013, 597, Eq 8-57
%
%  [aPert, aPert1] = fullgeop ( recef, order, gravarr);
% ----------------------------------------------------------------------------

function [aPert, aPert1] = GravAccelGTDS ( recef, order, gravarr)
    constastro;

    aPert = zeros(order+3,order+3);
    aPert1 = zeros(order+3,order+3);

    sumM1 = 0.0;
    sumM2 = 0.0;
    sumM3 = 0.0;

    % --------------------find latgc and lon---------------------- }
    [latgc, latgd, lon, hellp] = ecef2ll(recef);

    % ---------------------Find Legendre polynomials -------------- }
    [legarrMU, legarrGU, legarrMN, legarrGN] = legpolyn(latgc, order+3);

    % compare values
    %double[,] diff = new double[6, 6];
    %diff = legarrMU - legarrGU;

    [trigarr, VArr, WArr] = trigpoly(recef, latgc, lon, order+3);

    magr = mag(recef);

    % now do the partials
    %     temp0 = (recef(1)^2 + recef(2)^2);
    %     temp1 = 1.0 / magr;
    %     temp2 = recef(3) / (magr^2*sqrt(recef(1)^2 + recef(2)^2));
    %     temp3 = 1.0 / temp0;
    %     temp4 = sqrt(temp0)/magr^2;
    %
    %     agrav(1) = (temp1* partrr - temp2*partrlat) * recef(1) - (temp3*partrlon)*recef(2);
    %     agrav(2) = (temp1* partrr - temp2*partrlat) * recef(2) - (temp3*partrlon)*recef(1);
    %     agrav(3) = (temp1* partrr) * recef(3) + temp4*partrlon;

    oor = re / magr;
    distPartR = 0.0;
    distPartPhi = 0.0;
    distPartLon = 0.0;
    temp = oor;

    for L = 2: order
        Li = L + 1;

        % will do the power as each L is indexed
        temp = temp * oor;
        sumM1 = 0.0;
        sumM2 = 0.0;
        sumM3 = 0.0;

        for m = 0:L
            mi = m + 1;
            sumM1 = sumM1 + legarrGU(Li, mi) * (gravarr.cNor(Li, mi) * trigarr(mi, 1+1) + gravarr.sNor(Li, mi) * trigarr(mi, 0+1));
            sumM2 = sumM2 + (legarrGU(Li, mi ++ 1) - trigarr(mi, 2+1) * legarrGU(Li, mi)) * ...
                (gravarr.cNor(Li, mi) * trigarr(mi, 1+1) + gravarr.sNor(Li, mi) * trigarr(mi, 0+1));
            sumM3 = sumM3 + m * legarrGU(Li, mi) * (gravarr.sNor(Li, mi) * trigarr(mi, 1+1) - gravarr.cNor(Li, mi) * trigarr(mi, 0+1));
        end % for m

        distPartR = distPartR + temp * (Li) * sumM1;
        distPartPhi = distPartPhi + temp * sumM2;
        distPartLon = distPartLon + temp * sumM3;
    end % for L

    distPartR = -oor * oor * sumM1;
    distPartPhi = oor * sumM2;
    distPartLon = oor * sumM3;

    fprintf(1,'----------Non - spherical perturbative acceleration ------------ \n');
    RDelta = sqrt(recef(1) * recef(1) + recef(2) * recef(2));
    oordelta = 1.0 / RDelta;
    temp = oor * distPartR - recef(3) * oor * oor * oordelta * distPartPhi;
    tmp = mum / (magr^3);

    aPert(1) = temp * recef(1) - oordelta * distPartLon * recef(1) - tmp * recef(1);
    aPert(2) = temp * recef(2) + oordelta * distPartLon * recef(1) - tmp * recef(2);
    aPert(3) = oor * distPartR * recef(3) + oor * oor * RDelta * distPartPhi - tmp * recef(3);


    % --------------- montenbruck approach
    magrecef = mag(recef);
    for L = 2: order
        Li = L + 1;
        % will do the power as each L is indexed }
        temp = mu / (magrecef * magrecef);

        for m = 0:L
            mi = m + 1;
            temp1 = factorial(L - m + 2) / factorial(L - m);

            % unnormalized
            % if (m == 0)
            %     aPert1(1) = aPert1(1) + temp * (-gravData.c(Li, mi) * VArr(L + 1, 1));
            %     aPert1(2) = aPert1(2) + temp * (-gravData.c(Li, mi) * WArr(L + 1, 1));
            % else
            %     aPert1(1) = aPert1(1) + 0.5 * temp * ((-gravData.c(Li, mi) * VArr(L + 1, m + 1) - gravData.s(Li, mi) * WArr(L + 1, m + 1)) +...
            %         temp1 * (gravData.c(Li, mi) * VArr(L + 1, m - 1) + gravData.s(L, m) * WArr(L + 1, m - 1)));
            %     aPert1(2) = aPert1(2) + 0.5 * temp * ((-gravData.c(Li, mi) * WArr(L + 1, m + 1) + gravData.s(Li, mi) * VArr(L + 1, m + 1)) +...
            %         temp1 * (-gravData.c(Li, mi) * WArr(L + 1, m - 1) + gravData.s(Li, mi) * VArr(L + 1, m - 1)));
            % end
            % aPert1(3) = aPert1(3) + temp * ((L - m + 1) * (-gravData.c(Li, mi) * VArr(L + 1, m) - gravData.s(Li, mi) * WArr(L + 1, m)));

            % normalized
            if (m == 0)
                aPert1(1) = aPert1(1) + temp * (-gravarr.cNor(Li, mi) * VArr(Li, 1));
                aPert1(2) = aPert1(2) + temp * (-gravarr.cNor(Li, mi) * WArr(Li, 1));
            else
                aPert1(1) = aPert1(1) + 0.5 * temp * ((-gravarr.cNor(Li, mi) * VArr(Li, mi) - gravarr.sNor(Li, mi) * WArr(Li, mi)) + ...
                    temp1 * (gravarr.cNor(L, m) * VArr(Li, mi - 1) + gravarr.sNor(Li, mi) * WArr(Li, mi - 1)));
                aPert1(2) = aPert1(2) + 0.5 * temp * ((-gravarr.cNor(Li, mi) * WArr(Li, mi) + gravarr.sNor(Li, mi) * VArr(Li, mi)) + ...
                    temp1 * (-gravarr.cNor(Li, mi) * WArr(Li, mi - 1) + gravarr.sNor(Li, mi) * VArr(Li, mi - 1)));
            end
            aPert1(3) = aPert1(3) + temp * ((L - m + 1) * (-gravarr.cNor(Li, mi) * VArr(Li, mi) - gravarr.sNor(Li, mi) * WArr(Li, mi)));
        end % for m
    end

    %             int n;                           % Loop counters
    %             double r_sqr, rho, Fac;               % Auxiliary quantities
    %             double x0, y0, z0;                      % Normalized coordinates
    %             double ax, ay, az;                      % Acceleration vector
    %             double C, S;                           % Gravitational coefficients
    %             double(,) V = new double(order+2, order + 2);
    %             double(,) W = new double(order+2, order + 2);

    % now try montenbruck code approach
    % Body-fixed position
    % r_bf = E * r;
    % Auxiliary quantities
    r_sqr = dot(recef, recef);               % Square of distance
    rho = re * re / r_sqr;
    x0 = re * recef(1) / r_sqr;          % Normalized
    y0 = re * recef(2) / r_sqr;          % coordinates
    z0 = re * recef(3) / r_sqr;

    % Evaluate harmonic functions
    %   V_nm = (R_ref/r)^(n+1) * P_nm(sin(phi)) * cos(m*lambda)
    % and
    %   W_nm = (R_ref/r)^(n+1) * P_nm(sin(phi)) * sin(m*lambda)
    % up to degree and order n_max+1
    % Calculate zonal terms V(n,0); set W(n,0)=0.0

    V(0+1, 0+1) = re / sqrt(r_sqr);
    W(0+1, 0+1) = 0.0;
    V(1+1, 0+1) = z0 * V(0+1, 0+1);
    W(1+1, 0+1) = 0.0;
    for n = 2: order
        ni = n + 1;
        V(ni, 0+1) = ((2 * n - 1) * z0 * V(ni - 1, 0+1) - (n - 1) * rho * V(ni - 2, 0+1)) / n;
        W(ni, 0+1) = 0.0;
    end

    % Calculate tesseral and sectorial terms
    for m = 1:order+1
        mi = m + 1;
        % Calculate V(m,m) .. V(n_max+1,m)
        V(mi, mi) = (2 * m - 1) * (x0 * V(mi - 1, mi - 1) - y0 * W(mi - 1, mi - 1));
        W(mi, mi) = (2 * m - 1) * (x0 * W(mi - 1, mi - 1) + y0 * V(mi - 1, mi - 1));
        if (m <= order)
            V(mi, mi) = (2 * mi) * z0 * V(mi, mi);
            W(mi, mi) = (2 * mi) * z0 * W(mi, mi);
        end
        for n = m + 2: order + 1
            ni = n + 1;
            V(ni, mi) = ((2 * n - 1) * z0 * V(ni - 1, mi) - (n + m - 1) * rho * V(ni - 2, mi)) / (n - m);
            W(ni, mi) = ((2 * n - 1) * z0 * W(ni - 1, mi) - (n + m - 1) * rho * W(ni - 2, mi)) / (n - m);
        end
    end

    % Calculate accelerations ax,ay,az
    for m = 0:order
        mi = m + 1;
        for (n = m:order)
            ni = n + 1;
            if (m == 0)
                C = gravarr.cNor(ni, 1);   % = C_n,0
                aPert1(1) = aPert1(1) - C * V(ni, 1);
                aPert1(2) = aPert1(2) - C * W(ni, 1);
                aPert1(3) = aPert1(3) - (n + 1) * C * V(ni, 0+1);
            else
                C = gravarr.cNor(ni, mi);   % = C_n,m
                S = gravarr.sNor(mi - 1, ni); % = S_n,m
                Fac = 0.5 * (n - m + 1) * (n - m + 2);
                aPert1(1) = aPert1(1) + 0.5 * (-C * V(ni, mi) - S * W(ni, mi))...
                    + Fac * (+C * V(ni, mi - 1) + S * W(ni, mi - 1));
                aPert1(2) = aPert1(2) + 0.5 * (-C * W(ni, mi) + S * V(ni, mi))...
                    + Fac * (-C * W(ni, mi - 1) + S * V(ni, mi - 1));
                aPert1(3) = aPert1(3)+(n - m + 1) * (-C * V(ni, mi) - S * W(ni, mi));
            end
        end
    end

    % Body-fixed acceleration
    aPert2(1) = (mu / (re * re)) * aPert1(1);
    aPert2(2) = (mu / (re * re)) * aPert1(2);
    aPert2(3) = (mu / (re * re)) * aPert1(3);

    % Inertial acceleration
    %return Transp(E) * a_bf;

end  % GravAccelGTDS
