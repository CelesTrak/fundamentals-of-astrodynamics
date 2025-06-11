% ----------------------------------------------------------------------------
%
%                           function GravAccelMont
%
%   this function finds the acceleration for the gravity field. the acceleration is
%   found in the body fixed frame. rotation back to inertial is done after this
%   routine. this is the montenbruck approach amnd it uses unnormalzied coefficients.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs        description                                   range / units
%    recef       - position vector ECEF                          km
%    gravarr     - gravitational coefficients
%    normArr     - normalization coefficients
%    degree      - degree of gravity field                        1..85
%    order       - order of gravity field                         1..85
%
%  outputs       :
%    apert       - ecef perturbation acceleration                 km / s^2
%
%  locals :
%    L, m        - degree and order indices
%    VArr        - array of trig terms
%    WArr        - array of trig terms
%
%  coupling      :
%   none
%
%  references :
%    vallado       2022, 600
%
%  [aPert] = GravAccelMont ( recef, gravarr, normArr, degree, order);
% ----------------------------------------------------------------------------

function [aPert] = GravAccelMont ( recef, gravarr, normArr, degree, order)
    constastro;

    % find normalized values - can be done ahead of time
    [normArr] = gravnorm(degree + 1);

    % --------------------find latgc and lon----------------------
    [latgc, latgd, lon, hellp] = ecef2ll(recef);

    aPert(1) = 0.0;
    aPert(2) = 0.0;
    aPert(3) = 0.0;

    % Body-fixed position
    % Auxiliary quantities
    magr2 = mag(recef);              
    magr2 = magr2 * magr2;

    rho = re * re / magr2;
    x0 = re * recef(1) / magr2;           
    y0 = re * recef(2) / magr2;     
    z0 = re * recef(3) / magr2;

    % Calculate zonal terms V(n,0); set W(n,0)=0.0
    V(0+1, 0+1) = re / sqrt(magr2);
    W(0+1, 0+1) = 0.0;
    V(1+1, 0+1) = z0 * V(0+1, 0+1);
    W(1+1, 0+1) = 0.0;
    for L = 2: degree + 1
        Li = L + 1;
        %conv = unnormArr(L, m);
        V(Li, 0+1) = ((2 * L - 1) * z0 * V(Li - 1, 0+1) - (L - 1) * rho * V(Li - 2, 0+1)) / L;
        W(Li, 0+1) = 0.0;
    end

    % Calculate tesseral and sectorial terms
    for m = 1:order + 1
        mi = m + 1;
        % note V and W not formulated for normalized
        % use un-normalized for now
        V(mi, mi) = (2 * m - 1) * (x0 * V(mi - 1, mi - 1) - y0 * W(mi - 1, mi - 1));
        W(mi, mi) = (2 * m - 1) * (x0 * W(mi - 1, mi - 1) + y0 * V(mi - 1, mi - 1));
        if m <= degree
            V(mi + 1, mi) = (2 * m + 1) * z0 * V(mi, mi);
            W(mi + 1, mi) = (2 * m + 1) * z0 * W(mi, mi);
        end
        for L = m + 2: degree + 1
            Li = L + 1;
            V(Li, mi) = ((2 * L - 1) * z0 * V(Li - 1, mi) - (L + m - 1) * rho * V(Li - 2, mi)) / (L - m);
            W(Li, mi) = ((2 * L - 1) * z0 * W(Li - 1, mi) - (L + m - 1) * rho * W(Li - 2, mi)) / (L - m);
        end
    end

    % Calculate accelerations, note the order is switched here
    for m = 0:order
        mi = m + 1;
        for L = m:degree
            Li = L + 1;
            if (m == 0)
                C = gravarr.cNor(Li, 1)*normArr(Li, 1);   % = C_n,0
                aPert(1) = aPert(1) - C * V(Li + 1, 2);
                aPert(2) = aPert(2) - C * W(Li + 1, 2);
                aPert(3) = aPert(3) - (L + 1) * C * V(Li + 1, 0+1);
            else
                C = gravarr.cNor(Li, mi)*normArr(Li, mi);   % = C_n,m
                S = gravarr.sNor(Li, mi)*normArr(Li, mi);   % = S_n,m
                Fac = 0.5 * (L - m + 1) * (L - m + 2);

                aPert(1) = aPert(1) + 0.5  * (-C * V(Li + 1, mi + 1) - S * W(Li + 1, mi + 1))...
                    + Fac * (+C * V(Li + 1, mi - 1) + S * W(Li + 1, mi - 1));
                aPert(2) = aPert(2) + 0.5  * (-C * W(Li + 1, mi + 1) + S * V(Li + 1, mi + 1))...
                    + Fac * (-C * W(Li + 1, mi - 1) + S * V(Li + 1, mi - 1));
                aPert(3) = aPert(3) + (L - m + 1) * (-C * V(Li + 1, mi) - S * W(Li + 1, mi));
            end
        end
    end

    % Body-fixed acceleration
    temp = mu / (re * re);
    aPert(1) = temp * aPert(1);
    aPert(2) = temp * aPert(2);
    aPert(3) = temp * aPert(3);

end  % GravAccelMont
