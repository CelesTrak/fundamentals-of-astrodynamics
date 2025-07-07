% --------------------------------------------------------------------------
%                           function lambhodograph
%
% this function accomplishes 180 deg transfer(and 360 deg) for lambert problem.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    r1    - ijk position vector 1                            km
%    r2    - ijk position vector 2                            km
%    v1    - intiial ijk velocity vector 1                    km/s
%    p     - semiparamater of transfer orbit                  km
%    ecc   - eccentricity of transfer orbit                   km
%    dnu   - true anomaly delta for transfer orbit            rad
%    dtsec - time between r1 and r2                           s
%    dnu - true anomaly change                                rad
%
%  outputs       :
%    v1t - ijk transfer velocity vector                       km/s
%    v2t - ijk transfer velocity vector                       km/s
%
%  references :
%    Thompson JGCD 2013 v34 n6 1925
%    Thompson AAS GNC 2018
%
% [v1t, v2t] = lambhodograph( r1, v1, r2, p, a, ecc, dnu, dtsec );
% ------------------------------------------------------------------------------

function [v1t, v2t] = lambhodograph( r1, v1, r2, p, ecc, dnu, dtsec )
    constastro;
    twopi = 2.0 * pi;

    magr1 = mag(r1);
    magr2 = mag(r2);

    oop = 1.0 / p;

    magr1 = mag(r1);
    oomagr1 = 1.0 / magr1;
    magr2 = mag(r2);
    oomagr2 = 1.0 / magr2;
    eps = 0.001 / magr2;  % 1e-14

    a = mu * (1.0 * oomagr1 - 1.0 * oop);  % not the semi - major axis
    b = (mu * ecc * oop)^2 - a * a;

    if b <= 0.0
        x1 = 0.0;
    else
        x1 = -sqrt(b);
    end
    % 180 deg, and multiple 180 deg transfers
    if abs(sin(dnu)) < eps
        rcrv = cross(r1, v1);
        nvec = rcrv / mag(rcrv);
        if ecc < 1.0
            ptx = twopi * sqrt(p^3/( mu*(1.0-ecc^2)^3));
            if (mod(dtsec,ptx) > ptx*0.5)
                x1 = -x1;
            end
        end
        fprintf(1,'less than\n');
    else
        % more common path?
        y2a = mu*oop - x1*sin(dnu) + a*cos(dnu);
        y2b = mu*oop + x1*sin(dnu) + a*cos(dnu);
        if abs(mu/magr2 - y2b) < abs(mu/magr2 - y2a)
            x1 = -x1;
        end
        % depending on the cross product, this will be normal or in plane,
        % or could even be a fan
        rcrr = cross(r1, r2);
        nvec = rcrr / mag(rcrr); % if this is r1, v1, the transfer is coplanar!
        if (mod(dnu, twopi) > pi)
            nvec = -nvec;
        end
        %       fprintf(1,'gtr than\n');
    end

    v1t = (sqrt(mu*p) * oomagr1) * ((x1/mu)*r1 + cross(nvec,r1)*oomagr1 );
    x2  = x1*cos(dnu) + a*sin(dnu);
    v2t = (sqrt(mu*p) * oomagr2) * ((x2/mu)*r2 + cross(nvec,r2)*oomagr2 );

end

