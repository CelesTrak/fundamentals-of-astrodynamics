%  ------------------------------------------------------------------------------
%
%                           procedure lambertTmaxrp
%
%  this procedure solves lambert's problem and finds the TOF for maximum rp
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    r1          - ijk position vector 1                    km
%    r2          - ijk position vector 2                    km
%    dm          - direction of motion                      'L', 'S'
%    de          - orbital energy                           'L', 'H'
%    nrev        - number of revs to complete               0, 1, 2, 3,
%
%  outputs       :
%    tmin        - minimum time of flight                   sec
%    tminp       - minimum parabolic tof                    sec
%    tminenergy  - minimum energy tof                       sec
%
%  locals        :
%    i           - index
%    loops       -
%    cosdeltanu  -
%    sindeltanu  -
%    dnu         -
%    chord       -
%    s           -
%
%  coupling      :
%    mag         - magnitude of a vector
%    dot         - dot product
%
%  references    :
%    thompson       2019
%
%  [tmaxrp, v1t] = lamberttmaxrp( r1, r2, dm, nrev );
%  ----------------------------------------------------------------------------

function [tmaxrp, v1t] = lamberttmaxrp( r1, r2, dm, nrev )
    mu = 3.986004415e5;

    v1t(1) = 0.0;
    v1t(2) = 0.0;
    v1t(3) = 0.0;

    magr1 = mag(r1);   % km
    magr2 = mag(r2);
    cosdeltanu = dot(r1, r2) / (magr1 * magr2);
    % make sure it's not more than 1.0
    if (abs(cosdeltanu) > 1.0)
        cosdeltanu = 1.0 * sgn(cosdeltanu);
    end
    rcrossr = cross(r1, r2);
    if (dm == 'S')
        sindeltanu = mag(rcrossr) / (magr1 * magr2);
    else
        sindeltanu = -mag(rcrossr) / (magr1 * magr2);
    end
    % these are the same
    chord = sqrt(magr1 * magr1 + magr2 * magr2 - 2.0 * magr1 * magr2 * cosdeltanu);
    %chord = MathTimeLibr.mag(r2 - r1);

    %s = (magr1 + magr2 + chord) * 0.5;

    y1 = mu / magr1;
    y2 = mu / magr2;

    % ------------- nearly circular endpoints
    tempvec = r1 + r2;
    % circular orbit  to within 1 meter
    if (mag(tempvec) < 0.001)
        c = sqrt(y1);
        x1 = 0.0;
        x2 = 0.0;
        r = 0.0;
        tmaxrp = (2.0 * nrev * pi + atan2(sindeltanu, cosdeltanu)) * (mu / (c^3));
    else
        %
        if (magr1 < magr2)
            c = sqrt((y2 - y1 * cosdeltanu) / (1.0 - cosdeltanu));
            x1 = 0.0;
            x2 = (y1 - c * c) * sindeltanu;
        else
            c = sqrt((y1 - y2 * cosdeltanu) / (1.0 - cosdeltanu));
            x1 = (-y2 + c * c) * sindeltanu;
            x2 = 0.0;
        end
        r = (sqrt(x1 * x1 + (y1 - c * c)^2)) / c;

        % check if acos is larger than 1
        temp = c * (r * r + y1 - c * c) / (r * y1);
        if (abs(temp) > 1.0)
            e1 = sign(temp) * acos(1.0);
        else
            e1 = acos(temp);
        end
        if (x1 < 0.0)
            e1 = 2.0 * pi - e1;
        end
        temp = c * (r * r + y2 - c * c) / (r * y2);
        if (abs(temp) > 1.0)
            e2 = sign(temp) * acos(1.0);
        else
            e2 = acos(temp);
        end
        if (x2 < 0.0)
            e2 = 2.0 * pi - e2;
        end
        if (e2 < e1)
            e2 = e2 + 2.0 * pi;
        end
        k = (e2 - e1) - sin(e2 - e1);

        tmaxrp = mu * ( ...
            (2.0 * nrev * pi + k) / (abs(c * c - r * r))^1.5 + ...
            (c * sindeltanu) / (y1 * y2) );

        sme = (r * r - c * c) * 0.5;
        % close to 180 deg transfer case
        nunit = unit(rcrossr);
        if (magr2 * sindeltanu > 0.001)  % 1 meter
            for i = 1 : 3
                nunit(i) = sign(sindeltanu) * nunit(i);
            end
            nr1 = cross(nunit, r1);
            nr2 = cross(nunit, r2);

            for i=1:3
                v1t(i) = (x1 / c) * r1(i) / magr1 + (y1 / c) * (nr1(i) / magr1);
                v2t(i) = (x1 / c) * r2(i) / magr2 + (y2 / c) * (nr2(i) / magr2);
            end
        end

        fprintf(1,'%c  %i  %f   \n',dm, nrev, tmaxrp);  % v1t(1), v1t(2), v1t(3)
    end

end