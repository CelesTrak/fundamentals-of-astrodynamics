%  ------------------------------------------------------------------------------
%
%                           procedure lambertminT
%
%  this procedure solves lambert's problem and finds the miniumum time for
%  multi-revolution cases.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    r1          - ijk position vector 1                   km
%    r2          - ijk position vector 2                   km
%    dm          - direction of motion                     'L', 'S'
%    de          - orbital energy                          'L', 'H'
%    nrev        - number of revs to complete              0, 1, 2, 3,
%
%  outputs       :
%    tmin        - minimum time of flight                  sec
%    tminp       - minimum parabolic tof                   sec
%    tminenergy  - minimum energy tof                      sec
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
%    vallado       2022, 481, Alg 57, ex 7-5
%    prussing      JAS 2000
%
%  [tmin, tminp, tminenergy] = lambertminT(r1, r2, dm, de, nrev);
%  ------------------------------------------------------------------------------

function [tmin, tminp, tminenergy] = lambertminT(r1, r2, dm, de, nrev)
    mu = 3.986004415e5;

    magr1 = mag(r1);
    magr2 = mag(r2);
    cosdeltanu = dot(r1, r2) / (magr1 * magr2);
    % make sure it's not more than 1.0
    if (abs(cosdeltanu) > 1.0)
        cosdeltanu = 1.0 * Math.Sign(cosdeltanu);
    end

    rcrossr = cross(r1, r2);
    if (de == 'L')
        sindeltanu = mag(rcrossr) / (magr1 * magr2);
    else
        sindeltanu = -mag(rcrossr) / (magr1 * magr2);
    end
    dnu = atan2(sindeltanu, cosdeltanu);
    % the angle needs to be positive to work for the long way
    if (dnu < 0.0)
        dnu = 2.0 * pi + dnu;
    end

    % ------------- try prussing test his numbers are wrong -------
    %dnu = 75.0 / (180.0 / pi);
    %cosdeltanu = cos(dnu);
    %magr1 = 1.0;
    %magr2 = 1.524;
    %mu = 4.0 * pi * pi;
    %dnu = 90.0 / (180.0 / pi);
    %cosdeltanu = cos(dnu);
    %magr1 = 1.0;  % er
    %magr2 = 1.0;
    %mu = 1.0;

    % these are the same
    %    if (de == 'L')
    chord = sqrt(magr1 * magr1 + magr2 * magr2 - 2.0 * magr1 * magr2 * cosdeltanu);
    %    else
    %    chord = -sqrt(magr1 * magr1 + magr2 * magr2 - 2.0 * magr1 * magr2 * cosdeltanu);
    %    end
    %chord = mag(r2 - r1);

    s = (magr1 + magr2 + chord) * 0.5;

    xi = 0.0;
    eta = 0.0;

    % ------------- calc tmin parabolic tof to see if the orbit is possible
    % ----- no ellitpical orbits exist below this --------
    if (dm == 'S')
        tminp = (1.0 / 3.0) * sqrt(2.0 / mu) * ((s^1.5) - (s - chord)^1.5);
    else
        tminp = (1.0 / 3.0) * sqrt(2.0 / mu) * ((s^1.5) + (s - chord)^1.5);
    end

    % could do this just for nrev cases, but you can also find these for any nrev if (nrev > 0)
    % ------------- this is the min energy ellipse tof
    amin = 0.5 * s;
    %alpha = pi;
    beta = 2.0 * asin(sqrt((s - chord) / s));
    if (dm == 'S')
        tminenergy = (amin^1.5) * ((2.0 * nrev + 1.0) * pi - beta + sin(beta)) / sqrt(mu);
    else
        tminenergy = (amin^1.5) * ((2.0 * nrev + 1.0) * pi + beta - sin(beta)) / sqrt(mu);
    end

    % -------------- calc min tof ellipse (prussing 1992 aas, 2000 jas, stern 1964 pg 230)
    an = 1.001 * amin;
    fa = 10.0;
    i = 1;
    rad = 180.0 / pi;

    while (abs(fa) > 0.00001 && i <= 20)
        a = an;
        alp = 1.0 / a;
        alpha = 2.0 * asin(sqrt(0.5 * s * alp));
        if (de == 'L')
            beta = 2.0 * asin(sqrt(0.5 * (s - chord) * alp));
        else
            beta = - 2.0 * asin(sqrt(0.5 * (s - chord) * alp));  % fix quadrant
        end
        xi = alpha - beta;
        eta = sin(alpha) - sin(beta);
        fa = (6.0 * nrev * pi + 3.0 * xi - eta) * (sin(xi) + eta) - 8.0 * (1.0 - cos(xi));

        fadot = ((6.0 * nrev * pi + 3.0 * xi - eta) * (cos(xi) + cos(alpha)) + ...
            (3.0 - cos(alpha)) * (sin(xi) + eta) - 8.0 * sin(xi)) * (-alp * tan(0.5 * alpha)) ...
            + ((6.0 * nrev * pi + 3.0 * xi - eta) * (-cos(xi) - cos(alpha)) + ...
            (-3.0 - cos(beta)) * (sin(xi) + eta) + 8.0 * sin(xi)) * (-alp * tan(0.5 * beta));
        del = fa / fadot;
        an = a - del;
        %        fprintf(1,'%2i %8.4f %11.5f %11.5f  %11.5f  %11.5f  %11.5f  %11.5f  %11.5f \n',i, dnu*rad, alpha*rad, beta*rad, xi, eta, fa, fadot, an);
        i = i + 1;
    end
    fprintf(1,'iter %2i ',i);
    % could update beta one last time with alpha too????
    if (dm == 'S')
        tmin = (an^1.5) * (2.0 * pi * nrev + xi - eta) / sqrt(mu);
    else
        tmin = (an^1.5) * (2.0 * pi * nrev + xi + eta) / sqrt(mu);
    end

    %      dm = 'S';
    %      de = 'L';
    %      nrev = 0;
    %      blair
    %      r1 = [6778.136300000, 0.000000, 0.000000 ];
    %      r2 = [-6694.857334274, -1180.483980026, 0.000000 ];
    %      v1 = [0.000000, 7.668558568, 0.000000 ];
    %      moving
    %      r1 = [ -6175.1034, 2757.0706, 1626.6556 ];
    %      v1 = [ 2.376641, 1.139677, 7.078097];
    %      r2 = [ -1078.007289, 8796.641859, 1890.7135 ];
    %
    %      [tmin, tminp, tminenergy] = lambertminT(r1, r2, dm, de, nrev);
    %
    fprintf(1,'%c  %c %i  %f  %f  %f  \n',dm, de, nrev, tmin, tminp, tminenergy);

end