% -----------------------------------------------------------------------------
%
%                           function findfandg
%
%  this function calculates the f and g functions for use in various applications.
%  several methods are available. the values are in normal (not canonical) units.
%  note that not all the input parameters are needed for each case. also, the step
%  size dtsec should be small, perhaps on the order of 60-120 secs!
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    r1          - position vector                          km
%    v1          - velocity vector                          km/s
%    r2          - position vector                          km
%    v2          - velocity vector                          km/s
%    x           - universal variable
%    c2          - stumpff function
%    c3          - stumpff function
%    dtsec       - step size                                sec (SMALL time steps only!!)
%    opt         - calculation method                       pqw, series, c2c3
%
%  outputs       :
%    f, g        - f and g functions
%    fdot, gdot  - fdot and gdot functions
%
%  locals        :
%                -
%  coupling      :
%
%  references    :
%    vallado       2022, 83, 87, 813
%
%  [f, g, fdot, gdot] = findfandg (r1, v1, r2, v2, dtsec, x, c2, c3, z, opt);
% -----------------------------------------------------------------------------

function [f, g, fdot, gdot] = findfandg (r1, v1, r2, v2, dtsec, x, c2, c3, z, opt)
    constastro;
    f = 0.0;
    g = 0.0;
    gdot = 0.0;

    magr1 = mag(r1);
    magv1 = mag(v1);

    % -------------------------  implementation    % ----------------
    switch (opt)

        case 'pqw'
            [hbar] = cross(r1, v1);
            h = mag(hbar);
            % find vectors in PQW frame
            [rpqw1, vpqw1] = rv2pqw(r1, v1);
            % find vectors in PQW frame
            [rpqw2, vpqw2] = rv2pqw(r2, v2);

            % normal units
            f = (rpqw2(1) * vpqw1(2) - vpqw2(1) * rpqw1(2)) / h;
            g = (rpqw1(1) * rpqw2(2) - rpqw2(1) * rpqw1(2)) / h;
            gdot = (rpqw1(1) * vpqw2(2) - vpqw2(1) * rpqw1(2)) / h;
            fdot = (vpqw2(1) * vpqw1(2) - vpqw2(2) * vpqw1(1)) / h;
        case 'series'
            u = mu / (magr1 * magr1 * magr1);
            p = dot(r1, v1) / (magr1 * magr1);
            q = (magv1 * magv1 - magr1 * magr1 * u) / (magr1 * magr1);
            p2 = p * p;
            p3 = p2 * p;
            p4 = p3 * p;
            p5 = p4 * p;
            p6 = p5 * p;
            u2 = u * u;
            u3 = u2 * u;
            u4 = u3 * u;
            u5 = u4 * u;
            q2 = q * q;
            q3 = q2 * q;
            dt2 = dtsec * dtsec;
            dt3 = dt2 * dtsec;
            dt4 = dt3 * dtsec;
            dt5 = dt4 * dtsec;
            dt6 = dt5 * dtsec;
            dt7 = dt6 * dtsec;
            dt8 = dt7 * dtsec;
            f = 1.0 - 0.5 * u * dt2 + 0.5 * u * p * dt3 ...
                + u / 24.0 * (-15 * p2 + 3.0 * q + u) * dt4 ...
                + p * u / 8.0 * (7.0 * p2 - 3 * q - u) * dt5 ...
                + u / 720.0 * (-945.0 * p4 + 630.0 * p2 * q + 210 * u * p2 - 45 * q2 - 24 * u * q - u2) * dt6 ...
                + p * u / 80.0 * (165 * p4 - 150 * p2 * q - 50 * u * p2 + 25 * q2 + 14.0 * u * q + u2) * dt7 ...
                + u / 40320.0 * (-135135.0 * p6 + 155925.0 * p4 * q + 51975.0 * u * p4 - 42525 * p2 * q2 ...
                - 24570.0 * u * p2 * q - 2205 * u2 * p2 + 1575 * q3 + 1107.0 * u * q2 + 117.0 * u2 * q + u3) * dt8;

            g = dtsec - 1.0 / 6.0 * u * dt3 + 0.25 * u * p * dt4 ...
                + u / 120.0 * (-45 * p2 + 9.0 * q + u) * dt5 ...
                + p * u / 24.0 * (14.0 * p2 - 6 * q - u) * dt6 ...
                + u / 5040.0 * (-4725 * p4 + 3150 * p2 * q + 630 * u * p2 - 225 * q2 - 54 * u * q - u2) * dt7 ...
                + p * u / 320.0 * (495 * p4 - 450 * p2 * q - 100 * u * p2 + 75 * q2 + 24.0 * u * q + u2) * dt8;

            fdot = -u * dtsec + 1.5 * u * p * dt2 ...
                + u / 6.0 * (-15 * p2 + 3 * q + u) * dt3 ...
                + 5.0 * p * u / 8.0 * (7.0 * p2 - 3 * q - u) * dt4 ...
                + u / 120.0 * (-945 * p4 + 630 * p2 * q + 210 * u * p2 - 45 * q2 - 24 * u * q - u2) * dt5 ...
                + 7.0 * p * u / 80.0 * (165 * p4 - 150 * p2 * q - 50 * u * p2 + 25 * q2 + 14 * u * q + u2) * dt6 ...
                + u / 5040.0 * (-135135.0 * p6 + 155925.0 * p4 * q + 51975.0 * u * p4 - 42525.0 * p2 * q2 - 24570.0 * u * p2 * q ...
                - 2205.0 * u2 * p2 + 1575 * q3 + 1107 * u * q2 + 117 * u2 * q + u3) * dt7;

            gdot = 1.0 - 0.5 * u * dt2 + u * p * dt3 ...
                + u / 24 * (-45 * p2 + 9 * q + u) * dt4 ...
                + p * u / 4 * (14 * p2 - 6 * q - u) * dt5 ...
                + u / 720 * (-4725.0 * p4 + 3150 * p2 * q + 630 * u * p2 - 225 * q2 - 54 * u * q - u2) * dt6 ...
                + p * u / 40 * (495 * p4 - 450 * p2 * q - 100 * u * p2 + 75 * q2 + 24 * u * q + u2) * dt7;

            %f = dt8* u* (-135135 * p6 + 155925 * p4* q + 51975 * p4* u - 42525 * p2 * q2 - 24570 * p2 * q * u - 2205 * p2 * u2  ...
            %+ 1575 * q3 + 1107 * q2 * u + 117 * q * u2 + u3) / 40320 + dt7 * p * u * (165 * p4 - 150 * p2 * q - 50 * p2 * u  ...
            %+ 25 * q2 + 14 * q * u + u2) / 80 + dt6* u * (-945 * p4 + 630 * p2 * q + 210 * p2 * u - 45 * q2 - 24 * q * u - u2) / ...
            %720 + dt5 * p * u * (7 * p2 - 3 * q - u) / 8 + dt4 * u * (-15 * p2 + 3 * q + u) / 24 + dt3 * p * u / 2 - dt2 * u / 2 + 1;
            %g = dt8 * p * u * (495 * p4 - 450 * p2 * q - 100 * p2 * u + 75 * q2 + 24 * q * u + u2) / 320 + dt7 * u * (-4725 * p4  ...
            %+ 3150 * p2 * q + 630 * p2 * u - 225 * q2 - 54 * q * u - u2) / 5040 + dt6 * p * u * (14 * p2 - 6 * q - u) / 24 + dt5 * u * (-45 * p2 + 9 * q + u) / 120  ...
            %+ dt4 * p * u / 4 - dt3 * u / 6 + dtsec;
            %fdot = dt7 * u * (-135135 * p6 + 155925 * p4 * q + 51975 * p4 * u - 42525 * p2 * q2 - 24570 * p2 * q * u - 2205 * p2 * u2 + 1575 * q3 + 1107 * q2 * u
            %+ 117 * q * u2 + u3) / 5040 + 7 * dt6 * p * u * (165 * p4 - 150 * p2 * q - 50 * p2 * u + 25 * q2 + 14 * q * u + u2) / 80 + dt5 * u * (-945 * p4 + 630 * p2 * q
            %+ 210 * p2 * u - 45 * q2 - 24 * q * u - u2) / 120 + 5 * dt4 * p * u * (7 * p2 - 3 * q - u) / 8 + dt3 * u * (-15 * p2 + 3 * q + u) / 6 + 3 * dt2 * p * u / 2 - dtsec * u;
            %gdot = dt7 * p * u * (495 * p4 - 450 * p2 * q - 100 * p2 * u + 75 * q2 + 24 * q * u + u2) / 40 + dt6 * u * (-4725 * p4 + 3150 * p2 * q + 630 * p2 * u - 225 * q2 - 54 * q * u - u2) / 720
            %+ dt5 * p * u * (14 * p2 - 6 * q - u) / 4 + dt4 * u * (-45 * p2 + 9 * q + u) / 24 + dt3 * p * u - dt2 * u / 2 + 1;

        case 'c2c3'
            xsqrd = x * x;
            magr2 = mag(r2);
            f = 1.0 - (xsqrd * c2 / magr1);
            g = dtsec - xsqrd * x * c3 / sqrt(mu);
            gdot = 1.0 - (xsqrd * c2 / magr2);
            fdot = (sqrt(mu) * x / (magr1 * magr2)) * (z * c3 - 1.0);
        otherwise
            f = 0.0;
            g = 0.0;
            fdot = 0.0;
            gdot = 0.0;
    end
    % g = g * tusec;  / to canonical if needed, fdot/tusec too
    %fdot = (f * gdot - 1.0) / g;

end  %  findfandg

