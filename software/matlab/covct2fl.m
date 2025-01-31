% ----------------------------------------------------------------------------
%
%                           function covct2fl
%
%  this function transforms a six by six covariance matrix expressed in cartesian elements
%    into one expressed in flight parameters
%
%  author        : david vallado                  719-573-2600   21 jun 2002
%
%  inputs          description                    range / units
%    cartcov     - 6x6 cartesian covariance matrix
%    cartstate   - 6x1 cartesian orbit state      (x y z vx vy vz)
%    anom        - anomaly                        'latlon', 'radec'
%    ttt         - julian centuries of tt         centuries
%    jdut1       - julian date of ut1             days from 4713 bc
%    lod         - excess length of day           sec
%    xp          - polar motion coefficient       arc sec
%    yp          - polar motion coefficient       arc sec
%    terms       - number of terms for ast calculation 0,2
%
%  outputs       :
%    flcov       - 6x6 flight covariance matrix
%    tm          - transformation matrix
%
%  locals        :
%    r           - matrix of partial derivatives
%    x,y,z       - components of position vector  km
%    vx,vy,vz    - components of position vector  km/s
%    magr        - eci position vector magnitude  km
%    magv        - eci velocity vector magnitude  km/sec
%    d           - r dot v
%    h           - angular momentum vector
%    hx,hy,hz    - components of angular momentum vector
%    hcrossrx,y,z- components of h cross r vector
%    p1,p2       - denominator terms for the partials
%
%  coupling      :
%    ecef2eci    - convert eci vectors to ecef
%
%  references    :
%    Vallado and Alfano 2015
%
%   [flcov,tm] = covct2fl( cartcov,cartstate, anomflt, iau80arr, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps);
% ----------------------------------------------------------------------------

function [flcov,tm] = covct2fl( cartcov,cartstate, anomflt, iau80arr, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps)

    % -------- parse the input vectors into cartesian components
    rx = cartstate(1) * 1000;  % keep all in m, m/s
    ry = cartstate(2) * 1000;  % this is eci always
    rz = cartstate(3) * 1000;
    vx = cartstate(4) * 1000;
    vy = cartstate(5) * 1000;
    vz = cartstate(6) * 1000;

    if strcmp(anomflt, 'latlon') == 1
        % -------- convert r to eci
        reci = [rx;ry;rz]/1000;
        veci = [vx;vy;vz]/1000;
        aeci = [0.0; 0.0; 0.0];
        % now find ecef coordinates
        [recef, vecef, aecef] = eci2ecef(reci, veci, aeci, iau80arr, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps );

        recef = recef * 1000;  % m
        vecef = vecef * 1000;
        rxf = recef(1);
        ryf = recef(2);
        rzf = recef(3);
        % [rxf; ryf; rzf]
    else
        if strcmp(anomflt, 'radec') == 1
        end
    end

    % -------- calculate common quantities
    magr = sqrt(rx^2 + ry^2 + rz^2);   % in m
    magv = sqrt(vx^2 + vy^2 + vz^2);
    rdotv  = rx*vx + ry*vy + rz*vz;

    h  = sqrt((rx*vy-ry*vx)^2 + (rz*vx-rx*vz)^2 + (ry*vz-rz*vy)^2);
    %          hx = ry*vz - rz*vy;  %
    %          hy = rz*vx - rx*vz;
    %          hz = rx*vy - ry*vx;
    %          hcrossrx = (rz*hy - ry*hz);
    %          hcrossry = (rx*hz - rz*hx);
    %          hcrossrz = (ry*hx - rx*hy);

    % ---------------- calculate matrix elements ------------------
    if strcmp(anomflt, 'latlon') == 1
        % partial of lon wrt (x y z vx vy vz)
        p0 = 1.0 / (rxf^2 + ryf^2);
        tm(1,1) = -p0*ryf;
        tm(1,2) =  p0*rxf;
        tm(1,3) = 0.0;
        tm(1,4) = 0.0;
        tm(1,5) = 0.0;
        tm(1,6) = 0.0;

        % partial of latgc wrt (x y z vx vy vz)
        p0 = 1.0 / (magr^2*sqrt(rxf^2 + ryf^2));
        tm(2,1) = -p0*(rxf*rzf);
        tm(2,2) = -p0*(ryf*rzf);
        tm(2,3) = sqrt(rxf^2 + ryf^2) / magr^2;
        tm(2,4) = 0.0;
        tm(2,5) = 0.0;
        tm(2,6) = 0.0;
    else
        if strcmp(anomflt, 'radec') == 1
            % partial of lon wrt (x y z vx vy vz)
            p0 = 1.0 / (rx^2 + ry^2);
            tm(1,1) = -p0*ry;
            tm(1,2) =  p0*rx;
            tm(1,3) = 0.0;
            tm(1,4) = 0.0;
            tm(1,5) = 0.0;
            tm(1,6) = 0.0;

            % partial of latgc wrt (x y z vx vy vz)
            p0 = 1.0 / (magr^2*sqrt(rx^2 + ry^2));
            tm(2,1) = -p0*(rx*rz);
            tm(2,2) = -p0*(ry*rz);
            tm(2,3) = sqrt(rx^2 + ry^2) / magr^2;
            tm(2,4) = 0.0;
            tm(2,5) = 0.0;
            tm(2,6) = 0.0;
        end
    end

    % partial of fpa wrt (x y z vx vy vz)
    rdot = rdotv / magr;  % (r dot v) / r
    %        p1 = -1.0 / (magr*sqrt(magv^2 - rdot^2));
    %         tm(3,1) = p1*( rdot*rx/magr - vx );
    %         tm(3,2) = p1*( rdot*ry/magr - vy );
    %         tm(3,3) = p1*( rdot*rz/magr - vz );
    %         tm(3,4) = p1*( rdot*magr*vx/magv^2 - rx );
    %         tm(3,5) = p1*( rdot*magr*vy/magv^2 - ry );
    %         tm(3,6) = p1*( rdot*magr*vz/magv^2 - rz );
    % Sal from mathcad matches previous with - sign on previous
    p0 = 1.0 / (magr^2*h);
    p1 = 1.0 / (magv^2*h);
    tm(3,1) = p0*( vx*(ry^2 + rz^2) - rx*(ry*vy + rz*vz) );
    tm(3,2) = p0*( vy*(rx^2 + rz^2) - ry*(rx*vx + rz*vz) );
    tm(3,3) = p0*( vz*(rx^2 + ry^2) - rz*(rx*vx + ry*vy) );
    tm(3,4) = p1*( rx*(vy^2 + vz^2) - vx*(ry*vy + rz*vz) );
    tm(3,5) = p1*( ry*(vx^2 + vz^2) - vy*(rx*vx + rz*vz) );
    tm(3,6) = p1*( rz*(vx^2 + vy^2) - vz*(rx*vx + ry*vy) );

    % partial of az wrt (x y z vx vy vz)
    %         p2 = 1.0 / ((magv^2 - rdot^2) * (rx^2 + ry^2));
    %         tm(4,1) = p2*( vy*(magr*vz - rz*rdot) - (rx*vy - ry*vx) * (rx*vz - rz*vx + rx*rz*rdot/magr) * (1.0 / magr) );
    %         tm(4,2) = p2*( -vx*(magr*vz - rz*rdot) + (rx*vy - ry*vx) * (ry*vz - rz*vy + ry*rz*rdot/magr) * (1.0 / magr) );
    %         p2 = 1.0 / (magr^2 * (magv^2 - rdot^2));
    %         tm(4,3) = p2 * rdot * (rx*vy - ry*vx);
    %         p2 = 1.0 / (magr * (magv^2 - rdot^2));
    %         tm(4,4) = -p2 * (ry*vz - rz*vy);
    %         tm(4,5) =  p2 * (rx*vz - rz*vx);
    %         tm(4,6) = -p2 * (rx*vy - ry*vx);
    % sal from mathcad
    p2 = 1.0 / ((magv^2 - rdot^2)*(rx^2 + ry^2));
    k1 = sqrt(rx^2 + ry^2 + rz^2)*(rx*vy - ry*vx);
    k2 = ry*(ry*vz - rz*vy) + rx*(rx*vz - rz*vx);
    tm(4,1) = p2 * ( vy*(magr*vz - rz*rdot) - (rx*vy - ry*vx)/magr * (rx*vz - rz*vx + rx*ry*rdot/magr) );
    p2 = 1.0 / (magr*(k1^2 + k2^2));
    tm(4,1) = p2 * (k1*magr*(rz*vx - 2*rx*vz) + k2*(-ry*vx*rx + vy*rx^2 + vy*magr^2));
    tm(4,2) = p2 * (k1*magr*(rz*vy - 2*ry*vz) + k2*(rx*vy*ry - vx*ry^2 - vx*magr^2));
    p2 = k1 / (magr^2*(k1^2 + k2^2));
    tm(4,3) = p2 * (k2*rz + (rx*vx + ry*vy)*magr^2);
    p2 = 1.0 / (k1^2 + k2^2);
    tm(4,4) = p2 * (k1*rx*rz - k2*ry*magr);
    tm(4,5) = p2 * (k1*ry*rz + k2*rx*magr);
    tm(4,6) = -p2 * (k1*(rx^2 + ry^2));

    % partial of r wrt (x y z vx vy vz)
    p0 = 1.0 / magr;
    tm(5,1) = p0*rx;
    tm(5,2) = p0*ry;
    tm(5,3) = p0*rz;
    tm(5,4) = 0.0;
    tm(5,5) = 0.0;
    tm(5,6) = 0.0;

    % partial of v wrt (x y z vx vy vz)
    p0 = 1.0 / magv;
    tm(6,1) = 0.0;
    tm(6,2) = 0.0;
    tm(6,3) = 0.0;
    tm(6,4) = p0*vx;
    tm(6,5) = p0*vy;
    tm(6,6) = p0*vz;

    % ---------- calculate the output covariance matrix -----------
    [flcov]= tm*cartcov*tm';

