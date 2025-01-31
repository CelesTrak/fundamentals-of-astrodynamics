% ----------------------------------------------------------------------------
%
%                           function covfl2ct
%
%  this function transforms a six by six covariance matrix expressed in
%    flight elements into one expressed in cartesian elements.
%
%  author        : david vallado                  719-573-2600   27 may 2003
%
%  revisions
%    vallado     - major update                                  26 aug 2015
%
%  inputs          description                    range / units
%    flcov       - 6x6 flight covariance matrix
%    flstate     - 6x1 flight orbit state         (r v latgc lon fpa az)
%    anom        - anomaly                        'latlon', 'radec'
%    ttt         - julian centuries of tt         centuries
%    jdut1       - julian date of ut1             days from 4713 bc
%    lod         - excess length of day           sec
%    xp          - polar motion coefficient       arc sec
%    yp          - polar motion coefficient       arc sec
%    terms       - number of terms for ast calculation 0,2
%
%  outputs       :
%    cartcov     - 6x6 cartesian covariance matrix
%    tm          - transformation matrix
%
%  locals        :
%    tm           - matrix of partial derivatives
%    magr        - eci position vector magnitude  km
%    magv        - eci velocity vector magnitude  km/sec
%    latgc       - geocentric latitude            rad
%    lon         - longitude                      rad
%    fpa         - sat flight path angle          rad
%    az          - sat flight path az             rad
%    fpav        - sat flight path anglefrom vert rad
%    xe,ye,ze    - ecef position vector componentskm
%
%  coupling      :
%    ecef2eci    - convert eci vectors to ecef
%
%  references    :
%    Vallado and Alfano 2015
%
% [cartcov,tm] = covfl2ct( flcov,flstate, anomflt, iau80arr, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps);
% ----------------------------------------------------------------------------

function [cartcov,tm] = covfl2ct( flcov,flstate, anomflt, iau80arr, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps)

    small  = 0.00000001;

    % -------- parse the input vectors into components
    lon  = flstate(1); % these will come in as either lon/lat or rtasc/decl depending on anom1
    latgc= flstate(2);
    fpa  = flstate(3);
    az   = flstate(4);
    magr = flstate(5);  % already converted to m in setcov
    magv = flstate(6);

    cfpa = cos(fpa);
    sfpa = sin(fpa);
    caz = cos(az);
    saz = sin(az);

    % --------- determine which set of variables is in use ---------
    % need to get eci vector
    % seems like it's easier to simply convert the lat-lon to rtasc
    % decl
    if strcmp(anomflt, 'latlon') == 1
        craf = cos(lon);  % earth fixed needed for the lon lat partials only
        sraf = sin(lon);
        cdf = cos(latgc);
        sdf = sin(latgc);
        recef(1) = magr*0.001*cos(latgc)*cos(lon);  % in km
        recef(2) = magr*0.001*cos(latgc)*sin(lon);
        recef(3) = magr*0.001*sin(latgc);
        % -------- convert r to eci
        % this vel is wrong but not needed except for special case ahead
        vecef(1) = magv*0.001*( -cos(lon)*sin(latgc)*caz*cfpa - sin(lon)*saz*cfpa + cos(lon)*cos(latgc)*sfpa ); % m/s
        vecef(2) = magv*0.001*( -sin(lon)*sin(latgc)*caz*cfpa + cos(lon)*saz*cfpa + sin(lon)*cos(latgc)*sfpa );
        vecef(3) = magv*0.001*( sin(lon)*sfpa + cos(latgc)*caz*cfpa );
        aecef = [0.0; 0.0; 0.0];
        [reci, veci, aeci] = ecef2eci(recef', vecef', aecef, iau80arr, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps );
        reci = reci*1000.0;  % in m
        veci = veci*1000.0;  % in m/s

        temp= sqrt( reci(1)*reci(1) + reci(2)*reci(2) );
        if ( temp < small )
            rtasc= atan2( veci(2) , veci(1) );
            temp
        else
            rtasc= atan2( reci(2) , reci(1) );
        end
        %decl = atan2( reci(3) , sqrt(reci(1)^2 + reci(2)^2) )
        decl = asin( reci(3)/magr );
        cra = cos(rtasc);
        sra = sin(rtasc);
        cd = cos(decl);
        sd = sin(decl);
    else
        if strcmp(anomflt, 'radec') == 1
            rtasc = lon;  % these come in as rtasc decl in this case
            decl = latgc;
            reci(1) = magr*cos(decl)*cos(rtasc);
            reci(2) = magr*cos(decl)*sin(rtasc);
            reci(3) = magr*sin(decl);
            veci(1) = magv*( -cos(rtasc)*sin(decl)*caz*cfpa - sin(rtasc)*saz*cfpa + cos(rtasc)*cos(decl)*sfpa ); % m/s
            veci(2) = magv*( -sin(rtasc)*sin(decl)*caz*cfpa + cos(rtasc)*saz*cfpa + sin(rtasc)*cos(decl)*sfpa );
            veci(3) = magv*( sin(decl)*sfpa + cos(decl)*caz*cfpa );
            cra = cos(rtasc);
            sra = sin(rtasc);
            cd = cos(decl);
            sd = sin(decl);
        end
    end

    % ---------------- calculate matrix elements ------------------
    % ---- partials of rx wrt (lon latgc fpa az r v)
    if strcmp(anomflt, 'radec') == 1
        tm(1,1) = -magr*cd*sra;
        tm(1,2) = -magr*sd*cra;
    else  % latlon
        tm(1,1) = -magr*cdf*sraf;
        tm(1,2) = -magr*sdf*craf;
    end
    tm(1,3) = 0.0;
    tm(1,4) = 0.0;
    tm(1,5) = cd*cra;
    tm(1,6) = 0.0;

    % ---- partials of ry wrt (lon latgc fpa az r v)
    if strcmp(anomflt, 'radec') == 1
        tm(2,1) =   magr*cd*cra;
        tm(2,2) =  -magr*sd*sra;
    else   % latlon
        tm(2,1) =   magr*cdf*craf;
        tm(2,2) =  -magr*sdf*sraf;
    end
    tm(2,3) =  0.0;
    tm(2,4) =  0.0;
    tm(2,5) =  cd*sra;
    tm(2,6) =  0.0;

    % ---- partials of rz wrt (lon latgc fpa az r v)
    if strcmp(anomflt, 'radec') == 1
        tm(3,1) =  0.0;
        tm(3,2) =  magr*cd;
    else   % latlon
        tm(3,1) =  0.0;
        tm(3,2) =  magr*cdf;
    end
    tm(3,3) =  0.0;
    tm(3,4) =  0.0;
    tm(3,5) =  sd;
    tm(3,6) =  0.0;

    % ---- partials of vx wrt (lon latgc fpa az r v)
    if strcmp(anomflt, 'radec') == 1
        tm(4,1) =  -magv*( -sra*caz*sd*cfpa + cra*saz*cfpa + cd*sra*sfpa );
        %  tm(4,1) = -vy;
        tm(4,2) =  -cra*magv*( sd*sfpa + cd*caz*cfpa );
        %  tm(4,2) = -vz*cra;
    else   % latlon

        tm(4,1) =  -magv*( -sraf*caz*sdf*cfpa + craf*saz*cfpa + cdf*sraf*sfpa );
        %  tm(4,1) = -vy;
        tm(4,2) =  -craf*magv*( sdf*sfpa + cdf*caz*cfpa );
        %  tm(4,2) = -vz*cra;
    end
    tm(4,3) =  magv*(  cra*caz*sd*sfpa + sra*saz*sfpa + cd*cra*cfpa );
    tm(4,4) =  magv*(  cra*saz*sd*cfpa - sra*caz*cfpa );
    tm(4,5) =  0.0;
    tm(4,6) =  -cra*caz*sd*cfpa - sra*saz*cfpa + cd*cra*sfpa;

    % ---- partials of vy wrt (lon latgc fpa az r v)
    if strcmp(anomflt, 'radec') == 1
        tm(5,1) =  magv*(-cra*caz*sd*cfpa - sra*saz*cfpa + cd*cra*sfpa);
        %  tm(5,1) = vx;
        tm(5,2) =  -sra*magv*( sd*sfpa + cd*caz*cfpa );
        %   tm(5,2) = -vz*sra;
    else   % latlon
        tm(5,1) =  magv*(-craf*caz*sdf*cfpa - sraf*saz*cfpa + cdf*craf*sfpa);
        %  tm(5,1) = vx;
        tm(5,2) =  -sraf*magv*( sdf*sfpa + cdf*caz*cfpa );
        %   tm(5,2) = -vz*sra;
    end
    tm(5,3) =  magv*( sra*caz*sd*sfpa - cra*saz*sfpa + cd*sra*cfpa);
    tm(5,4) =  magv*( sra*saz*sd*cfpa + cra*caz*cfpa);
    tm(5,5) =  0.0;
    tm(5,6) =  -sra*caz*sd*cfpa + cra*saz*cfpa + cd*sra*sfpa;

    % ---- partials of vz wrt (lon latgc fpa az r v)
    if strcmp(anomflt, 'radec') == 1
        tm(6,1) =  0.0;
        tm(6,2) =  magv*(cd*sfpa - sd*caz*cfpa);
    else   % latlon
        tm(6,1) =  0.0;
        tm(6,2) =  magv*(cdf*sfpa - sdf*caz*cfpa);
    end
    tm(6,3) =  magv*(sd*cfpa - cd*caz*sfpa);
    tm(6,4) =  -magv*cd*saz*cfpa;
    tm(6,5) =  0.0;
    tm(6,6) =  sd*sfpa + cd*caz*cfpa;

    % ---------- calculate the output covariance matrix -----------
    [cartcov] =  tm*flcov*tm';


