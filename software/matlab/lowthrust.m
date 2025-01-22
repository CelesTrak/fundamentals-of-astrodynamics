% ------------------------------------------------------------------------------
%
%                           procedure lowthrust
%
%  this procedure calculates the delta v's for a change in semimajor axis and
%  inclination. the orbits are assumed to be circular. finding the lambda
%  value is difficult.
%
%  author        : sal alfano, david vallado      719-573-2600  27 jul 2012
%
%  inputs          description                    range / units
%    ainit       - initial semimajor axis         km
%    afinal      - final semimajor axis           km
%    incl1       - initial inclination            rad
%    incl2       - final inclination              rad
%    isp         - motor isp                      sec
%    mdot        - specific mass flow rate        /s
%    accelinit   - initial acceleration           km/s
%
%  outputs       :
%    deltav      - total change in velocity       km/s
%    tof         - time of flight                 day
%
%  locals        :
%
%  coupling      :
%    lowuz       - analytical method for finding cv
%    pkepler     - simple j2 secular propagation
%
%  references    :
%    vallado       2007, chap 6
% ex 6-12  only a for test
%  [deltav, tof] = lowthrust(1.05*6378.135, 42159.400, 28.5, 28.5, 3800.0, -2.48e-7, 4.0e-3, -0.54);
% ex lowa
%  [deltav, tof] = lowthrust(6378.135+0, 42159.470, 90.0, 179.0, 3800.0, -2.5318398E-08, 9.62099125E-04, -0.54);
% ex 6-13 chg inc so >, both a and i test
%  [deltav, tof] = lowthrust(6378.136+200.000, 42164.000, 28.5, 0.0, 5000.0, -2.0e-7, 1.0e-5, -0.54);
% [deltav, tof] = lowthrust(rinit, rfinal, iinit, ifinal, isp, mdot, accelinit, lambda );
% ----------------------------------------------------------------------------- }

 function [ deltav, tof] = lowthrust(ainit, afinal, incl1, incl2, isp, mdot, accelinit, lambda );
    rad = 180.0 / pi;
    re = 6378.137;
    vkmps = 7.905365719014;
    mu = 3.986004418e5;  %  km^3/s2

    %       numthrusters = 4;
    %       isp = 3800.0; % sec
    %       F = 165.0*0.001*numthrusters;  % N = (m/s) kg / sec
    %       g0 = 10.0; % m/s2
    %       mdotkgs = -F / (isp*g0); % kg/s
    %       mikg = 250.0 + 436.0; % kg
    %       mdotspecific = mdotkgs / mikg; % /s
    % mdot = mdotspecific;
    %       aim = 6378135.0;  % m
    %       afm = 6.61*aim;  % DU of satellite
    %       Acci_ms2 = F/mikg; % m/s2
    % accelinit = Acci_ms2;

    %      fprintf(1,' isp      %11.7f sec  F        %11.7f  N = (m/s) kg / sec \n', isp, F );
    %      fprintf(1,' mdotkgs  %11.7f kg/s aim      %11.7f m \n', mdotkgs, aim );
    %      fprintf(1,' mikg     %11.7f kg   afm      %11.7f DU \n', mikg, afm );
    %      fprintf(1,' mdotspec %11.7f /s   Acci_ms2 %11.7f m/s2 \n', mdotspecific, Acci_ms2 );

    fprintf(1,' isp      %11.7f sec  \n', isp );
    fprintf(1,' mdotspec %11.7f /s   Acci_ms2 %11.7f m/s2 \n', mdot, accelinit );

    iinit = incl1;
    ifinal = incl2;
    ratio = afinal / ainit;
    atrans = (afinal + ainit) * 0.5;  % km
    etrans = (atrans - ainit) / atrans;
    deltai = incl2 - incl1;

    % ---- setup canonical units
    du = ainit;
    tu = sqrt(ainit^3/mu);
    duptu = sqrt(mu/ainit);
        
    atrans = atrans / du; % du
    fprintf(1,' iinit %11.7f ifinal %11.7f ratio %11.7f atrans %11.7f km etrans %11.7f di %11.7f  \n', ...
        incl1, incl2, ratio, atrans*du, etrans, deltai );

    % find lambda
    % lambda = -0.54; % initial guess only!!
    m = 686; % kg
    steps = 360;  % steps per orbit
    % tu = sqrt(ainit^3/mu);

    % form position vector, assuming circular orbit
    r1(1) = ainit;  % km
    r1(2) = 0.0;
    r1(3) = 0.0;
    magr = mag(r1);
    v1(1) = 0.0;
    v1(2) = sqrt(mu/magr) * cos(iinit);  % km/s
    v1(3) = sqrt(mu/magr) * sin(iinit);

    [hbar] = cross( r1,v1 );
    magh= mag( hbar );
    %    hk= hbar(3)/magh;
    %    incl= acos( hk );

    period = 2.0*pi * sqrt(magr^3 / mu);  % in sec
    dtsec = period / steps;

    % start this at zero (the node) this is the angle from the node (perigee)
    % to the current position
    cosf = 1.0;
    cosfold = 0.0;

    afinal = afinal / ainit;  % du
    a1 = 1.0;  % normalize to initial radius du

    t = 0.0; % initial time
    ktr = 1;
    incl1 = incl1/rad;
    rev = 0.0;
nodevec = unit(r1);
%    [p1,a21,ecc1,incl1,omega1,argp1,nu1,m1,arglat1,truelon1,lonper1 ] = rv2coe (r1, v1, mu);
    %    cosforig = truelon1;
  %  cosforig*rad
%a1 < afinal
%incl1 >= ifinal
    while  (a1 < afinal)% && (incl1 >= ifinal) )
        cv = 1.0 / (4.0 * lambda * lambda * a1 * a1 + 1.0);  % canonical

        cosf = lowuz(cv);

        steering = atan( cosf / sqrt(1.0 / (1.0/cv - 1.0 )) );  

        accel = accelinit / (1.0 + mdot*t); % km/s2

        avec = accel * [ 0.0   cos(steering)  sin(steering)];  % switch vnc to ntw

        [rntw1, vntw1, tmntw] = rv2ntw(r1, v1);
        acc1 = tmntw' * avec';  % ntw to  eci

        vtot = v1 + acc1' * dtsec;

        [r1n,v1n] = kepler( r1,vtot, dtsec );  % ndot, nddot
        %        [r1n,v1n] = pkepler( r1,vtot, dtsec, 0.0, 0.0 );  % ndot, nddot

        % use complete classical orbital element transformation
        % you can abbreviate this
    %    [p1,a1,ecc1,incl1,omega1,argp1,nu1,m1,arglat1,truelon1,lonper1 ] = rv2coe (r1n, v1n, mu);

        magr= mag( r1n ); % needed for cosf too
        magv= mag( v1n );
        sme= ( magv*magv*0.5  ) - ( mu /magr );
        a1= -mu  / (2.0 *sme);
        [hbar] = cross( r1n,v1n );
        % sal keeps the oringal hmag, but the current hk component. This lets the
        % inc go down...
%        [hbar] = cross( r1,v1 );
        magh = mag(hbar);
        hk= hbar(3)/magh;
        incl1 = acos( hk );
       
%(argp1+nu1)*rad
        %cosf = cos((argp1+nu1)- cosforig); % switch back   
        cosf = dot(r1n,nodevec)/magr;
        if (cosf*cosfold < 0.0)
            rev = rev + 0.5;
%            fprintf(1,'rev %11.5f a %11.5f DU* i %11.5f \n', rev, a1/ainit, incl1*rad );
        end
        % update total time
        t = t + dtsec;

        r1 = r1n;
%        magr = mag(r1);
        v1 = v1n;
        period = 2.0*pi * sqrt(a1^3/mu);  % in sec
        dtsec = period / steps;

        % alt sals approach, is the same.
        %dtsec = 2*pi*sqrt((a1/du)^3) *tu/steps;  
        cosfold = cosf;

        % save position vectors for plotting
%        r(ktr,1) = r1(1);
%        r(ktr,2) = r1(2);
%        r(ktr,3) = r1(3);

        ktr = ktr + 1;
       if rem(ktr,5000) == 0
           fprintf(1,'%11d  a %11.7f i %11.7f  %11.3f sec  \n',ktr, a1, incl1*rad, t);
       end
        a1 = a1 / du; % get back in du
    end  % while

    deltav = (1.0 - sqrt(1.0 / ratio)) * duptu;  % get in km/s
    tof = t / 86400.0;  % sec

    fprintf(1,' deltav %11.7f km/s  %11.7f du/tu  \n', deltav, deltav/duptu );
    fprintf(1,' tof    %11.7f day \n', tof );
    fprintf(1,' a      %11.7f km \n', a1*du );
    fprintf(1,' incl   %11.7f deg \n', ifinal*rad );

    fprintf(1,' mpfinal %11.7f m final %11.7f kg \n', t*mdot, m + t*mdot );

    %    plot(r(:,1),r(:,2));
    %hold on;
    %plot(x,z);

    %steering(1:100)*rad


