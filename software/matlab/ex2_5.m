%     -----------------------------------------------------------------
%
%                              Ex2_5.m
%
%  this file demonstrates example 2-5. it also includes some stressing
%  cases for the coe and rv conversions for all orbit types.
%
%                          companion code for
%             fundamentals of astrodynamics and applications
%                                 2013
%                            by david vallado
%
%     (h)               email davallado@gmail.com
%     (w) 719-573-2600, email dvallado@agi.com
%
%     *****************************************************************
%
%  current :
%            16 feb 19  david vallado
%                         update for new constants
%  changes :
%            13 feb 07  david vallado
%                         original baseline
%
%     *****************************************************************

    rad = 180.0 /pi;
    
    r = [ -605.79221660; -5870.22951108; 3493.05319896 ];
    v = [ -1.56825429; -3.70234891; -6.47948395 ];
    [p,a,ecc,incl,raan,argp,nu,m,arglat,truelon,lonper ] = rv2coe(r,v);
    
    [r, v] = coe2rv(p, ecc, incl, raan, argp, nu, arglat, truelon, lonper);
    
    fprintf(1,'start %15.9f %15.9f %15.9f',r );
    fprintf(1,' v %15.10f %15.10f %15.10f\n',v );
    
    [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe (r,v);
    
    fprintf(1,'          p km       a km      ecc      incl deg     raan deg     argp deg      nu deg      m deg      arglat   truelon    lonper\n');
    fprintf(1,'coes %11.4f %11.4f %13.9f %13.7f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f\n',...
        p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad, ...
        arglat*rad,truelon*rad,lonper*rad );
    
    
    
    fprintf(1,'\ncoe test ----------------------------\n' );
    r=[ 6524.834;6862.875;6448.296];
    v=[ 4.901327;5.533756;-1.976341];
    
    fprintf(1,'start %15.9f %15.9f %15.9f',r );
    fprintf(1,' v %15.10f %15.10f %15.10f\n',v );
    
    % --------  rv2coe       - position and velocity vectors to classical elements
    
    [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe(r,v);
    
    fprintf(1,'          p km       a km      ecc      incl deg     raan deg     argp deg      nu deg      m deg      arglat   truelon    lonper\n');
    fprintf(1,'coes %11.4f %11.4f %13.9f %13.7f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f\n',...
        p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad, ...
        arglat*rad,truelon*rad,lonper*rad );
    
    
    


