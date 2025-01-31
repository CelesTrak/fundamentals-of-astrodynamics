% -----------------------------------------------------------------------------
%
%                           procedure gibbs
%
%  this procedure performs the gibbs method of orbit determination.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    r1          -  position vector #1                        km
%    r2          -  position vector #2                        km
%    r3          -  position vector #3                        km
%
%  outputs       :
%    v2          -  velocity vector for r2                    km / s
%    theta       - angle between vectors                         rad
%    error       - flag indicating success                       'ok',...
%
%  locals        :
%    tover2      -
%    l           -
%    small       - tolerance for roundoff errors
%    r1mr2       - MathTimeLibr.magnitude of r1 - r2
%    r3mr1       - MathTimeLibr.magnitude of r3 - r1
%    r2mr3       - MathTimeLibr.magnitude of r2 - r3
%    p           - p vector     r2 x r3
%    q           - q vector     r3 x r1
%    w           - w vector     r1 x r2
%    d           - d vector     p + q + w
%    n           - n vector (r1)p + (r2)q + (r3)w
%    s           - s vector
%                    (r2-r3)r1+(r3-r1)r2+(r1-r2)r3
%    b           - b vector     d x r2
%    theta1      - temp angle between the vectors                rad
%    pn          - p unit vector
%    r1n         - r1 unit vector
%    dn          - d unit vector
%    nn          - n unit vector
%
%  coupling      :
%    mag         - magnitude of a vector
%    cross       - cross product of two vectors
%    dot         - dot product of two vectors
%    add3vec     - add three vectors
%    norm        - creates a unit vector
%    angle       - angle between two vectors
%
%  references    :
%    vallado       2022, 460, alg 54, ex 7-3
%
% [v2, theta,theta1,copa, error] = function gibbs( r1,r2,r3);
% ------------------------------------------------------------------------------

function [v2, theta,theta1,copa, error] = gibbs( r1,r2,r3)

    % -------------------------  implementation   -------------------------
    constmath;
    constastro;

    small= 0.000001;
    theta= 0.0;
    error = '          ok';
    theta1= 0.0;

    magr1 = mag( r1 );
    magr2 = mag( r2 );
    magr3 = mag( r3 );
    for i= 1 : 3
        v2(i)= 0.0;
    end

    p = cross( r2,r3 );
    q = cross( r3,r1 );
    w = cross( r1,r2 );
    pn = unit( p );
    r1n = unit( r1 );
    copa=  asin( dot( pn,r1n ) );

    if ( abs( copa ) > 0.017452406 )
        error= 'not coplanar';
    end

    % --------------- | don't continue processing --------------
    d = p + q + w;
    magd = mag(d);
    n = magr1*p + magr2*q + magr3*w;
    magn = mag(n);
    nn = unit( n );
    dn = unit( d );

    % -------------------------------------------------------------
    %       determine if  the orbit is possible.  both d and n must be in
    %         the same direction, and non-zero.
    % -------------------------------------------------------------
    if ( ( abs(magd)<small ) | ( abs(magn)<small ) | ...
            ( dot(nn,dn) < small ) )
        error= '  impossible';
    else
        theta  = angl( r1,r2 );
        theta1 = angl( r2,r3 );

        % ----------- perform gibbs method to find v2 -----------
        r1mr2= magr1-magr2;
        r3mr1= magr3-magr1;
        r2mr3= magr2-magr3;
        s  = r1mr2*r3 + r3mr1*r2 + r2mr3*r1;
        b  = cross( d,r2 );
        l  = sqrt(mu / (magd*magn) );
        tover2= l / magr2;
        v2 = tover2 * b + l * s;
    end

    %    fprintf(1,'p     %11.7f   %11.7f  %11.7f km2 \n',p);
    %    fprintf(1,'n     %11.7f   %11.7f  %11.7f km3 \n',n);
    %    fprintf(1,'d     %11.7f   %11.7f  %11.7f km3 \n',d);
    %    fprintf(1,'s     %11.7f   %11.7f  %11.7f km2 \n',s);
    %    fprintf(1,'theta     %11.7f   %11.7f deg \n',theta*180/pi, theta1*180/pi );
    %    fprintf(1,'b     %11.7f   %11.7f  %11.7f km3 \n',b);
    %    fprintf(1,'l     %11.7f  /kms \n',l);

end