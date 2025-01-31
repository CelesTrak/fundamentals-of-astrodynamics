% ----------------------------------------------------------------------------
%
%                           function iau06xys
%
%  this function calculates the transformation matrix that accounts for the
%    effects of precession-nutation in the iau2006 cio theory.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    jdtt        - julian date in tt                           days from 4713 bc
%    jdftt       - fractional julian date in tt                     days
%    ddx         - delta x correction to gcrf                       rad
%    ddy         - delta y correction to gcrf                       rad
%    interp      - interpolation type (x for full series)           x, n, l, s
%                  none, linear, spline
%    fArgs06     - fundamental arguments in an array
%    iau06arr    - iau06 array of values
%    xysarr      - iau06 array of xys values for interpolation
%
%  outputs       :
%    iau06xys    - transformation matrix for itrf-gcrf
%    x           - coordinate of cip                                rad
%    y           - coordinate of cip                                rad
%    s           - coordinate                                       rad
%
%  locals        :
%    ttt         - julian centuries of tt
%
%  coupling      :
%
%  references    :
%    vallado       2022, 214, 221
%
% [x, y, s, pn] = iau06xys (iau06arr, fArgs, ttt, ddx, ddy, opt1);
% ----------------------------------------------------------------------------

function [x, y, s, pn] = iau06xys (iau06arr, fArgs06, xys06table, ttt, ddx, ddy, opt1)
    sethelp;

    rad = 180.0 / pi;

    if (contains(opt1, 'x'))
        [x ,y, s] = iau06xysS (iau06arr, fArgs06, ttt );
    else
        % not sure we need to pass this in? or calc
        jdtt = ttt*36525.0 + 2451545.0;
        [x, y, s] = findxysparam(jdtt, 0.0, 's', xys06table);
    end

    %x =  0.000390429583103574; y =3.52648563069063E-05; s =-1.4673151940162E-08;

    % add corrections if available
    x = x + ddx;
    y = y + ddy;

    % ---------------- now find a
    a = 0.5 + 0.125*(x*x + y*y);    % units take on whatever x and y are

    %if iauhelp == 'x'
    fprintf(1,'06xys  x  %14.12f y  %14.12f s %14.12f a %14.12f rad \n',x, y, s, a );
    fprintf(1,'06xys  x  %14.12f y  %14.12f s %14.12f a %14.12f deg \n',x*rad,y*rad,s*rad,a*rad );
    fprintf(1,'06xys  x  %14.12f" y  %14.12f" s %14.12f" a %14.12fdeg \n',x*rad*3600,y*rad*3600,s*rad*3600,a*rad );
    %end;

    % ----------------- find nutation matrix ----------------------
    nut1(1,1) = 1.0 - a*x*x;
    nut1(1,2) = -a*x*y;
    nut1(1,3) = x;
    nut1(2,1) = -a*x*y;
    nut1(2,2) = 1.0 - a*y*y;
    nut1(2,3) = y;
    nut1(3,1) = -x;
    nut1(3,2) = -y;
    nut1(3,3) = 1.0 - a*(x*x + y*y);
    %nut1

    nut2 = eye(3);
    nut2(1,1) =  cos(s);
    nut2(2,2) =  cos(s);
    nut2(1,2) =  sin(s);
    nut2(2,1) = -sin(s);

    pn = nut1 * nut2;

    %       the matrix appears to be orthogonal now, so the extra processing is not needed.
    %        if (x ~= 0.0) && (y ~= 0.0)
    %            e = atan2(y,x);
    %          else
    %            e = 0.0;
    %          end;
    %        d = atan( sqrt((x^2 + y^2) / (1.0-x^2-y^2)) );
    %        nut1 = rot3mat(-e)*rot2mat(-d)*rot3mat(e+s)

end
