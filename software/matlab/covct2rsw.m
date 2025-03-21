
% ----------------------------------------------------------------------------
%
%                           function covct2rsw
%
%  this function transforms a six by six covariance matrix expressed in cartesian
%    into one expressed in orbit plane, rsw frame
%
%  author        : david vallado                  719-573-2600   20 may 2003
%
%  revisions
%    vallado     - fix indices                                   16 jul 2003
%    vallado     - send out tm                                   21 jul 2003
%
%  inputs          description                                 range / units
%    cartcov     - 6x6 cartesian covariance matrix
%    cartstate   - 6x1 cartesian orbit state                   (x y z vx vy vz)
%
%  outputs       :
%    covoprsw    - 6x6 orbit plane rsw covariance matrix
%
%  locals        :
%    r           - position vector                             m
%    v           - velocity vector                             m/s
%    tm          - transformation matrix
%    temv        - temporary vector
%
%  coupling      :
%    none
%
%  references    :
%    none
%
%  [covoprsw, tm] = covct2rsw( cartcov, cartstate )
% ----------------------------------------------------------------------------

function [covoprsw, tm] = covct2rsw( cartcov, cartstate )
    tm = zeros(6,6);

    x  = cartstate(1);
    y  = cartstate(2);
    z  = cartstate(3);
    vx = cartstate(4);
    vy = cartstate(5);
    vz = cartstate(6);
    r = [x;y;z];
    v = [vx;vy;vz];

    rv = unit(r);
    temv = cross(r,v);
    wv = unit(temv);
    sv = cross(wv,rv);

    for i = 1:6
        for j = 1:6
            tm(i,j) = 0.0;
        end
    end

    tm(1,1) = rv(1);
    tm(1,2) = rv(2);
    tm(1,3) = rv(3);
    tm(2,1) = sv(1);
    tm(2,2) = sv(2);
    tm(2,3) = sv(3);
    tm(3,1) = wv(1);
    tm(3,2) = wv(2);
    tm(3,3) = wv(3);

    tm(4,4) = rv(1);
    tm(4,5) = rv(2);
    tm(4,6) = rv(3);
    tm(5,4) = sv(1);
    tm(5,5) = sv(2);
    tm(5,6) = sv(3);
    tm(6,4) = wv(1);
    tm(6,5) = wv(2);
    tm(6,6) = wv(3);

    covoprsw = tm * cartcov * tm';

end
