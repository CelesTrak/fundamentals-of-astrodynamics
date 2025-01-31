% -----------------------------------------------------------------------------
%
%                           function rv2ntw
%
%  this function converts position and velocity vectors into normal, in-
%    track, and cross-track coordinates.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    r           - position vector                           km
%    v           - velocity vector                           km/s
%
%  outputs       :
%    rntw        - position vector                           km
%    vntw        - velocity vector                           km/s
%    transmat    - transformation matrix
%
%  locals        :
%    tempvec     - temporary vector
%    nvec,tvec,wvec - direction Cosines
%
%  coupling      :
%
%  references    :
%    vallado       2022, 166
%
% [rntw,vntw,transmat] = rv2ntw( r,v );
% ------------------------------------------------------------------------------

function [rntw, vntw, transmat] = rv2ntw(r, v)
    % compute satellite velocity vector magnitude
    vmag = mag(v);

    % in order to work correctly each of the components must be
    %  unit vectors !
    % in-velocity component
    tvec = v / vmag;

    % cross-track component
    wvec = cross(r, v);
    wvec = unit( wvec );

    % along-radial component
    nvec = cross(tvec, wvec);
    nvec = unit( nvec );

    % assemble transformation matrix from to ntw frame (individual
    %  components arranged in row vectors)
    transmat(1,1) = nvec(1);
    transmat(1,2) = nvec(2);
    transmat(1,3) = nvec(3);
    transmat(2,1) = tvec(1);
    transmat(2,2) = tvec(2);
    transmat(2,3) = tvec(3);
    transmat(3,1) = wvec(1);
    transmat(3,2) = wvec(2);
    transmat(3,3) = wvec(3);

    rntw = matvecmult(transmat,r,3);
    vntw = matvecmult(transmat,v,3);

end