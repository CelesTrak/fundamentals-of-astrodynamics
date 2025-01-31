% ------------------------------------------------------------------------------
%
%                           function rv2sez
%
%  this function converts position and velocity vectors into site topcentric 
%    coordinates. this supports the hill transformations. the reverse
%    values are found using the transmat transpose. 
%
%  author        : david vallado                  719-573-2600  13 aug 2010
%
%  revisions
%                -
%  inputs          description                    range / units
%    reci        - position vector                km
%    veci        - velocity vector                km/s
%    lat         - rel lat of 2nd satellite wrt 1st
%    lon         - rel lon of 2nd satellite wrt 1st
%  outputs       :
%    rsez        - position vector                km
%    vsez        - velocity vector                km/s
%
%  locals        :
%    temp        - temporary position vector
%
%  coupling      :
%    mag         - magnitude of a vector
%
%  references    :
%    vallado       2004, 162
%
% [rsez,vsez,transmat] = rv2sez( reci,veci, lat, lon );
% ------------------------------------------------------------------------------

function [rsez,vsez,transmat] = rv2sez( reci,veci, lat, lon )

        % each of the components must be unit vectors
        % zenith component
        rs(1) = mag(reci)*cos(lat)*cos(lon);
        rs(2) = mag(reci)*cos(lat)*sin(lon);
        rs(3) = mag(reci)*sin(lat);
        zvec = unit(rs);
         
        % east component
        k = [0 0 1];
        evec    = cross(k,zvec);
        evec    = unit( evec );

        % south component
        svec    = cross(evec,zvec);
        svec    = unit( svec );

        % assemble transformation matrix from to rsw frame (individual
        % components arranged in row vectors)
        transmat(1,1) = svec(1);
        transmat(1,2) = svec(2);
        transmat(1,3) = svec(3);
        transmat(2,1) = evec(1);
        transmat(2,2) = evec(2);
        transmat(2,3) = evec(3);
        transmat(3,1) = zvec(1);
        transmat(3,2) = zvec(2);
        transmat(3,3) = zvec(3);

        rsez = transmat*reci';
        vsez = transmat*veci';

%   alt approach
%       rrsw(1) = mag(reci);
%       rrsw(2) = 0.0;
%       rrsw(3) = 0.0;
%       vrsw(1) = dot(reci,veci)/rrsw(1);
%       vrsw(2) = sqrt(veci(1)^2 + veci(2)^2 + veci(3)^2 - vrsw(1)^2);
%       vrsw(3) = 0.0;

