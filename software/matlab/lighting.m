% ------------------------------------------------------------------------------
%
%                           function lighting
%
%  this function determines if a spacecraft is sunlit or in the dark at a
%    particular time.  an oblate earth and cylindrical shadow is assumed.
%
%  author        : david vallado                  719-573-2600   27 may 2002
%
%  revisions
%                -
%
%  inputs          description                    range / units
%    r           - position vector of sat         er
%    jd          - julian date at desired time    days from 4713 bc
%    whichkind   - spherical or ellipsoidal earth 's', 'e'*default
%
%  outputs       :
%    vis         - visibility flag                'yes','no '
%
%  locals        :
%    rtasc       - suns right ascension           rad
%    decl        - suns declination               rad
%    rsun        - sun vector                     au
%    auer        - conversion from au to er
%
%  coupling      :
%    sun         - position vector of sun
%    lncom1      - multiple a vector by a constant
%    sight       - does line-of-sight exist beteen vectors
%
%  references    :
%    vallado       2001, 291-295, alg 35, ex 5-6
%
% [lit] = lighting ( r, jd, whichkind );
% ------------------------------------------------------------------------------

function [lit] = lighting ( r, jd, whichkind )

        auer= 149597870.0 /6378.1363;

        % -------------------------  implementation   -------------------------
        [rsun,~,~] = sun( jd );
        rsun= auer*rsun;

        % ------------ is the satellite in the shadow? ----------------
        [lit] = sight( rsun,r,whichkind );

