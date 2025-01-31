% -----------------------------------------------------------------------------
%                           procedure lstime
%
%  this procedure finds the local sidereal time at a given location. gst is from iau-82.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    lon         - site longitude (west -)                  -2pi to 2pi rad
%    jdut1       - julian date in ut1                       days from 4713 bc
%
%  outputs       :
%    lst         - local sidereal time                      0.0 to 2pi rad
%    gst         - greenwich sidereal time                  0.0 to 2pi rad
%
%  locals        :
%    none.
%
%  coupling      :
%    gstime        finds the greenwich sidereal time
%
%  references    :
%    vallado       2022, 190, Alg 15, ex 3-5
%
% [lst,gst] = lstime ( lon, jd );
% -----------------------------------------------------------------------------

function [lst,gst] = lstime ( lon, jd )

    twopi  = 2.0*pi;

    % ------------------------  implementation   ------------------
    [gst] = gstime( jd );
    lst = lon + gst;

    % ----------------------- check quadrants ---------------------
    lst = rem( lst,twopi );
    if ( lst < 0.0 )
        lst= lst + twopi;
    end

end