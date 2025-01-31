% -----------------------------------------------------------------------------
%
%                           function gstime
%
%  this function finds the greenwich sidereal time (iau-82).
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    jdut1       - julian date in ut1                       days from 4713 bc
%
%  outputs       :
%    gstime      - greenwich sidereal time                  0 to 2pi rad
%
%  locals        :
%    temp        - temporary variable for doubles           rad
%    tut1        - julian centuries from the
%                  jan 1, 2000 12 h epoch (ut1)
%
%  coupling      :
%    none
%
%  references    :
%    vallado       2022, 189, eq 3-48
%
% gst = gstime(jdut1);
% -----------------------------------------------------------------------------

function gst = gstime(jdut1)

    twopi      = 2.0*pi;
    deg2rad    = pi/180.0;

    % ------------------------  implementation   ------------------
    tut1= ( jdut1 - 2451545.0 ) / 36525.0;

    temp = - 6.2e-6 * tut1 * tut1 * tut1 + 0.093104 * tut1 * tut1  ...
        + (876600.0 * 3600.0 + 8640184.812866) * tut1 + 67310.54841;

    % 360/86400 = 1/240, to deg, to rad
    temp = rem( temp*deg2rad/240.0,twopi );

    % ------------------------ check quadrants --------------------
    if ( temp < 0.0 )
        temp = temp + twopi;
    end

    gst = temp;

end