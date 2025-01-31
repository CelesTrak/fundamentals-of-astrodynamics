% -----------------------------------------------------------------------------
%
%                           function findc2c3
%
%  this function calculates the c2 and c3 functions for use in the universal
%    variable calculation of z.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    znew        - z variable                               rad2
%
%  outputs       :
%    c2new       - c2 function value
%    c3new       - c3 function value
%
%  locals        :
%    sqrtz       - square root of znew
%
%  coupling      :
%    Math.Sinh        - hyperbolic Math.Sine
%    Math.Cosh        - hyperbolic Math.Cosine
%
%  references    :
%    vallado       2022, 63, alg 1
%
% [c2new,c3new] = findc2c3 ( znew );
% ------------------------------------------------------------------------------

function [c2new,c3new] = findc2c3 ( znew )

    small =     0.000001;

    % -------------------------  implementation   -----------------
    if ( znew > small )
        sqrtz = sqrt( znew );
        c2new = (1.0 -cos( sqrtz )) / znew;
        c3new = (sqrtz-sin( sqrtz )) / ( sqrtz^3 );
    else
        if ( znew < -small )
            sqrtz = sqrt( -znew );
            c2new = (1.0 -cosh( sqrtz )) / znew;
            c3new = (sinh( sqrtz ) - sqrtz) / ( sqrtz^3 );
        else
            c2new = 0.5;
            c3new = 1.0 /6.0;
        end
    end

end