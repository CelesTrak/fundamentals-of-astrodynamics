% ----------------------------------------------------------------------------
%
%                           procedure newtone
%
%  this procedure finds the mean and true anomaly given the eccentric anomaly.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    ecc         - eccentricity                   0.0 to
%    eccanom     - eccentric anomaly              0.0 to 2pi rad
%
%  outputs       :
%    m           - mean anomaly                   0.0 to 2pi rad
%    nu          - true anomaly                   0.0 to 2pi rad
%
%  locals        :
%    e1          - eccentric anomaly, next value  rad
%    sinv        - sine of nu
%    cosv        - cosine of nu
%
%  coupling      :
%
%  references    :
%    vallado       2022, 65, alg 2, ex 2-1
%
% [m,nu] = newtone ( ecc,e0 );
% ------------------------------------------------------------------------------

function [m,nu] = newtone ( ecc,e0 )

    % -------------------------  implementation   -----------------
    small =     0.00000001;

    % ------------------------- circular --------------------------
    if ( abs( ecc ) < small )
        m = e0;
        nu= e0;
    else

        % ----------------------- elliptical ----------------------
        if ( ecc < 0.999  )
            m= e0 - ecc*sin(e0);
            sinv= ( sqrt( 1.0 -ecc*ecc ) * sin(e0) ) / ( 1.0 -ecc*cos(e0) );
            cosv= ( cos(e0)-ecc ) / ( 1.0  - ecc*cos(e0) );
            nu  = atan2( sinv,cosv );
        else

            % ---------------------- hyperbolic  ------------------
            if ( ecc > 1.0001  )
                m= ecc*sinh(e0) - e0;
                sinv= ( sqrt( ecc*ecc-1.0  ) * sinh(e0) ) / ( 1.0  - ecc*cosh(e0) );
                cosv= ( cosh(e0)-ecc ) / ( 1.0  - ecc*cosh(e0) );
                nu  = atan2( sinv,cosv );
            else

                % -------------------- parabolic ------------------
                m= e0 + (1.0 /3.0 )*e0*e0*e0;
                nu= 2.0 *datan(e0);
            end
        end
    end

end