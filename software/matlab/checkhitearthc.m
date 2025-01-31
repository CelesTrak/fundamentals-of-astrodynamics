% ------------------------------------------------------------------------------
%
%                           function checkhitearthc
%
%  this function checks to see if the trajectory hits the earth during the
%    transfer. Calc in canonical units.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    altPadc     - pad for alt above surface                 er
%    r1c         - initial position vector of int            er
%    v1tc        - initial velocity vector of trns           er/tu
%    r2c         - final position vector of int              er
%    v2tc        - final velocity vector of trns             er/tu
%    nrev        - number of revolutions                     0, 1, 2,
%
%  outputs       :
%    hitearth    - is earth was impacted                     'y' 'n'
%    hitearthstr - is earth was impacted                     "y - radii" "no"
%    a           - semimajor axis of transfer                er
%
%  locals        :
%    sme         - specific mechanical energy
%    rp          - radius of perigee                         er
%    ecc         - eccentricity of transfer
%    p           - semi-paramater of transfer                er
%    hbar        - angular momentum vector of
%                  transfer orbit
%    radiuspadc  - radius including user pad                 er
%
%  coupling      :
%    dot         - dot product of vectors
%    mag         - magnitude of a vector
%    cross       - cross product of vectors
%
%  references    :
%    vallado       2022, 483, alg 58
%
% [hitearth, hitearthstr] = checkhitearthc ( altpadc, r1c, v1tc, r2c, v2tc, nrev );
% ------------------------------------------------------------------------------

function [hitearth, hitearthstr] = checkhitearthc ( altpadc, r1c, v1tc, r2c, v2tc, nrev )

    % --------------------------  implementation   -----------------
    show = 'n';

    hitearth = 'n';
    hitearthstr = 'no';

    magr1c = mag(r1c);
    magr2c = mag(r2c);
    rpadc = 1.0 + altpadc;

    %fprintf(1,'mr1 %11.7f mr2 %11.7f ',magr1c, magr2c);
    % check whether Lambert transfer trajectory hits the Earth
    % this check may not be needed depending on input data
    if (magr1c < rpadc || magr2c < rpadc)
        % hitting earth already at start or stop point
        hitearth = 'y';
        hitearthstr = strcat(hitearth,' initradii');
        if show == 'y'
            fprintf( 1, 'hitearth? %s \n',hitearthstr );
        end
    else
        rdotv1c = dot(r1c, v1tc);
        rdotv2c = dot(r2c, v2tc);

        % Solve for a
        ainv = 2.0 / magr1c - mag(v1tc)^2;

        % Find ecos(E)
        ecosea1 = 1.0 - magr1c * ainv;
        ecosea2 = 1.0 - magr2c * ainv;
        %fprintf(1,'ec1 %11.7f ec2 %11.7f ainv %11.7f ',ecosea1, ecosea2, ainv);

        % Determine radius of perigee
        % 4 distinct cases pass thru perigee
        % nrev > 0 you have to check
        if (nrev > 0)
            a = 1.0 / ainv;
            % elliptical orbit
            if (a > 0.0)
                esinea1 = rdotv1c / sqrt(a);
                ecc = sqrt(ecosea1 * ecosea1 + esinea1 * esinea1);
            else
                % hyperbolic orbit
                esinea1 = rdotv1c / sqrt(abs(-a));
                ecc = sqrt(ecosea1 * ecosea1 - esinea1 * esinea1);
            end
            rp = a * (1.0 - ecc);
            if (rp < rpadc)
                hitearth = 'y';
                hitearthstr = strcat(hitearth,' Sub_Earth_nrev');
            end
            % nrev = 0, 3 cases:
            % heading to perigee and ending after perigee
            % both headed away from perigee, but end is closer to perigee
            % both headed toward perigee, but start is closer to perigee
        else
            if ((rdotv1c < 0.0 && rdotv2c > 0.0) || (rdotv1c > 0.0 && rdotv2c > 0.0 && ecosea1 < ecosea2) || ...
                    (rdotv1c < 0.0 && rdotv2c < 0.0 && ecosea1 > ecosea2))
                % parabola
                if (abs(ainv) <= 1.0e-10)
                    hbar = cross(r1c, v1tc);
                    magh = mag(hbar); % find h magnitude
                    rp = magh * magh * 0.5 / mu;
                    if (rp < rpadc)
                        hitearth = 'y';
                        hitearthstr = strcat(hitearth, ' Sub_Earth_para');
                    end
                else
                    % for both elliptical & hyperbolic
                    a = 1.0 / ainv;
                    esinea1 = rdotv1c / sqrt(abs(a));
                    if (ainv > 0.0)
                        ecc = sqrt(ecosea1 * ecosea1 + esinea1 * esinea1);
                    else
                        ecc = sqrt(ecosea1 * ecosea1 - esinea1 * esinea1);
                    end
                    if (ecc < 1.0)
                        rp = a * (1.0 - ecc);
                        if (rp < rpadc)
                            hitearth = 'y';
                            hitearthstr = strcat(hitearth, ' Sub_Earth_ell');
                        end
                    else
                        % hyperbolic heading towards the earth
                        if (rdotv1c < 0.0 && rdotv2c > 0.0)
                            rp = a * (1.0 - ecc);
                            if (rp < rpadc)
                                hitearth = 'y';
                                hitearthstr = strcat(hitearth, ' Sub_Earth_hyp');
                            end
                        end
                    end
                    %   fprintf( 1, 'hitearth? %s rp %11.7f  %11.7f km \n',hitearthstr, rp*6378.137, rpad*6378.137 );
                end % nrev = 0 checks
            end
        end  % end of nrev = 0 checks
        if show == 'y'
            fprintf( 1, 'hitearth? %s rp %11.7f km \n',hitearth, rp*6378.0 );
        end
    end % checkhitearthc

end