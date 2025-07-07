% ------------------------------------------------------------------------------
%
%                           function checkhitearth
%
%  this function checks to see if the trajectory hits the earth during the
%    transfer.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    altPad      - pad for alt above surface                km
%    r1          - initial position vector of int           km
%    v1t         - initial velocity vector of trns          km/s
%    r2          - final position vector of int             km
%    v2t         - final velocity vector of trns            km/s
%    nrev        - number of revolutions                    0, 1, 2,
%
%  outputs       :
%    hitearth    - is earth was impacted                    'y' 'n'
%    hitearthstr - is earth was impacted                    "y - radii" "no"
%    a           - semimajor axis of transfer               km
%
%  locals        :
%    sme         - specific mechanical energy
%    rp          - radius of perigee                        km
%    ecc         - eccentricity of transfer
%    p           - semi-paramater of transfer               km
%    hbar        - angular momentum vector of
%                  transfer orbit
%    radiuspad   - radius including user pad                km
%
%  coupling      :
%    dot         - dot product of vectors
%    mag         - magnitude of a vector
%    MathTimeLibr.cross       - MathTimeLibr.cross product of vectors
%
%  references    :
%    vallado       2022, 483, alg 58
%
% [hitearth, hitearthstr] = checkhitearth ( altpad, r1, v1t, r2, v2t, nrev );
% ------------------------------------------------------------------------------

function [hitearth, hitearthstr] = checkhitearth ( altpad, r1, v1t, r2, v2t, nrev )
    constastro;

    % --------------------------  implementation   -----------------
    show = 'n';

    hitearth = 'n';
    hitearthstr = 'no';

    magr1 = mag(r1);
    magr2 = mag(r2);
    rpad = re + altpad;

    %fprintf(1,'mr1 %11.7f mr2 %11.7f ',magr1, magr2);
    % check whether Lambert transfer trajectory hits the Earth
    % this check may not be needed depending on input data
    if (magr1 < rpad || magr2 < rpad)
        % hitting earth already at start or stop point
        hitearth = 'y';
        hitearthstr = 'y initradii';
        if show == 'y'
            fprintf( 1, 'hitearth? %s \n',hitearthstr );
        end
    else
        rdotv1 = dot(r1, v1t);
        rdotv2 = dot(r2, v2t);

        % Solve for a
        ainv = 2.0 / magr1 - mag(v1t)^2 / mu;

        % Find ecos(E)
        ecosea1 = 1.0 - magr1 * ainv;
        ecosea2 = 1.0 - magr2 * ainv;
        %fprintf(1,'ec1 %11.7f ec2 %11.7f ainv %11.7f ',ecosea1, ecosea2, ainv);

        % Determine radius of perigee
        % 4 distinct cases pass thru perigee
        % nrev > 0 you have to check
        if (nrev > 0)
            a = 1.0 / ainv;
            % elliptical orbit
            if (a > 0.0)
                esinea1 = rdotv1 / sqrt(mu * a);
                ecc = sqrt(ecosea1 * ecosea1 + esinea1 * esinea1);
            else
                % hyperbolic orbit
                esinea1 = rdotv1 / sqrt(mu * abs(-a));
                ecc = sqrt(ecosea1 * ecosea1 - esinea1 * esinea1);
            end
            rp = a * (1.0 - ecc);
            if (rp < rpad)
                hitearth = 'y';
                hitearthstr = 'y Sub_Earth_nrev';
            end
            % nrev = 0, 3 cases:
            % heading to perigee and ending after perigee
            % both headed away from perigee, but end is closer to perigee
            % both headed toward perigee, but start is closer to perigee
        else
            if ((rdotv1 < 0.0 && rdotv2 > 0.0) || (rdotv1 > 0.0 && rdotv2 > 0.0 && ecosea1 < ecosea2) || ...
                    (rdotv1 < 0.0 && rdotv2 < 0.0 && ecosea1 > ecosea2))
                % parabola
                if (abs(ainv) <= 1.0e-10)
                    hbar = cross(r1, v1t);
                    magh = mag(hbar); % find h magnitude
                    rp = magh * magh * 0.5 / mu;
                    if (rp < rpad)
                        hitearth = 'y';
                        hitearthstr = 'y Sub_Earth_parb';
                    end
                else
                    % for both elliptical & hyperbolic
                    a = 1.0 / ainv;
                    esinea1 = rdotv1 / sqrt(mu*abs(a));
                    if (ainv > 0.0)
                        ecc = sqrt(ecosea1 * ecosea1 + esinea1 * esinea1);
                    else
                        ecc = sqrt(ecosea1 * ecosea1 - esinea1 * esinea1);
                    end
                    if (ecc < 1.0)
                        rp = a * (1.0 - ecc);
                        if (rp < rpad)
                            hitearth = 'y';
                            hitearthstr = 'y Sub_Earth_ecc';
                        end
                    else
                        % hyperbolic heading towards the earth
                        if (rdotv1 < 0.0 && rdotv2 > 0.0)
                            rp = a * (1.0 - ecc);
                            if (rp < rpad)
                                hitearth = 'y';
                                hitearthstr = 'y Sub_Earth_hyp';
                            end
                        end
                    end
                    %   fprintf( 1, 'hitearth? %s rp %11.7f  %11.7f km \n',hitearthstr, rp*6378.137, rpad*6378.137 );
                end % nrev = 0 checks
            end % end of perigee check

            if show == 'y'
                fprintf( 1, 'hitearth? %s rp %11.7f km \n',hitearth, rp*6378.0 );
            end
        end
    end % checkhitearth

end