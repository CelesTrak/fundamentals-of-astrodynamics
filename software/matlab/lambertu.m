%  ------------------------------------------------------------------------------
%
%                           procedure lambertuniv
%
%  this procedure solves the lambert problem for orbit determination and returns
%    the velocity vectors at each of two given position vectors.  the solution
%    uses universal variables for calculation and a bissection technique for
%    updating psi. note that v1 (the initial orbit velocity and not the transfer 
%    orbit velocity at the start) is not formally part of the lambert solution, 
%    but here, it's needed for the hodegraph implementation for transfers very 
%    near 180 deg. it can be set to zero if you're not near the 180 deg transfer. 
%
%  algorithm     : setting the initial bounds:
%                  using -8pi and 4pi2 will allow single rev solutions
%                  using -4pi2 and 8pi2 will allow multi-rev solutions
%                  the farther apart the initial guess, the more iterations
%                    because of the iteration
%                  inner loop is for special cases. must be sure to exit both!
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    r1          - ijk position vector 1                     km
%    r2          - ijk position vector 2                     km
%    v1          - ijk velocity vector 1 if avail            km/s
%    dm          - direction of motion                       'L', 'S'
%    de          - orbital energy                            'L', 'H'
%                  only affects nrev >= 1 upper/lower bounds
%    dtsec       - time between r1 and r2                    sec
%    nrev        - number of revs to complete                0, 1, 2, 3,
%    kbi         - time value for min
%    show        - control output don't output for speed      'y', 'n'
%
%  outputs       :
%    v1t         - ijk transfer velocity vector              km/s
%    v2t         - ijk transfer velocity vector              km/s
%    detailAll   - detail all output 
%
%  locals        :
%    vara        - variable of the iteration,
%                  not the semi or axis!
%    y           - area between position vectors
%    upper       - upper bound for z
%    lower       - lower bound for z
%    cosdeltanu  - cosine of true anomaly change             rad
%    f           - f expression
%    g           - g expression
%    gdot        - g dot expression
%    xold        - old universal variable x
%    xoldcubed   - xold cubed
%    zold        - old value of z
%    znew        - new value of z
%    c2new       - c2(z) function
%    c3new       - c3(z) function
%    timenew     - new time                                  sec
%    small       - tolerance for roundoff errors
%    i, j        - index
%
%  coupling
%    mag         - magnitude of a vector
%    dot         - dot product of two vectors
%    findc2c3    - find c2 and c3 functions
%
%  references    :
%    vallado       2022, 499, alg 60, ex 7-5
%
% [v1tu, v2tu, detailAll] = lambertu(r1, r2, v1, dm, de, nrev, dtsec, kbi, show );
% ------------------------------------------------------------------------------

function [v1tu, v2tu, detailAll] = lambertu(r1, r2, v1, dm, de, nrev, dtsec, kbi, show )

    % -------------------------  implementation   -------------------------
    constastro;
    small = 0.000001; % can affect cases where znew is multiples of 2pi^2
    numiter= 40;
    detailAll = 'lambertu\n';
    for i= 1 : 3
        v1tu(i) = 0.0;
        v2tu(i) = 0.0;
    end

    dtold = 0.0;

    % try canonical units for testing
    %constastro;
    %mu = 1.0;
    %r1 = r1/re;
    %r2 = r2/re;
    %dtsec = dtsec / tusec;

    % ---- find parameters that are constant for the initial geometry
    magr1 = mag(r1);
    magr2 = mag(r2);

    % this value stays constant in all calcs, vara changes with df
    cosdeltanu = dot(r1,r2) / (magr1*magr2);
    if abs(cosdeltanu) > 1.0
        cosdeltanu = sign(cosdeltanu) * 1.0;
    end
    if ( dm == 'L' )  %dm == 'l'
        vara = -sqrt( magr1*magr2*(1.0 + cosdeltanu) );
    else
        vara =  sqrt( magr1*magr2*(1.0 + cosdeltanu) );
    end

    % setup variables for speed
    oomu = 1.0 / sqrt(mu);

    % --------- set up initial bounds for the bissection ----------
    if ( nrev == 0 )
        lower = -16.0*pi*pi;  % allow hyperbolic and parabolic solutions
        upper =  4.0*pi*pi;  % could be negative infinity for all cases
    else
        % set absolute limits for multi-rev cases
        lower = 4.0*nrev^2*pi*pi;
        upper = 4.0*(nrev + 1.0)^2*pi*pi;
        % adjust based on long or short way if dm == 'l'
        %if ((dm == 'l') && (df == 'd')) || ((dm == 's') && (df == 'r'))
        %if ((df == 'r') && (dm == 's')) || ((df == 'd') && (dm == 'l'))
        if (de == 'H') %  && (dm == 'L')) || ((de == 'L') && (dm == 'L'))
            upper = kbi;
        else
            lower = kbi;
        end
    end

    % ---------------  form initial guesses   ---------------------
    dtdpsi = 0.0;
    xold   = 0.0;
    psinew = 0.0;
    if (nrev == 0)
        % use log to get initial guess
        % empirical relation here from 10000 random draws
        % 10000 cases up to 85000 dtsec  0.11604050x + 9.69546575
        psiold = (log(dtsec) - 9.61202327)/0.10918231;
        if psiold > upper
            psiold = upper - pi;
        end
        % try arora-russel 2010 approach tbi( 1 = psi, 2 = tof in sec)
        %         k1 = 1.0;
        %         k2 = 5.712388980384687;
        %         yx = magr1 + magr2 - k1*vara;
        %         taux = sqrt(yx)/sqrt(mu) * (k2*yx + vara);
    else
        %if (de == 'L')  % dm == 's' seemed better other way?????????????
        psiold = lower + (upper - lower)*0.5;
        % else
        %            psiold = lower + (upper - lower)*0.6;
        %         end
        % try arora-russel 2010 approach tbi( 1 = psi, 2 = tof in sec)
        %         if dm == 'S'
        %             psiaux = 0.3 * tbi(nrev,1) + 0.74 * nrev^2 * pi^2;
        %         else
        %             psiaux = 0.3 * tbi(nrev,1) + 0.74 * (nrev + 1)^2 * pi^2;
        %         end
        %         [c2new,c3new] = findc2c3( psiaux );
        %         y= magr1 + magr2 - ( vara*(1.0-psiaux*c3new) / sqrt(c2new) );
        %         taux    = (xold^3*c3new + vara*sqrt(y)) * oomu;
        %         a = (log(taux) - log(tbi(nrev,2))) / (psiaux - tbi(nrev,1));
        %         b = -2.0 * a * tbi(nrev,1);
        %         c = log(tbi(nrev,2)) + a * tbi(nrev,1)^2;
        %         a = a /log(dtsec);
        %         b = b /log(dtsec);
        %         c = c /log(dtsec);
        %         discrim = b*b - 4.0 *a*c;
        %         if ( discrim > 0.0  )
        %             psiold1 = ( -b + sqrt(discrim) ) / ( 2.0 *a );
        %                psiold = ( -b - sqrt(discrim) ) / ( 2.0 *a );
        %psiold = lower + (upper - lower)*0.3;
    end

    [c2new,c3new] = findc2c3( psiold );

    oosqrtmu = 1.0 / sqrt(mu);

    % find initial dtold from psiold
    if (abs(c2new) > small)
        y = magr1 + magr2 - (vara * (1.0 - psiold * c3new) / sqrt(c2new));
    else
        y = magr1 + magr2;
    end
    if (abs(c2new) > small)
        xold = sqrt(y / c2new);
    else
        xold = 0.0;
    end
    xoldcubed = xold * xold * xold;
    dtold = (xoldcubed * c3new + vara * sqrt(y)) * oosqrtmu;

    % -------  determine if  the orbit is possible at all ---------
    if ( abs( vara ) > 0.2 )   % not exactly zero
        loops  = 0;
        ynegktr= 1;  % y neg ktr
        dtnew = -10.0;

        while ((abs(dtnew-dtsec) >= small) && (loops < numiter) && (ynegktr <= 10))
            % fprintf(1,'%3i  dtnew-dtsec %11.7f yneg %3i \n',loops,dtnew-dtsec,ynegktr );
            if ( abs(c2new) > small )
                y= magr1 + magr2 - ( vara*(1.0 - psiold*c3new)/sqrt(c2new) );
            else
                y= magr1 + magr2;
            end
            % ----------- check for negative values of y ----------
            if ( (vara > 0.0) && (y < 0.0) )  % ( vara > 0.0 ) &
                ynegktr= 1;
                while (( y < 0.0 ) && ( ynegktr < 10 ))
                    psinew = 0.8*(1.0 / c3new)*( 1.0 - (magr1 + magr2)*sqrt(c2new)/vara  );
                    % -------- find c2 and c3 functions -----------
                    [c2new,c3new] = findc2c3( psinew );
                    psiold = psinew;
                    lower  = psiold;
                    if ( abs(c2new) > small )
                        y= magr1 + magr2 - ( vara*(1.0-psiold*c3new) / sqrt(c2new) );
                    else
                        y= magr1 + magr2;
                    end
                    % outfile
                    if show == 'y'
                        detailAll = sprintf('%s %s',detailAll, sprintf('yneg %3i  y %11.7f lower %11.7f c2 %11.7f psinew %11.7f yneg %3i \n', ...
                            loops,y,lower,c2new,psinew,ynegktr) );
                    end
                    ynegktr = ynegktr + 1;
                end % while
            end  % if  y neg

            loops = loops + 1;

            if ( ynegktr < 10 )
                if ( abs(c2new) > small )
                    xold= sqrt( y / c2new );
                else
                    xold= 0.0;
                end
                xcubed= xold^3;
                dtnew    = (xcubed*c3new + vara*sqrt(y)) * oomu;

                % try newton rhapson iteration to update psi
                if abs(psiold) > 1e-5
                    c2dot = 0.5/psiold * (1.0 - psiold*c3new - 2.0*c2new);
                    c3dot = 0.5/psiold * (c2new - 3.0*c3new);
                else  % case for parabolic orbit
                    c2dot = -1.0/factorial(4) + 2.0*psiold/factorial(6) - 3.0*psiold^2/factorial(8) + 4.0*psiold^3/factorial(10) - 5.0*psiold^4/factorial(12);
                    c3dot = -1.0/factorial(5) + 2.0*psiold/factorial(7) - 3.0*psiold^2/factorial(9) + 4.0*psiold^3/factorial(11) - 5.0*psiold^4/factorial(13);
                end
                dtdpsi = (xcubed*(c3dot - 3.0*c3new*c2dot/(2.0*c2new)) + 0.125*vara * (3.0*c3new*sqrt(y)/c2new + vara/xold)) * oomu;
                % Newton iteration test to see if it keeps within the bounds
                psinew =  psiold - (dtnew - dtsec)/dtdpsi;

                % check if newton guess for psi is outside bounds (too steep a slope)
                if psinew > upper || psinew < lower
                    % --------  readjust upper and lower bounds -------
                    if (de == 'L' || (nrev == 0))
                        if (dtold < dtsec)
                            lower = psiold;
                        else
                            upper = psiold;
                        end
                    else
                        if (dtold < dtsec)
                            upper = psiold;
                        else
                            lower = psiold;
                        end
                    end
                    psinew= (upper+lower) * 0.5;
                    psilast = psinew;
                end

                % ------------- find c2 and c3 functions ----------
                [c2new,c3new] = findc2c3( psinew );
                %if nrev > 0
                    if (show == 'y')
                        detailAll = sprintf('%s %s', detailAll, sprintf('%3i  y %11.7f x %11.7f %11.7f dtnew %11.7f %11.7f %11.7f psinew %11.7f %11.7f \n', ...
                            loops,y,xold,dtsec, dtnew, lower, upper, psinew, dtdpsi) ); %(dtnew - dtsec)/dtdpsi );  % c2dot, c3dot
                    end
                %end
                psilast = psiold;  % keep previous iteration
                psiold = psinew;
                dtold = dtnew;

                % --- make sure the first guess isn't too close ---
                if ( (abs(dtnew - dtsec) < small) && (loops == 1) )
                    dtnew= dtsec - 1.0;
                end
            end  % if  ynegktr < 10

            %              fprintf(1,'#%3i  y %11.7f x %11.7f dtnew %11.7f psinew %11.7f \n',loops,y,x,dtnew,psinew );
            %              fprintf(1,'%3i  y %11.7f x %11.7f dtnew %11.7f psinew %11.7f \n',loops,y/re,x/sqrt(re),dtnew/tusec,psinew );
            %              fprintf(1,'%3i  y %11.7f x %11.7f dtnew %11.7f psinew %11.7f \n',loops,y/re,x/sqrt(re),dtnew/60.0,psinew );
            ynegktr = 1;
        end % while loop

        if ( (loops >= numiter) || (ynegktr >= 10) )
            detailAll= sprintf('%s %s',detailAll, 'gnotconv dtsec dtnew delta (s) ', dtsec, ' ', num2str(abs(dtnew - dtsec)));
            if ( ynegktr >= 10 )
                detailAll= sprintf('%s %s',detailAll, 'y negati');
            end
        else
            % --- use f and g series to find velocity vectors -----
            f   = 1.0 - y/magr1;
            gdot= 1.0 - y/magr2;
            g   = 1.0 / (vara*sqrt( y/mu ));  % 1 over g
            %fprintf(1,'%11.7f  %11.7f  %11.7f \n',f, gdot, g);

            %  fdot = sqrt(mu*y)*(-magr2-magr1 + y)/(magr1*magr2*vara);
            %  f*gdot - fdot*g
            for i= 1 : 3
                v1tu(i) = ( r2(i) - f*r1(i) )*g;
                v2tu(i) = ( gdot*r2(i) - r1(i) )*g;
            end
        end   % if  the answer has converged
    else
        detailAll= sprintf('%s %s',detailAll, 'impos180');

        %call Battin...

        % do hohmann but in 3-d...
        % can't do bissection because w series is not accurate
        mum = 3.986004418e5;   % 14 m3/s2
        atx = (mum*(dtsec/(1.0*pi))^2)^(1.0/3.0);  % 1pi since half period
        v1tmag = sqrt(2.0*mum/(magr1) - mum/atx);
        v2tmag = sqrt(2.0*mum/(magr2) - mum/atx);
        wx = cross(r1,v1);
        wxu = unit(wx);
        v1dir = cross(r1,wxu); % get retro direc
        v2dir = cross(r2,wxu); % get retro direc
        v1diru = unit(v1dir);
        v2diru = unit(v2dir);
        v1t = -v1tmag*v1diru;
        v2t = -v2tmag*v2diru;
        fprintf(1,'lamb univ imposs r1  %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f \n',r1, v1t);
        fprintf(1,'lamb univ imposs r2  %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f \n',r2, v2t);

        v1tu = v1t;% / velkmps;
        v2tu = v2t;% / velkmps;

        tof = dtsec;
        %pause;
        %             % use JGCD 2011 v34 n6 1925 to solve 180 deg case
        %             p = 2.0*magr1*magr2 / (magr1 + magr2);
        %             ecc = sqrt(1.0 - 4.0*magr1*magr2 / ((magr1+magr2)^2) );
        %             dt = sqrt(pi^2/mu * (p / (1.0 - ecc^2))^3);
        %             dnu = acos(cosdeltanu);
        %      fprintf(1,'hodo %11.6f   %11.6f  %11.6f  %11.6f %14.10f \n',p, ecc, dt, dtsec, vara);
        %             if abs(dt-dtsec) < 160
        %                 [v1dv, v2dv] = lambhodograph( r1, v1, r2, p, ecc, dnu, dtsec );
        %             end;
    end  % if  var a > 0.0

    if show == 'y'
        if (strcmp(detailAll, '      ok') ~= 0)
            [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe (r1,v1tu);
            detailAll= sprintf('%s %s',detailAll, sprintf('%10s %3i %3i %2s %2s %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f case %11.7f %11.7f %11.7f %11.7f %11.7f ', ...
                detailAll, loops, nrev, dm, de, dtnew, y, xold, v1tu(1), v1tu(2), v1tu(3), v2tu(1), v2tu(2), v2tu(3), lower, upper, psinew, dtdpsi, ecc) ); %(dtnew - dtsec)/dtdpsi, ecc );  % c2dot, c3dot
        end
    end
end

