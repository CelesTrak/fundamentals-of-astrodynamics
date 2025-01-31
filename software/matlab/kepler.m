% ------------------------------------------------------------------------------
%
%                           function kepler
%
%  this function solves keplers problem for orbit determination and returns a
%    future geocentric equatorial (ijk) position and velocity vector.  the
%    solution uses universal variables.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    r1          - ijk position vector - initial            km
%    vo          - ijk velocity vector - initial            km / s
%    dtsec       - length of time to propagate              s
%
%  outputs       :
%    r           - ijk position vector                      km
%    v           - ijk velocity vector                      km / s
%    error       - error flag                               'ok',
%
%  locals        :
%    f           - f expression
%    g           - g expression
%    fdot        - f dot expression
%    gdot        - g dot expression
%    xold        - old universal variable x
%    xoldsqrd    - xold squared
%    xnew        - new universal variable x
%    xnewsqrd    - xnew squared
%    znew        - new value of z
%    c2new       - c2(psi) function
%    c3new       - c3(psi) function
%    dtsec       - change in time                           s
%    timenew     - new time                                 s
%    rdotv       - result of r1 dot vo
%    a           - semi or axis                             km
%    alpha       - reciprocol  1/a
%    sme         - specific mech energy                     km2 / s2
%    period      - time period for satellite                s
%    s           - variable for parabolic case
%    w           - variable for parabolic case
%    h           - angular momentum vector
%    temp        - temporary real*8 value
%    i           - index
%
%  coupling      :
%    mag         - magnitude of a vector
%    findc2c3    - find c2 and c3 functions
%
%  references    :
%    vallado       2022, 94, alg 8, ex 2-4
%
% [r, v] =  kepler(ro, vo, dtseco);
% ------------------------------------------------------------------------------

function [r, v] =  kepler(ro, vo, dtseco)

    % -------------------------  implementation   -----------------
    % set constants and intermediate printouts
    constmath;
    constastro;
    show = 'n';
    numiter    =    50;

    if show =='y'
        fprintf(1,' ro %16.8f %16.8f %16.8f km \n',ro );
        fprintf(1,' vo %16.8f %16.8f %16.8f km/s \n',vo );
        fprintf(1,' ro %16.8f %16.8f %16.8f ER \n',ro/re );
        fprintf(1,' vo %16.8f %16.8f %16.8f ER/TU \n',vo/velkmps );
    end

    % --------------------  initialize values   -------------------
    ktr   = 0;
    xold  = 0.0;
    znew  = 0.0;
    errork = '      ok';
    dtsec = dtseco;
    nrev = 0;
    smu = sqrt(mu);

    if ( abs( dtseco ) > small )
        magro = mag( ro );
        magvo = mag( vo );
        rdotv= dot( ro,vo );

        % -------------  find sme, alpha, and a  ------------------
        sme= ( (magvo^2)*0.5  ) - ( mu /magro );
        alpha= -sme*2.0/mu;

        if ( abs( sme ) > small )
            a= -mu / ( 2.0 *sme );
        else
            a= infinite;
        end
        if ( abs( alpha ) < small )   % parabola
            alpha= 0.0;
        end

        if show =='y'
            fprintf(1,' sme %16.8f  a %16.8f alp  %16.8f ER \n',sme/(mu/re), a/re, alpha*re );
            fprintf(1,' sme %16.8f  a %16.8f alp  %16.8f km \n',sme, a, alpha );
            fprintf(1,' ktr      xn        psi           r          xn+1        dtn \n' );
        end

        % ------------   setup initial guess for x  ---------------
        % -----------------  circle and ellipse -------------------
        if ( alpha >= small )
            period= twopi * sqrt( abs(a)^3.0/mu  );
            % ------- next if needed for 2body multi-rev ----------
            if ( abs( dtseco ) > abs( period ) )
                % including the truncation will produce vertical lines that are parallel
                % (plotting chi vs time)
                dtsec = rem( dtseco, period );
                nrev = floor(dtseco/period);
            end
            xold = smu *dtsec * alpha;
        else
            % --------------------  parabola  ---------------------
            if ( abs( alpha ) < small )
                h = cross( ro,vo );
                magh = mag(h);
                p= magh*magh/mu;
                s= 0.5  * (halfpi - atan( 3.0 *sqrt( mu / (p*p*p) )* dtsec ) );
                w= atan( tan( s )^(1.0 /3.0 ) );
                xold = sqrt(p) * ( 2.0 *cot(2.0 *w) );
                alpha= 0.0;
            else
                % ------------------  hyperbola  ------------------
                temp= -2.0 * mu * dtsec / ...
                    ( a*( rdotv + sign(dtsec)*sqrt(-mu*a)* ...
                    (1.0 -magro*alpha) ) );
                xold= sign(dtsec) * sqrt(-a) *log(temp);
            end
        end

        ktr= 1;
        dtnew = -10.0;
        % conv for dtsec to x units
        tmp = 1.0 / smu;

        while ((abs(dtnew*tmp - dtsec) >= small) && (ktr < numiter))
            xoldsqrd = xold*xold;
            znew     = xoldsqrd * alpha;

            % ------------- find c2 and c3 functions --------------
            [c2new, c3new] = findc2c3( znew );

            % ------- use a newton iteration for new values -------
            rval = xoldsqrd*c2new + rdotv*tmp *xold*(1.0 -znew*c3new) + ...
                magro*( 1.0  - znew*c2new );
            dtnew= xoldsqrd*xold*c3new + rdotv*tmp*xoldsqrd*c2new + ...
                magro*xold*( 1.0  - znew*c3new );

            % ------------- calculate new value for x -------------
            temp1 = ( dtsec*smu - dtnew ) / rval;
            xnew = xold + temp1;

            % ----- check if the univ param goes negative. if so, use bissection
            if (xnew < 0.0 && dtsec > 0.0)
                xnew = xold*0.5;
            end

            if show =='y'
                fprintf(1,'kep %3i %11.7f %11.7f %11.7f %11.7f %11.7f \n', ...
                    ktr,xold,znew,rval,xnew,dtnew*tmp);
                fprintf(1,'%3i %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f \n', ...
                    ktr,xold/sqrt(re),znew,c2new, c3new, rval/re,xnew/sqrt(re),dtnew/sqrt(mu));
            end

            ktr = ktr + 1;
            xold = xnew;
        end

        if ( ktr >= numiter )
            errork= 'knotconv';
            fprintf(1,'kep not conv in %2i iter %11.3f \n',numiter, dtseco );
            for i= 1 : 3
                v(i)= 0.0;
                r(i)= v(i);
            end
        else
            % --- find position and velocity vectors at new time --
            xnewsqrd = xnew*xnew;
            f = 1.0  - ( xnewsqrd*c2new / magro );
            g = dtsec - xnewsqrd*xnew*c3new/sqrt(mu);

            for i= 1 : 3
                r(i)= f*ro(i) + g*vo(i);
            end
            magr = mag( r );
            gdot = 1.0  - ( xnewsqrd*c2new / magr );
            fdot = ( sqrt(mu)*xnew / ( magro*magr ) ) * ( znew*c3new-1.0  );
            for i= 1 : 3
                v(i)= fdot*ro(i) + gdot*vo(i);
            end

            temp= f*gdot - fdot*g;
            if ( abs(temp-1.0 ) > 0.00001  )
                errork= 'fandg';
            end

            if show =='y'
                fprintf(1,'f %16.8f g %16.8f fdot %16.8f gdot %16.8f \n',f, g, fdot, gdot );
                tusec = sqrt(re^3/mu);
                fprintf(1,'f %16.8f g %16.8f fdot %16.8f gdot %16.8f \n',f, g/tusec, fdot*tusec, gdot );
                fprintf(1,'r1 %16.8f %16.8f %16.8f km \n',r );
                fprintf(1,'v1 %16.8f %16.8f %16.8f km/s \n',v );
                fprintf(1,'r1 %16.8f %16.8f %16.8f ER \n',r/re );
                fprintf(1,'v1 %16.8f %16.8f %16.8f ER/TU \n',v/velkmps );
            end
        end  % if
    else
        % ----------- set vectors to incoming since 0 time --------
        for i=1:3
            r(i)= ro(i);
            v(i)= vo(i);
        end
    end

end