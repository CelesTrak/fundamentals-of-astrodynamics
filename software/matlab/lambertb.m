%  ------------------------------------------------------------------------------
%
%                           procedure lamberbattin
%
%  this procedure solves lambert's problem using battins method. the method is
%    developed in battin (1987).
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    r1          - ijk position vector 1                      km
%    r2          - ijk position vector 2                      km
%    v1          - ijk velocity vector 1 if avail             km/s
%    dm          - direction of motion                       'L', 'S'
%    de          - orbital energy                            'L', 'H'
%                  only affects nrev >= 1 solutions
%    dtsec       - time between r1 and r2                     sec
%    nrev        - number of revs to complete                 0, 1, 2, 3,
%    altpad      - altitude pad for hitearth calc             km
%    show        - control output don't output for speed      'y', 'n'
%
%  outputs       :
%    v1t         - ijk transfer velocity vector               km/s
%    v2t         - ijk transfer velocity vector               km/s
%    hitearth    - flag if hti or not                         'y', 'n'
%    errorsum    - error flag                                 'ok',
%    errorout    - text for iterations / last loop
%
%  references    :
%    vallado       2022, 505, Alg 61, ex 7-5
%    thompson      AAS GNC 2018
%
% [v1dv, v2dv, errorb] = lambertb ( r1, r2, v1, dm, de, nrev, dtsec );
% ------------------------------------------------------------------------------

function [v1t, v2t, errorb] = lambertb ( r1, r2, v1, dm, de, nrev, dtsec )
    constmath;
    constastro;

    errorb = '      ok';
    show = 'n';
    y = 0.0;
    k2 = 0.0;
    u = 0.0;
    v1t = [1000,1000,1000];
    v2t = [1000,1000,1000];

    magr1 = mag(r1);
    magr2 = mag(r2);

    cosdeltanu= dot(r1,r2)/(magr1*magr2);
    % make sure it's not more than 1.0
    if (abs(cosdeltanu) > 1.0)
        cosdeltanu = 1.0 * sign(cosdeltanu);
    end

    rcrossr = cross( r1,r2 );
    magrcrossr = mag(rcrossr);
    if dm == 'S'  % old dm == 's'  or df == 'd'
        sindeltanu= magrcrossr/(magr1*magr2);
    else
        sindeltanu= -magrcrossr/(magr1*magr2);
    end
    dnu   = atan2( sindeltanu, cosdeltanu );
    % the angle needs to be positive to work for the long way
    if dnu < 0.0
        dnu = 2.0*pi + dnu;
    end

    % these are the same
    chord= sqrt( magr1*magr1 + magr2*magr2 - 2.0*magr1*magr2*cosdeltanu );
    % chord= mag(r2-r1);

    s    = (magr1 + magr2 + chord)*0.5;
    ror   = magr2/magr1;
    eps   = ror - 1.0;

    lam = 1.0/s * sqrt(magr1*magr2) * cos(dnu*0.5);
    L = ((1.0 - lam)/(1.0 + lam))^2;
    m = 8.0*mu*dtsec*dtsec / (s^3*(1.0 + lam)^6);
    %        tan2w = 0.25*eps*eps / (sqrt(ror) + ror * (2.0 + sqrt(ror) ) );
    %        rp    = sqrt(magr1*magr2)*( (cos(dnu*0.25))^2 + tan2w );
    %        if ( dnu < pi )
    %            L = ( (sin(dnu*0.25))^2 + tan2w ) / ( (sin(dnu*0.25))^2 + tan2w + cos(dnu*0.5) );
    %        else
    %            L = ( (cos(dnu*0.25))^2 + tan2w - cos(dnu*0.5) ) / ( (cos(dnu*0.25))^2 + tan2w );
    %        end
    %        m    = mu * dtsec*dtsec / ( 8.0*rp*rp*rp );
    % initial guess
    if (nrev > 0)
        xn = 1.0 + 4.0*L;
    else
        xn   = L;   %l    % 0.0 for par and hyp, l for ell
    end

    %    lim1 = sqrt(m/L);
    % alt approach for high energy (long way, retro multi-rev) case
    if (de == 'H') && (nrev > 0)  % old dm == 'l'
        xn = 1e-20;  % be sure to reset this here!!
        x    = 10.0;  % starting value
        loops = 1;
        while ((abs(xn-x) >= small) && (loops <= 20))
            x = xn;
            temp = 1.0 / (2.0*(L - x*x));
            temp1 = sqrt(x);
            temp2 = (nrev*pi*0.5 + atan(temp1)) / temp1;
            h1   = temp * (L + x) * (1.0 + 2.0*x + L);
            h2   = temp * m * temp1 * ((L - x*x) * temp2 - (L + x));

            b = 0.25*27.0*h2 / ( (temp1*(1.0+h1))^3 );
            if b < 0.0 % reset the initial condition
                f = 2.0 * cos(1.0/3.0 * acos(sqrt(b + 1.0)));
            else
                A = (sqrt(b) + sqrt(b+1.0))^(1.0/3.0);
                f = A + 1.0/A;
            end

            y  = 2.0/3.0 * temp1 * (1.0 + h1) *(sqrt(b + 1.0) / f + 1.0 );
            xn = 0.5 * ((m/(y*y) - (1.0 + L)) - sqrt((m/(y*y) - (1.0 + L))^2 - 4.0*L));

            if isnan(y)
                y = 75.0;
                xn = 1.0;
            end
            fprintf(1,' %3i yh %11.6f x %11.6f h1 %11.6f h2 %11.6f b %11.6f f %11.7f \n',loops, y, x, h1, h2, b, f );
            loops = loops + 1;
        end  % while

        fprintf(1,' %3i yh %11.6f x %11.6f h1 %11.6f h2 %11.6f b %11.6f f %11.7f \n',loops, y, x, h1, h2, b, f );
        x = xn;
        a = s*(1.0 + lam)^2*(1.0 + x)*(L + x) / (8.0*x);
        p = (2.0*magr1*magr2*(1.0 + x)*sin(dnu*0.5)^2) / (s*(1 + lam)^2 * (L + x));  % thompson (1.0 + x)*
        ecc = sqrt(1.0 - p/a);
        [v1t, v2t] = lambhodograph( r1, v1, r2, p, ecc, dnu, dtsec );
        fprintf(1,'high v1t %16.8f %16.8f %16.8f \n',v1t );
    else
        % standard processing
        % note that the dr nrev=0 case is not represented
        loops= 1;
        y1 = 0.0;
        x    = 10.0;  % starting value
        while ((abs(xn-x) >= small) && (loops <= 30))
            if (nrev > 0)
                x = xn;
                temp = 1.0 / ( (1.0 + 2.0*x + L) * (4.0*x^2) );
                temp1 = (nrev*pi*0.5 + atan(sqrt(x))) / sqrt(x);
                h1   = temp * (L + x)^2 * (3.0*(1.0 + x)^2 * temp1 - (3.0 + 5.0*x));
                h2   = temp * m * ((x*x - x*(1.0 + L) - 3.0*L) * temp1 + (3.0*L +x));
            else
                x    = xn;
                tempx  = seebatt(x);
                denom= 1.0 / ( (1.0 + 2.0*x + L) * (4.0*x + tempx*(3.0 + x) ) );
                h1   = (L + x)^2 * (1.0 + 3.0*x + tempx)*denom;
                h2   = m*(x - L + tempx)*denom;
            end

            % ----------------------- evaluate cubic ------------------
            b = 0.25*27.0*h2 / ((1.0+h1)^3 );

            %        if b < -1.0 % reset the initial condition
            %fprintf(1,'xx %11.6f  %11.6f  %11.6f \n',L, xn, b);
            %            xn = 1.0 - 2.0*L
            %        end
            %        else
            %            if y1 > lim1
            %                xn = xn * (lim1/y1)
            %            end
            %            else
            u  = 0.5*b / ( 1.0 + sqrt(1.0 + b) );
            k2 = kbatt(u);
            y  = ( (1.0 + h1) / 3.0 )*(2.0 + sqrt(1.0 + b) / (1.0 + 2.0*u*k2*k2) );
            xn= sqrt( ( (1.0 - L)*0.5 )^2 + m/(y*y) ) - (1.0 + L)*0.5;
            %                    xn = sqrt(l*l + m/(y*y)) - (1.0 - l); alt, doesn't seem to work
            %            end
            %        end

            y1=  sqrt( m/((L + x)*(1.0 + x)) );

            if isnan(y)
                y = 75.0;
                xn = 1.0;
            end

            loops = loops + 1;
            if show == 'y'
                fprintf(1,' %3i yb %11.6f x %11.6f k2 %11.6f b %11.6f u %11.6f y1 %11.7f \n',loops,y, x, k2, b, u, y1 );
            end
        end  % while
        if show == 'y'
            fprintf(1,' %3i yb %11.6f x %11.6f k2 %11.6f b %11.6f u %11.6f y1 %11.7f \n',loops,y, x, k2, b, u, y1 );
        end

        if (loops < 30)
            % blair approach use y from solution
            %       lam = 1.0/s * sqrt(magr1*magr2) * cos(dnu*0.5);
            %       m = 8.0*mu*dtsec*dtsec / (s^3*(1.0 + lam)^6);
            %       L = ((1.0 - lam)/(1.0 + lam))^2;
            %a = s*(1.0 + lam)^2*(1.0 + x)*(lam + x) / (8.0*x);
            % p = (2.0*magr1*magr2*(1.0 + x)*sin(dnu*0.5)^2)^2 / (s*(1 + lam)^2*(lam + x));  % loechler, not right?
            p = (2.0*magr1*magr2*y*y*(1.0 + x)^2*sin(dnu*0.5)^2) / (m*s*(1 + lam)^2);  % thompson
            ecc = sqrt( (eps^2 + 4.0*magr2/magr1*sin(dnu*0.5)^2*((L-x)/(L+x))^2) / (eps^2 + 4.0*magr2/magr1*sin(dnu*0.5)^2)  );
            [v1t, v2t] = lambhodograph( r1, v1, r2, p, ecc, dnu, dtsec );
            %            fprintf(1,'oldb v1t %16.8f %16.8f %16.8f %16.8f\n',v1dv, mag(v1dv) );
            %         r_180 = 0.001;  % 1 meter
            %         [v1dvh,v2dvh] = lambert_vel(r1, v1, r2, dnu, p, ecc, mu, dtsec, r_180);
            %         fprintf(1,'newb v1t %16.8f %16.8f %16.8f %16.8f\n',v1dvh, mag(v1dvh) );

            % Battin solution to orbital parameters (and velocities)
            % thompson 2011, loechler 1988
            % if dnu > pi
            %     lam = -sqrt((s-chord)/s);
            % else
            %     lam = sqrt((s-chord)/s);
            % end
            % %      x = xn;
            %
            % % loechler pg 21 seems correct!
            % v1dvl = 1.0/(lam*(1.0 + lam))*sqrt(mu*(1.0+x)/(2.0*s^3*(L + x)))*((r2-r1) + s*(1.0+lam)^2*(L + x)/(magr1*(1.0 + x))*r1);
            % % added v2
            % v2dvl = 1.0/(lam*(1.0 + lam))*sqrt(mu*(1.0+x)/(2.0*s^3*(L + x)))*((r2-r1) - s*(1.0+lam)^2*(L + x)/(magr2*(1.0 + x))*r2);
            % %fprintf(1,'loe v1t %16.8f %16.8f %16.8f %16.8f\n',v1dvl, mag(v1dvl) );
            %fprintf(1,'loe v2t %16.8f %16.8f %16.8f %16.8f\n',v2dvl, mag(v2dvl) );
        end  % if loops converged < 30
    end  % if de and nrev

end  % lambertb


% ------- two recursion algorithms needed by the lambertbattin routine
function kbatt = kbat( v )

    % -------------------------  implementation   -------------------------
    d(1) =     1.0 /    3.0;
    d(2) =     4.0 /   27.0;
    d(3) =     8.0 /   27.0;
    d(4) =     2.0 /    9.0;
    d(5) =    22.0 /   81.0;
    d(6) =   208.0 /  891.0;
    d(7) =   340.0 / 1287.0;
    d(8) =   418.0 / 1755.0;
    d(9) =   598.0 / 2295.0;
    d(10)=   700.0 / 2907.0;
    d(11)=   928.0 / 3591.0;
    d(12)=  1054.0 / 4347.0;
    d(13)=  1330.0 / 5175.0;
    d(14)=  1480.0 / 6075.0;
    d(15)=  1804.0 / 7047.0;
    d(16)=  1978.0 / 8091.0;
    d(17)=  2350.0 / 9207.0;
    d(18)=  2548.0 /10395.0;
    d(19)=  2968.0 /11655.0;
    d(20)=  3190.0 /12987.0;
    d(21)=  3658.0 /14391.0;

    % ----------------- process forwards ------------------------
    sum1    = d(1);
    delold  = 1.0;
    termold = d(1);
    i = 2;

    ktr = 21;
    while ((i <= ktr) && (abs(termold) > 0.00000001 ))
        del  = 1.0 / ( 1.0 + d(i)*v*delold );
        term = termold * ( del - 1.0 );
        sum1 = sum1 + term;
        i    = i + 1;
        delold = del;
        termold= term;
    end

    sum2  = 0.0;
    term2 = 1.0 + d(ktr)*v;
    for i=1:ktr-2
        sum2  = d(ktr-i)*v / term2;
        term2 = 1.0 + sum2;
    end

    kbatt = d(1) / term2;

    %            test=  d(1) / ...
    %                   (1 + (d(2)*v / ...
    %                         (1 + (d(3)*v / ...
    %                               (1 + (d(4)*v / ...
    %                                     (1 + (d(5)*v / ...
    %                                           (1 + (d(6)*v / ...
    %                                                 (1 + (d(7)*v / ...
    %                                                       (1 + (d(8)*v / ...
    %                                                             (1 + (d(9)*v / ...
    %                                                                  (1 + (d(10)*v / ...
    %                                                                        (1 + (d(11)*v )))))))) ...
    %                                                                         ))))))))))));
    %kbatt = test
    % form coefficients
    % for i = 1:12
    %     j= i*2;
    %     fprintf(1,'K %3d %14d %14d \n',j, 2.0*(3*i+1)*(6*i-1), 9.0*(4*i-1)*(4*i+1) );
    %     j= i*2+1;
    %     fprintf(1,'K %3d %14d %14d \n',j, 2.0*(3*i+2)*(6*i+1), 9.0*(4*i+1)*(4*i+3) );
    % end
    %
end


% ------- two recursion algorithms needed by the lambertbattin routine
function seebatt = seebatt( v )

    % -------------------------  implementation   -------------------------
    c(1) =    0.2;
    c(2) =    9.0 /  35.0;
    c(3) =   16.0 /  63.0;
    c(4) =   25.0 /  99.0;
    c(5) =   36.0 / 143.0;
    c(6) =   49.0 / 195.0;
    c(7) =   64.0 / 255.0;
    c(8) =   81.0 / 323.0;
    c(9) =  100.0 / 399.0;
    c(10)=  121.0 / 483.0;
    c(11)=  144.0 / 575.0;
    c(12)=  169.0 / 675.0;
    c(13)=  196.0 / 783.0;
    c(14)=  225.0 / 899.0;
    c(15)=  256.0 /1023.0;
    c(16)=  289.0 /1155.0;
    c(17)=  324.0 /1295.0;
    c(18)=  361.0 /1443.0;
    c(19)=  400.0 /1599.0;
    c(20)=  441.0 /1763.0;
    c(21)=  484.0 /1935.0;
    sqrtopv= sqrt(1.0 + v);
    eta    = v / ( 1.0 + sqrtopv )^2;

    % ------------------- process forwards ----------------------
    delold = 1.0;
    termold= c(1);   % * eta
    sum1   = termold;
    i= 2;
    while ((i <= 21) && (abs(termold) > 0.00000001 ))
        del  = 1.0 / ( 1.0 + c(i)*eta*delold );
        term = termold * (del - 1.0);
        sum1 = sum1 + term;
        i    = i + 1;
        delold = del;
        termold= term;
    end

    seebatt= 1.0/ ((1.0/(8.0*(1.0+sqrtopv))) * ( 3.0 + sum1 / ( 1.0+eta*sum1 ) ) );

    c(1) =    9.0 /   7.0;
    c(2) =   16.0 /  63.0;
    c(3) =   25.0 /  99.0;
    c(4) =   36.0 / 143.0;
    c(5) =   49.0 / 195.0;
    c(6) =   64.0 / 255.0;
    c(7) =   81.0 / 323.0;
    c(8) =  100.0 / 399.0;
    c(9) =  121.0 / 483.0;
    c(10)=  144.0 / 575.0;
    c(11)=  169.0 / 675.0;
    c(12)=  196.0 / 783.0;
    c(13)=  225.0 / 899.0;
    c(14)=  256.0 /1023.0;
    c(15)=  289.0 /1155.0;
    c(16)=  324.0 /1295.0;
    c(17)=  361.0 /1443.0;
    c(18)=  400.0 /1599.0;
    c(19)=  441.0 /1763.0;
    c(20)=  484.0 /1935.0;


    ktr   = 20;
    sum2  = 0.0;
    term2 = 1.0 + c(ktr)*eta;
    for i = 1 : ktr - 2
        sum2  = c(ktr-i)*eta / term2;
        term2 = 1.0 + sum2;
    end

    seebatt = 8.0*(1.0 + sqrtopv) / ...
        (3.0 + ...
        (1.0 / ...
        (5.0 + eta + ((9.0/7.0)*eta/term2 ) ) ) );

    % form coefficients
    % for i= 3:24
    %     j= i;
    %     fprintf(1,'E %3d %14d %14d \n',j, i*i*1.0, (2.0*i)^2-1 );
    % end

end
