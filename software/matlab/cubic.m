% ------------------------------------------------------------------------------
%
%                           function cubic
%
%  this function solves for the three roots of a cubic equation.  there are
%    no restrictions on the coefficients, and imaginary results are passed
%    out as separate values.  the general form is y = ax3 + bx2 + cx + d0.  note
%    that r1i will always be zero since there is always at least one real root.
%
%  author        : david vallado                  719-573-2600    1 mar 2001
%
%  revisions
%    vallado     - convert to matlab              719-573-2600   18 dec 2002
%
%  inputs          description                    range / units
%    a3          - coefficient of x cubed term
%    b2          - coefficient of x squared term
%    c1          - coefficient of x term
%    d0          - constant
%    opt         - option for output              I all roots including imaginary
%                                                 R only real roots
%                                                 U only unique real roots (no repeated)
%
%  outputs       :
%    r1r         - real portion of root 1
%    r1i         - imaginary portion of root 1
%    r2r         - real portion of root 2
%    r2i         - imaginary portion of root 2
%    r3r         - real portion of root 3
%    r3i         - imaginary portion of root 3
%
%  locals        :
%    temp1       - temporary value
%    temp2       - temporary value
%    p           - coefficient of x squared term where x cubed term is 1.0
%    q           - coefficient of x term where x cubed term is 1.0
%    r           - coefficient of constant term where x cubed term is 1.0
%    delta       - discriminator for use with cardans formula
%    e0          - angle holder for trigonometric solution
%    phi         - angle used in trigonometric solution
%    cosphi      - cosine of phi
%    sinphi      - sine of phi
%
%  coupling      :
%    quadric     - roots of second order polynomial
%
%  references    :
%    vallado       2007, 975
%
% [r1r,r1i,r2r,r2i,r3r,r3i] = cubic ( a3, b2, c1, d0, opt );
% ------------------------------------------------------------------------------

function [r1r,r1i,r2r,r2i,r3r,r3i] = cubic ( a3, b2, c1, d0, opt )
    % --------------------  implementation   ----------------------
    rad       = 180.0/pi;
    onethird  = 1.0 /3.0;
    small     = 0.0000000001;
    r1r  = 0.0;
    r1i  = 0.0;
    r2r  = 0.0;
    r2i  = 0.0;
    r3r  = 0.0;
    r3i  = 0.0;

    if (abs(a3) > small)
        % ----------- force coefficients into std form ----------------
        p= b2/a3;
        q= c1/a3;
        r= d0/a3;

        a3= onethird*( 3.0 *q - p*p );
        b2= (1.0 /27.0 )*( 2.0 *p*p*p - 9.0 *p*q + 27.0 *r );

        delta= (a3*a3*a3/27.0 ) + (b2*b2*0.25 );

        % ------------------ use cardans formula ----------------------
        if ( delta > small )
            temp1= (-b2*0.5 )+sqrt(delta);
            temp2= (-b2*0.5 )-sqrt(delta);
            temp1= sign(temp1)*abs(temp1)^onethird;
            temp2= sign(temp2)*abs(temp2)^onethird;
            r1r= temp1 + temp2 - p*onethird;

            if (opt=='I')
                r2r= -0.5 *(temp1 + temp2) - p*onethird;
                r2i= -0.5 *sqrt( 3.0  )*(temp1 - temp2);
                r3r= -0.5 *(temp1 + temp2) - p*onethird;
                r3i= -r2i;
            else
                r2r= 99999.9;
                r3r= 99999.9;
            end
        else
            % --------------- evaluate zero point ---------------------
            if ( abs( delta ) < small  )
                r1r= -2.0*sign(b2)*abs(b2*0.5)^onethird - p*onethird;
                r2r=      sign(b2)*abs(b2*0.5)^onethird - p*onethird;
                if (opt=='U')
                    r3r= 99999.9;
                else
                    r3r= r2r;
                end
            else
                % ------------ use trigonometric identities -----------
                e0     = 2.0 *sqrt(-a3*onethird);
                cosphi = (-b2/(2.0 *sqrt(-a3*a3*a3/27.0 )) );
                sinphi = sqrt( 1.0 -cosphi*cosphi );
                phi    = atan2( sinphi,cosphi );
                if (phi < 0.0)
                    phi = phi + 2.0*pi;
                end
                r1r= e0*cos( phi*onethird ) - p*onethird;
                r2r= e0*cos( phi*onethird + 120.0 /rad ) - p*onethird;
                r3r= e0*cos( phi*onethird + 240.0 /rad ) - p*onethird;
            end
        end
    else
        [r1r,r1i,r2r,r2i] =  quadric   ( b2,c1,d0,opt );
        r3r  = 99999.9;
        r3i  = 99999.9;
    end

