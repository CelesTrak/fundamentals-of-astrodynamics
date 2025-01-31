%  ------------------------------------------------------------------------------
%
%                           procedure angleslaplace
%
%  this procedure solves the problem of orbit determination using three
%    optical sightings and the method of laplace. the 8th order root is generally
%    the big point of discussion. A Halley iteration permits a quick solution to
%    find the correct root, with a starting guess of 20000 km. the general
%    formulation yields polynomial coefficients that are very large, and can easily
%    become overflow operations. Thus, canonical units are used only until teh root is found,
%    then regular units are resumed.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    tdecl1       - declination #1                               rad
%    tdecl2       - declination #2                               rad
%    tdecl3       - declination #3                               rad
%    trtasc1      - right ascension #1                           rad
%    trtasc2      - right ascension #2                           rad
%    trtasc3      - right ascension #3                           rad
%    jd1, jdf1    - julian date of 1st sighting                  days from 4713 bc
%    jd2, jdf2    - julian date of 2nd sighting                  days from 4713 bc
%    jd3, jdf3    - julian date of 3rd sighting                  days from 4713 bc
%    diffsites    - if sites are different (need better test)    'y', 'n'
%    rs1          - eci site position vector #1                  km
%    rs2          - eci site position vector #2                  km
%    rs3          - eci site position vector #3                  km
%
%  outputs        :
%    r2           -  position vector                             km
%    v2           -  velocity vector                             km / s
%    errstr       - output results for debugging
%
%  locals         :
%    l1           - line of sight vector for 1st
%    l2           - line of sight vector for 2nd
%    l3           - line of sight vector for 3rd
%    ldot         - 1st derivative of l2
%    lddot        - 2nd derivative of l2
%    rs2dot       - 1st derivative of rs2 - vel
%    rs2ddot      - 2nd derivative of rs2
%    t12t13       - (t1-t2) * (t1-t3)
%    t21t23       - (t2-t1) * (t2-t3)
%    t31t32       - (t3-t1) * (t3-t2)
%    i            - index
%    d            -
%    d1           -
%    d2           -
%    d3           -
%    d4           -
%    oldr         - previous iteration on r2
%    rho          - range from site to satellite at t2
%    rhodot       -
%    dmat         -
%    d1mat        -
%    d2mat        -
%    d3mat        -
%    d4mat        -
%    earthrate    - angular rotation of the earth
%    l2dotrs      - vector l2 dotted with rs
%    temp         - temporary vector
%    temp1        - temporary vector
%    small        - tolerance
%    roots        -
%
%  coupling       :
%    mag          - magnitude of a vector
%    determinant  - evaluate the determinant of a matrix
%    cross        - cross product of two vectors
%    norm         - normlize a matrix
%    sgnval       - sgn a value to a matrix
%    factor       - find the roots of a polynomial
%
%  references     :
%    vallado       2022, 441
%
% [r2, v2] = anglesl(decl1, decl2, decl3, rtasc1, rtasc2, ...
%        rtasc3, jd1, jdf1, jd2, jdf2, jd3, jdf3, diffsites, rs1, rs2, rs3, fid);
% ------------------------------------------------------------------------------

function [r2, v2] = anglesl(decl1, decl2, decl3, rtasc1, rtasc2, ...
        rtasc3, jd1, jdf1, jd2, jdf2, jd3, jdf3, diffsites, rs1, rs2, rs3, fid)

    constastro;
    % -------------------------  implementation   -------------------------
    %     omegaearth   =     0.000072921158553;  % earth rad/s
    %     omegaearth   =     0.017202791208627;  % sun rad/s
    %     omegaearth   =     2.0 * pi/365.24221897;  % au / day
    earthrate(1)= 0.0;
    earthrate(2)= 0.0;
    earthrate(3)= omegaearthradptu;
    %        tuday        =    58.132440906;
    %        mu           =     1.32712428e11;
    % need to switch these for interplanetary

    %        constant;
    small = 0.00000001;
    rad = 180.0 / pi;

    % ---------- set middle to 0, find deltas to others -----------

    % test problem%%%%%%%%%%%%%%%%%%%%%%%%%%%/
    %      los1= [-0.425365 0.777650 0.462953];  % just in km
    %      los2 =[-0.825702 0.259424 0.500914];
    %      los3 = [ -0.947067 -0.129576 0.293725];
    % los1 =[-0.425364592304 ,0.777650239833, 0.462952554914];
    % los2=[ -0.825702365309, 0.259423566604, 0.500914181287];
    % los3 =[-0.947067028031, -0.129575647726, 0.2937246941];
    % t1=-1200;  % sec
    % t2=0;
    % t3=1200;
    % tau12 = t1-t2;
    % tau13 = t1-t3;
    % tau32 = t3-t2;
    % test problem%%%%%%%%%%%%%%%%%%%%%%%%%%%/

    tau12 = (jd1 - jd2) * 86400.0 + (jdf1 - jdf2) * 86400.0; % days to sec
    tau13 = (jd1 - jd3) * 86400.0 + (jdf1 - jdf3) * 86400.0;
    tau32 = (jd3 - jd2) * 86400.0 + (jdf3 - jdf2) * 86400.0;

    % switch to cannonical for better accuracy in roots
    tau12c = tau12 / tusec;
    tau13c = tau13 / tusec;
    tau32c = tau32 / tusec;

    for i=1:3
        rs1c(i) = rs1(i) / re;
        rs2c(i) = rs2(i) / re;
        rs3c(i) = rs3(i) / re;
    end

    % --------------- find line of sight vectors ------------------
    % should be eci...
    los1(1)= cos(decl1)*cos(rtasc1);
    los1(2)= cos(decl1)*sin(rtasc1);
    los1(3)= sin(decl1);
    los2(1)= cos(decl2)*cos(rtasc2);
    los2(2)= cos(decl2)*sin(rtasc2);
    los2(3)= sin(decl2);
    los3(1)= cos(decl3)*cos(rtasc3);
    los3(2)= cos(decl3)*sin(rtasc3);
    los3(3)= sin(decl3);

    % same- they're both unit vectors
    %     l1
    %     unit(l1)

    % -------------------------------------------------------------
    %       using lagrange interpolation formula to derive an expression
    %       for l(t), substitute t=t2 and differentiate to obtain the
    %       derivatives of l.
    % -------------------------------------------------------------
    s1 = -tau32c / (tau12c * tau13c);
    s2 = (tau12c + tau32c) / (tau12c * tau32c); % be careful here! it's -t1-t3 which can be 0
    s3 = -tau12c / (-tau13c * tau32c);
    s4 = 2.0 / (tau12c * tau13c);
    s5 = 2.0 / (tau12c * tau32c);
    s6 = 2.0 / (-tau13c * tau32c);

    for i = 1:3
        ldot(i) = s1 * los1(i) + s2 * los2(i) + s3 * los3(i);  % rad / s
        lddot(i) = s4 * los1(i) + s5 * los2(i) + s6 * los3(i);  % rad / s^2
    end
    %    ldotmag = mag( ldot )
    %    lddotmag = mag( lddot )
    % should these unit vectors use a diff name????????%
    %    ldot = unit( ldot );
    %    lddot = unit( lddot );
    %    ldot
    %    lddot

    % ------------------- find 2nd derivative of rs ---------------
    %     temp = cross( rs1,rs2 );
    %     temp1 = cross( rs2,rs3 );
    %
    %     %      needs a different test xxxx%%
    %     if ( ( mag(temp) > small ) & ( mag( temp1) > small )  )
    %         %          fix this test here

    if diffsites == 'n'
        % ------------ all sightings from one site -----------------
        [rs2cdot] = cross( earthrate, rs2c );
        [rs2cddot] = cross( earthrate, rs2cdot );
    else
        % ---------- each sighting from a different site ----------
        for i= 1 : 3
            rs2cdot(i) = s1 * rs1c(i) + s2 * rs2c(i) + s3 * rs3c(i);
            rs2cddot(i) = s4 * rs1c(i) + s5 * rs2c(i) + s6 * rs3c(i);
        end
    end

    rs2cdot;
    rs2cddot;

    for i= 1 : 3
        dmat(i,1) = los2(i);
        dmat(i,2) = ldot(i);
        dmat(i,3) = lddot(i);
        dmat;
        % ----------------  position determinants -----------------
        dmat1(i,1) =los2(i);
        dmat1(i,2) =ldot(i);
        dmat1(i,3) =rs2cddot(i);
        dmat2(i,1) =los2(i);
        dmat2(i,2) =ldot(i);
        dmat2(i,3) =rs2c(i);

        % ------------  velocity determinants ---------------------
        dmat3(i,1) =los2(i);
        dmat3(i,2) =rs2cddot(i);
        dmat3(i,3) =lddot(i);
        dmat4(i,1) =los2(i);
        dmat4(i,2) =rs2c(i);
        dmat4(i,3) =lddot(i);
    end

    dmat1;
    dmat2;
    d = 2.0 * det(dmat);
    d1c= det(dmat1);
    d2c= det(dmat2);
    d3c= det(dmat3);
    d4c= det(dmat4);
    fprintf(1,'d, di  %11.6g %11.6g %11.6g %11.6g %11.6g \n', d, d1c, d2c, d3c, d4c);
    %
    % ---------------  iterate to find rho magnitude ----------------
    %     magr= 1.5   % first guess
    %     writeln( 'input initial guess for magr ' )
    %     readln( magr )
    %     i= 1
    %     repeat
    %         oldr= magr
    %         rho= -2.0*d1/d - 2.0*d2/(magr*magr*magr*d)
    %         magr= sqrt( rho*rho + 2.0*rho*l2dotrs + rs2(4)*rs2(4) )
    %         inc(i)
    %         magr= (oldr - magr ) / 2.0             % simple bissection
    %         writeln( fileout,'rho guesses ',i:2,'rho ',rho:14:7,' magr ',magr:14:7,oldr:14:7 )
    % seems to converge, but wrong numbers
    %         inc(i)
    %     until ( abs( oldr-magr ) < small ) | ( i .ge. 30 )

    if ( abs(d) > 1.0-14 )
        % --------------- solve eighth order poly -----------------
        l2dotrs= dot( los2,rs2c );
        poly( 1)=  1.0;  % r2^8th variable%%%%%%%%%%%%%%
        poly( 2)=  0.0;
        poly( 3)=  (l2dotrs*4.0*d1c/d - 4.0*d1c*d1c/(d*d) ...
            - mag(rs2c)*mag(rs2c) );
        poly( 4)=  0.0;
        poly( 5)=  0.0;
        poly( 6)=  (l2dotrs*4.0*d2c/d - 8.0*d1c*d2c/(d*d) );
        poly( 7)=  0.0;
        poly( 8)=  0.0;
        poly( 9)=  -4.0*d2c*d2c/(d*d);
        rootarr = roots( poly );

        x = rootarr(1,1);
        poly(1)*x^8 + poly(3)*x^6+poly(6)*x^3+poly(9);


        % ------------------ find correct (xx) root ----------------
        % bigr2= 0.0;
        % for j= 1 : 8
        %     %                if ( abs( roots(j,2) ) < small )
        %     %                    writeln( 'root ',j,roots(j,1),' + ',roots(j,2),'j' )
        %     %        temproot= roots(j,1)*roots(j,1)
        %     %        temproot= temproot*temproot*temproot*temproot +
        %     %                  poly(3)*temproot*temproot*temproot + poly(6)*roots(j,1)*temproot + poly(9)
        %     %                    writeln( fileout,'root ',j,roots(j,1),' + ',roots(j,2),'j  value = ',temproot )
        %     if ( rootarr(j,1) > bigr2 )
        %         bigr2= rootarr(j,1);
        %     end
        %     %                  end
        % end
        bigr2 = 100.0;

        kk = 0;
        bigr2c = 20000.0 / re; % er guess ~GPS altitude
        % bigr2nx = bigr2c;
        while (abs(bigr2 - bigr2c) > 8.0e-5 && kk < 15)  % do in er, 0.5 km

            bigr2 = bigr2c;
            %    bigr2x = bigr2nx;
            deriv = bigr2^8 + poly(3) * bigr2^6 + poly(6) * bigr2^3 + poly(9);
            deriv1 = 8.0 * bigr2^7 + 6.0 * poly(3) * bigr2^5 + 3.0 * poly(6) *bigr2^2;
            deriv2 = 56.0 * bigr2^6 + 30.0 * poly(3) * bigr2^4 + 6.0 * poly(6) * bigr2;
            % just use Halley
            %bigr2n = bigr2 - deriv * n2 / (deriv1 * (n2 - deriv * deriv2 * 0.5));
            % Halley iteration
            bigr2c = bigr2 - (2.0 * deriv * deriv1) / (2.0 * deriv1 * deriv1 - deriv * deriv2);
            kk = kk + 1;
        end

        % if (bigr2c < 0.0 || bigr2c * re > 50000.0)
        %     bigr2c = 35000.0 / re;  % simply set this to about GEO, allowing for less than that too.
        % end
        % fprintf(1,'bigr2 %11.7f  %11.7f er \n',bigr2, bigr2/re);
        %  fprintf(1,'keep this root ? ');
        %        input (bigr2);

        %rho= -2.0*d1/d - 2.0*mu*d2 / (bigr2*bigr2*bigr2*d);
        % still in canonical units
        rho = -2.0 * d1c / d - 2.0 * d2c / (bigr2c * bigr2c * bigr2c * d);
        %rho = -2.0 * d1 / d - gravConst.mu * 2.0 * d2 / (bigr2 * bigr2 * bigr2 * d);

        % ---- find the middle position vector ----
        for k=1:3
            r2c(k) = rho * los2(k) + rs2c(k);
        end
        fprintf(1,'rho %11.7f  %11.7f er \n',rho, rho/re);

        rhodot = -d3c / d - d4c / (bigr2c * bigr2c * bigr2c * d);
        %rhodot = -d3 / d - gravConst.mu * d4 / (bigr2 * bigr2 * bigr2 * d);

        % -----find middle velocity vector---- -
        for i=1:3
            v2c(i) = rhodot * los2(i) + rho * ldot(i) + rs2cdot(i);
        end
        % now convert back to normal units
        for i=1:3
            r2(i) = r2c(i) * re;
            v2(i) = v2c(i) * re / tusec;
        end

        % -------------- find middle velocity vector --------------
        for i= 1 : 3
            v2(i)= rhodot*los2(i) + rho*ldot(i) + rs2cdot(i);
        end
    else
        fprintf(1,'determinant value was zero %11.7f ',d );
    end

end