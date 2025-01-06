%
% ----------------------------------------------------------------------------
%
%                           function fundarg
%
%  this function calulates the delauany variables and planetary values for
%  several theories.
%
%  author        : david vallado                  719-573-2600   16 jul 2004
%
%  inputs          description                    range / units
%    ttt         - julian centuries of tt
%    opt         - method option                  '06', '96', '80'
%
%  outputs       :
%    l           - delaunay element               rad
%    ll          - delaunay element               rad
%    f           - delaunay element               rad
%    d           - delaunay element               rad
%    omega       - delaunay element               rad
%    planetary longitudes                         rad
%
%  locals        :
%    ttt2,ttt3,  - powers of ttt
%
%  coupling      :
%    none        -
%
%  references    :
%    vallado       vallado       2022, 210-212, 226
%
% [ l, l1, f, d, omega, ...
%   lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate 
%
%   [fArgs] = fundarg( ttt, opt );
% ----------------------------------------------------------------------------

function [fArgs] = fundarg( ttt, opt )

        sethelp;

        twopi = pi * 2.0;
        deg2rad = pi/180.0;
        % arcsec to deg
        oo3600 = 1.0 / 3600.0;

        % ---- determine coefficients for iau 2000 nutation theory ----
        % ---- iau 2006 theory
        if opt == '06'
                % ------ form the delaunay fundamental arguments in ", converted to deg
                l = ((((-0.00024470 * ttt + 0.051635) * ttt + 31.8792) * ttt + 1717915923.2178) * ttt ...
                    + 485868.249036) * oo3600;
                l1 = ((((-0.00001149 * ttt + 0.000136) * ttt - 0.5532) * ttt + 129596581.0481) * ttt ...
                    + 1287104.793048) * oo3600;
                f = ((((+0.00000417 * ttt - 0.001037) * ttt - 12.7512) * ttt + 1739527262.8478) * ttt ...
                    + 335779.526232) * oo3600;
                d = ((((-0.00003169 * ttt + 0.006593) * ttt - 6.3706) * ttt + 1602961601.2090) * ttt ...
                    + 1072260.703692) * oo3600;
                omega = ((((-0.00005939 * ttt + 0.007702) * ttt + 7.4722) * ttt - 6962890.5431) * ttt ...
                    + 450160.398036) * oo3600;

                % ------ form the planetary arguments in ", converted to deg
                %lonmer = (908103.259872 + 538101628.688982 * ttt) * oo3600;
                %lonven = (655127.283060 + 210664136.433548 * ttt) * oo3600;
                %lonear = (361679.244588 + 129597742.283429 * ttt) * oo3600;
                %lonmar = (1279558.798488 + 68905077.493988 * ttt) * oo3600;
                %lonjup = (123665.467464 + 10925660.377991 * ttt) * oo3600;
                %lonsat = (180278.799480 + 4399609.855732 * ttt) * oo3600;
                %lonurn = (1130598.018396 + 1542481.193933 * ttt) * oo3600;
                %lonnep = (1095655.195728 + 786550.320744 * ttt) * oo3600;
                %precrate = ((1.112022 * ttt + 5028.8200) * ttt) * oo3600;
                % these are close (all in rad) - usually 1e-10, but some are as high as 1e-06
                % these are from TN-36
                lonmer = rem((4.402608842 + 2608.7903141574 * ttt), twopi);
                lonven = rem((3.176146697 + 1021.3285546211 * ttt), twopi);
                lonear = rem((1.753470314 + 628.3075849991 * ttt), twopi);
                lonmar = rem((6.203480913 + 334.0612426700 * ttt), twopi);
                lonjup = rem((0.599546497 + 52.9690962641 * ttt), twopi);
                lonsat = rem((0.874016757 + 21.3299104960 * ttt), twopi);
                lonurn = rem((5.481293872 + 7.4781598567 * ttt), twopi);
                lonnep = rem((5.311886287 + 3.8133035638 * ttt), twopi);
                precrate = ((0.024381750 + 0.00000538691 * ttt) * ttt);
        end


        % ---- iau 1996 theory
        if opt == '96'
            % ------ form the delaunay fundamental arguments in deg
            l = ((((-0.00024470 * ttt + 0.051635) * ttt + 31.8792) * ttt + 1717915923.2178) * ttt) * oo3600 + 134.96340251;
            l1 = ((((-0.00001149 * ttt - 0.000136) * ttt - 0.5532) * ttt + 129596581.0481) * ttt) * oo3600 + 357.52910918;
            f = ((((+0.00000417 * ttt + 0.001037) * ttt - 12.7512) * ttt + 1739527262.8478) * ttt) * oo3600 + 93.27209062;
            d = ((((-0.00003169 * ttt + 0.006593) * ttt - 6.3706) * ttt + 1602961601.2090) * ttt) * oo3600 + 297.85019547;
            omega = ((((-0.00005939 * ttt + 0.007702) * ttt + 7.4722) * ttt - 6962890.2665) * ttt) * oo3600 + 125.04455501;
            % ------ form the planetary arguments in deg
            lonmer  = 0.0;
            lonven  = 181.979800853  +  58517.8156748  *ttt;   % deg
            lonear  = 100.466448494  +  35999.3728521  *ttt;
            lonmar  = 355.433274605  +  19140.299314   *ttt;
            lonjup  =  34.351483900  +   3034.90567464 *ttt;
            lonsat  =  50.0774713998 +   1222.11379404 *ttt;
            lonurn  = 0.0;
            lonnep  = 0.0;
            precrate=   1.39697137214*ttt + 0.0003086*ttt;
        end

        % ---- iau 1980 theory
        if opt == '80'
            % ------ form the delaunay fundamental arguments in deg
            l = ((((0.064) * ttt + 31.310) * ttt + 1717915922.6330) * ttt) * oo3600 ...
                + 134.96298139;
            l1 = ((((-0.012) * ttt - 0.577) * ttt + 129596581.2240) * ttt) * oo3600 ...
                + 357.52772333;
            f = ((((0.011) * ttt - 13.257) * ttt + 1739527263.1370) * ttt) * oo3600 ...
                + 93.27191028;
            d = ((((0.019) * ttt - 6.891) * ttt + 1602961601.3280) * ttt) * oo3600 ...
                + 297.85036306;
            omega = ((((0.008) * ttt + 7.455) * ttt - 6962890.5390) * ttt) * oo3600 ...
                + 125.04452222;
            % ------ form the planetary arguments in deg
            % iers tn13 shows no planetary
            % seidelmann shows these equations
            % circ 163 shows no planetary  ???????
            lonmer = 252.3 + 149472.0 * ttt;
            lonven = 179.9 + 58517.8 * ttt;
            lonear = 98.4 + 35999.4 * ttt;
            lonmar = 353.3 + 19140.3 * ttt;
            lonjup = 32.3 + 3034.9 * ttt;
            lonsat = 48.0 + 1222.1 * ttt;
            lonurn = 0.0;
            lonnep = 0.0;
            precrate = 0.0;
        end

        % ---- convert units to rad
        l    = rem( l,360.0  )     * deg2rad; % rad
        l1   = rem( l1,360.0  )    * deg2rad;
        f    = rem( f,360.0  )     * deg2rad;
        d    = rem( d,360.0  )     * deg2rad;
        omega= rem( omega,360.0  ) * deg2rad;
        if opt ~= '06'
            lonmer= rem( lonmer,360.0 ) * deg2rad;  % rad
            lonven= rem( lonven,360.0 ) * deg2rad;
            lonear= rem( lonear,360.0 ) * deg2rad;
            lonmar= rem( lonmar,360.0 ) * deg2rad;
            lonjup= rem( lonjup,360.0 ) * deg2rad;
            lonsat= rem( lonsat,360.0 ) * deg2rad;
            lonurn= rem( lonurn,360.0 ) * deg2rad;
            lonnep= rem( lonnep,360.0 ) * deg2rad;
            precrate= rem( precrate,360.0 ) * deg2rad;
        end
        %iauhelp='y';
        if iauhelp == 'y'
            fprintf(1,'fa %11.7f  %11.7f  %11.7f  %11.7f  %11.7f deg \n',l*180/pi,l1*180/pi,f*180/pi,d*180/pi,omega*180/pi );
            fprintf(1,'fa %11.7f  %11.7f  %11.7f  %11.7f deg \n',lonmer*180/pi,lonven*180/pi,lonear*180/pi,lonmar*180/pi );
            fprintf(1,'fa %11.7f  %11.7f  %11.7f  %11.7f deg \n',lonjup*180/pi,lonsat*180/pi,lonurn*180/pi,lonnep*180/pi );
            fprintf(1,'fa %11.7f  \n',precrate*180/pi );
        end


        fArgs(1) = l;   % delaunay variables
        fArgs(2) = l1;
        fArgs(3) = f;
        fArgs(4) = d;
        fArgs(5) = omega;
        fArgs(6) = lonmer;  % begin planetary longitudes
        fArgs(7) = lonven;
        fArgs(8) = lonear;
        fArgs(9) = lonmar;
        fArgs(10) = lonjup;
        fArgs(11) = lonsat;
        fArgs(12) = lonurn;
        fArgs(13) = lonnep;
        fArgs(14) = precrate;

        % test if they are equivalent
        % most around 1e-10, but some at 1e-6
        %         oo3600 = 1.0 / 3600.0;
        %         deg2rad = pi / 180.0;
        %         ttt = 0.34698738576;
        %         twopi = 2.0 * pi;
        %         lonmer = mod((908103.259872 + 538101628.688982 * ttt) * oo3600,360)*deg2rad;
        %         lonven = mod((655127.283060 + 210664136.433548 * ttt) * oo3600,360)*deg2rad;
        %         lonear = mod((361679.244588 + 129597742.283429 * ttt) * oo3600,360)*deg2rad;
        %         lonmar = mod((1279558.798488 + 68905077.493988 * ttt) * oo3600,360)*deg2rad;
        %         lonjup = mod((123665.467464 + 10925660.377991 * ttt) * oo3600,360)*deg2rad;
        %         lonsat = mod((180278.799480 + 4399609.855732 * ttt) * oo3600,360)*deg2rad;
        %         lonurn = mod((1130598.018396 + 1542481.193933 * ttt) * oo3600,360)*deg2rad;
        %         lonnep = mod((1095655.195728 + 786550.320744 * ttt) * oo3600,360)*deg2rad;
        %         precrate = ((1.112022 * ttt + 5028.8200) * ttt) * oo3600*deg2rad;
        %
        %         lonmer1 = mod (4.402608842 + 2608.7903141574 * ttt , twopi);
        %         lonven1 = mod (3.176146697 + 1021.3285546211 * ttt , twopi);
        %         lonear1 = mod(1.753470314 + 628.3075849991 * ttt , twopi);
        %         lonmar1 = mod(6.203480913 + 334.0612426700 * ttt , twopi);
        %         lonjup1 = mod(0.599546497 + 52.9690962641 * ttt , twopi);
        %         lonsat1 = mod(0.874016757 + 21.3299104960 * ttt , twopi);
        %         lonurn1 = mod(5.481293872 + 7.4781598567 * ttt , twopi);
        %         lonnep1 = mod(5.311886287 + 3.8133035638 * ttt , twopi);
        %         precrate1 = (0.024381750 + 0.00000538691 * ttt ) *ttt;
        %
        %         lonmer-lonmer1
        %         lonven-lonven1
        %         lonear-lonear1
        %         lonmar-lonmar1
        %         lonjup-lonjup1
        %         lonsat-lonsat1
        %         lonurn-lonurn1
        %         lonnep-lonnep1
        %         precrate-precrate1
        %

