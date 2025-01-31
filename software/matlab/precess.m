% -----------------------------------------------------------------------------
%
%                           function precess
%
%  this function calculates the transformation matrix that accounts for the effects
%    of precession. both the 1980 and 2006 theories are handled. note that the
%    required parameters differ a little.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    ttt         - julian centuries of tt
%    opt         - method option                           e80, e96, e06cio
%
%  outputs       :
%    psia        - cannonical precession angle                    rad    (00 only)
%    wa          - cannonical precession angle                    rad    (00 only)
%    epsa        - cannonical precession angle                    rad    (00 only)
%    chia        - cannonical precession angle                    rad    (00 only)
%    precess     - matrix converting from "mod" to gcrf
%
%  locals        :
%    zeta        - precession angle                               rad
%    z           - precession angle                               rad
%    theta       - precession angle                               rad
%    oblo        - obliquity value at j2000 epoch                 rad
%
%  coupling      :
%    none        -
%
%  references    :
%    vallado       2022, 219, 227-229
%
% [prec,psia,wa,ea,xa] = precess ( ttt, opt );
% ----------------------------------------------------------------------------

function [prec,psia,wa,ea,xa] = precess ( ttt, opt )
    sethelp;

    % " to rad
    convrt = pi / (180.0*3600.0);
    ttt2= ttt * ttt;
    ttt3= ttt2 * ttt;

    prec = eye(3);

    % ------------------- fk4 b1950 precession angles --------------------
    if (opt == '50')

        % ---- Seidelmann pg 107
        % for these calls, ttt will come in with the current jd
        %            t1= 0.0; %(ttt - 2433282.42345905)/365242.198782;  % set start as B1850, 0.0 to simplify
        %            t2= (ttt - 2396758.203)/365242.198782;  % uses B1850
        %            ttt = t2 - t1;
        %            ttt2 = ttt * ttt;
        %            ttt3 = ttt * ttt2;
        %            fprintf(1,'50prec %15.9f  \n',ttt);
        % exp supp 61 pg 38
        %            psia = 50.3708 + 0.0050 * ttt;
        %            wa = 0.0;
        %            ea = 0.0;
        %            xa = 0.1247 - 0.0188 * ttt;
        %            zeta =  (23035.545 + 139.720*t1 + 0.060 *t1*t1)*ttt + (30.240-0.270*t1)*ttt2 + 17.995*ttt3; % "
        %            theta=  (20051.12 - 85.29*t1 - 0.37 *t1*t1)*ttt + (-42.65-0.37*t1)*ttt2 - 41.80*ttt3;
        %            z    =  (23035.545 + 139.720*t1 + 0.060 *t1*t1)*ttt + (109.480+0.390*t1)*ttt2 + 18.325*ttt3;

        % ---- Newcomb Exp Supp 61 approach, but see GTDS pg 3-17
        % Exp Supp 61 says use 1900? but gtds says use 1950.
        % for these calls, ttt will come in with the current jd
        t1= 0.0; %(ttt - 2415020.31352)/36524.2198782;  % set start as B1900, 0.0 to simplify
        %            t2= (ttt - 2415020.31352)/36525;  % uses B1900
        t2= (ttt - 2433282.42345905)/36525;  % uses B1950
        %            ttt = t2 - t1;
        %            ttt2 = ttt * ttt;
        %            ttt3 = ttt * ttt2;
        % exp supp 61 pg 38
        psia = 50.3708 + 0.0050 * ttt;
        wa = 0.0; % not sure which one is which...
        ea = 84428.26 -   46.845*ttt - 0.00059*ttt2 + 0.00181*ttt3;
        xa = 0.1247 - 0.0188 * ttt;
        fprintf(1,'50prec %15.9f  \n',ttt);
        % seems like Exp supp 61 is ok with 1900 as epoch, and Seidlemann is ok with listed measr,
        %            zeta =  (2304.25 + 1.396*t1)*ttt + 0.302*ttt2 + 0.018*ttt3; % "
        %            theta=  (2004.682 - 0.853*t1)*ttt -0.426*ttt2 - 0.042*ttt3;
        %            z    =  (2304.25 + 1.396*t1)*ttt + 1.093*ttt2 + 0.018*ttt3;
        % GTDS pg 3-17 using days from 1950 - avoids long rpecession
        % constants...
        zeta =  2304.9969*ttt + 0.302*ttt2 + 0.01808*ttt3; % "
        theta=  2004.2980*ttt -0.425936*ttt2 - 0.0416*ttt3;
        z    =  2304.9969*ttt + 1.092999*ttt2 + 0.0192*ttt3;


        % tp-008 36-45
        % ttt is tropical centruies from 1950 36524.22 days
        prec(1,1) =  1.0 - 2.9696e-4 * ttt2 - 1.3e-7 * ttt3;
        prec(1,2) =  2.234941e-2 * ttt + 6.76e-6 * ttt2 - 2.21e-6 * ttt3;
        prec(1,3) =  9.7169e-3 * ttt - 2.07e-6 * ttt2 - 9.6e-7 * ttt3;
        prec(2,1) =  -prec(1,2);
        prec(2,2) =  1.0 - 2.4975e-4 * ttt2 - 1.5e-7 * ttt3;
        prec(2,3) =  - 1.0858e-4 * ttt2;
        prec(3,1) =  -prec(1,3);
        prec(3,2) =  prec(2,3);
        prec(3,3) =  1.0 - 4.721e-5 * ttt2;


        % pass these back out for testing
        psia = zeta;
        wa = theta;
        ea = z;
        % ------------------- iau 76 precession angles --------------------
    else
        if (opt == '80')
            %     fprintf(1,'80prec %15.9f  \n',ttt);
            psia =             5038.7784*ttt - 1.07259*ttt2 - 0.001147*ttt3; % "
            wa   = 84381.448                 + 0.05127*ttt2 - 0.007726*ttt3;
            ea   = 84381.448 -   46.8150*ttt - 0.00059*ttt2 + 0.001813*ttt3;
            xa   =               10.5526*ttt - 2.38064*ttt2 - 0.001125*ttt3;

            zeta =             2306.2181*ttt + 0.30188*ttt2 + 0.017998*ttt3; % "
            theta=             2004.3109*ttt - 0.42665*ttt2 - 0.041833*ttt3;
            z    =             2306.2181*ttt + 1.09468*ttt2 + 0.018203*ttt3;
            % ------------------ iau 06 precession angles -------------------
        else
            oblo =  84381.406; % "
            psia =  (((( -0.0000000951 * ttt + 0.000132851 ) * ttt - 0.00114045 ) * ttt - 1.0790069 ) * ttt + 5038.481507 ) * ttt; % "
            wa   =  ((((  0.0000003337 * ttt - 0.000000467 ) * ttt - 0.00772503 ) * ttt + 0.0512623 ) * ttt -    0.025754 ) * ttt + oblo;
            ea   =  (((( -0.0000000434 * ttt - 0.000000576 ) * ttt + 0.00200340 ) * ttt - 0.0001831 ) * ttt -   46.836769 ) * ttt + oblo;
            xa   =  (((( -0.0000000560 * ttt + 0.000170663 ) * ttt - 0.00121197 ) * ttt - 2.3814292 ) * ttt +   10.556403 ) * ttt;

            zeta =  (((( -0.0000003173 * ttt - 0.000005971 ) * ttt + 0.01801828 ) * ttt + 0.2988499 ) * ttt + 2306.083227 ) * ttt + 2.650545; % "
            theta=  (((( -0.0000001274 * ttt - 0.000007089 ) * ttt - 0.04182264 ) * ttt - 0.4294934 ) * ttt + 2004.191903 ) * ttt;
            z    =  ((((  0.0000002904 * ttt - 0.000028596 ) * ttt + 0.01826837 ) * ttt + 1.0927348 ) * ttt + 2306.077181 ) * ttt - 2.650545;

            % fukishima williams approach
            %gamb
            %     zeta = ((((0.0000000260 * ttt - 0.000002788) * ttt - 0.00031238) * ttt + 0.4932044) * ttt + 10.556378) * ttt + -0.052928; % "
            %phib
            %     theta = ((((-0.0000000176 * ttt - 0.000000440) * ttt + 0.00053289) * ttt + 0.0511268) * ttt + -46.811016) * ttt + 84381.412819;
            %psib
            %     z = ((((-0.0000000148 * ttt - 0.000026452) * ttt - 0.00018522) * ttt + 1.5584175) * ttt + 5038.481484) * ttt - 0.041775;

        end
    end

    % convert units to rad
    psia = psia  * convrt; % rad
    wa   = wa    * convrt;
    ea   = ea    * convrt;
    xa   = xa    * convrt;

    zeta = zeta  * convrt;
    theta= theta * convrt;
    z    = z     * convrt;
    %iauhelp='y';
    %fprintf(1,'pr %11.7f  %11.7f  %11.7fdeg \n',zeta*180/pi,theta*180/pi,z*180/pi );
    if (iauhelp == 'y')
        fprintf(1,'pr %11.7f  %11.7f  %11.7f %11.7fdeg \n',psia*180/pi,wa*180/pi,ea*180/pi,xa*180/pi );
        fprintf(1,'pr %11.7f  %11.7f  %11.7fdeg \n',zeta*180/pi,theta*180/pi,z*180/pi );
    end

    if ( strcmp(opt,'80')==1 || strcmp(opt,'06') == 1)
        coszeta  = cos(zeta);
        sinzeta  = sin(zeta);
        costheta = cos(theta);
        sintheta = sin(theta);
        cosz     = cos(z);
        sinz     = sin(z);

        % ----------------- form matrix  mod to j2000 -----------------
        prec(1,1) =  coszeta * costheta * cosz - sinzeta * sinz;
        prec(1,2) =  coszeta * costheta * sinz + sinzeta * cosz;
        prec(1,3) =  coszeta * sintheta;
        prec(2,1) = -sinzeta * costheta * cosz - coszeta * sinz;
        prec(2,2) = -sinzeta * costheta * sinz + coszeta * cosz;
        prec(2,3) = -sinzeta * sintheta;
        prec(3,1) = -sintheta * cosz;
        prec(3,2) = -sintheta * sinz;
        prec(3,3) =  costheta;

        % ----------------- do rotations instead ----------------------
        %             fprintf(1,'prec matrix \n');
        %             prec
        %             fprintf(1,'prec rotations z \n');
        %             p1 = rot3mat( z );
        %             p2 = rot2mat( -theta );
        %             p3 = rot3mat( zeta );
        %             prec1 = p3*p2*p1
        %
        %             fprintf(1,'prec rotations w \n');
        %             a4  = rot3mat(-xa);
        %             a5  = rot1mat(wa);
        %             a6  = rot3mat(psia);
        %             a7  = rot1mat(-ea);
        %             prec2 = a7*a6*a5*a4
    else  %if %(strcmp(opt,'50') ~= 1)
        %             oblo = oblo * convrt; % " to rad
        %             a4  = rot3mat(-xa);
        %             a5  = rot1mat(wa);
        %             a6  = rot3mat(psia);
        %             a7  = rot1mat(-oblo);
        %             prec3 = a7*a6*a5*a4
    end

end