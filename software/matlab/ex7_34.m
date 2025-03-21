% ------------------------------------------------------------------------------
%
%                              Ex7_34.m
%
%  this file demonstrates example 7-34.
%
%                          companion code for
%             fundamentals of astrodyanmics and applications
%                                 2022
%                            by david vallado
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
% ------------------------------------------------------------------------------

      fprintf(1,'-------------------- problem ex 7-3 \n');
      r1 = [0.0 0.0 6378.137];
      r2 = [0.0 -4464.696 -5102.509];
      r3 = [0.0 5740.323 3189.068];
      fprintf(1,' r1 %15.10f  %15.10f  %15.10f\n',r1 );
      fprintf(1,' r2 %15.10f  %15.10f  %15.10f\n',r2 );
      fprintf(1,' r3 %15.10f  %15.10f  %15.10f\n',r3 );

      [v2g, theta,theta1,copa, errorg] = gibbs( r1,r2,r3);

      fprintf(1,' v2g %15.9f   %15.9f   %15.9f km/s \n\n\n',v2g );
      rad = 180.0 / pi;
      fprintf(1,' theta %11.6f  theta1  %11.6f  copa %11.6f \n',theta*rad, theta1*rad,copa*rad);


      fprintf(1,'-------------------- problem ex 7-4 \n');
      r1 = [3419.85564  6019.82602  2784.60022];
      r2 = [2935.91195  6326.18324  2660.59584];
      r3 = [2434.95202  6597.38674  2521.52311];
      fprintf(1,' r1 %15.10f  %15.10f  %15.10f\n',r1 );
      fprintf(1,' r2 %15.10f  %15.10f  %15.10f\n',r2 );
      fprintf(1,' r3 %15.10f  %15.10f  %15.10f\n',r3 );

      % just needs to be in secs - so precision isn't lost
      jd1 = 0.0/86400.0;
      jd2 = (60.0 + 16.48)/86400.0;
      jd3 = (120.0 + 33.04)/86400.0;

      % use t1, t2, t3 as secs for greater accuracy
      [v2h, theta,theta1,copa, errorh ] =  hgibbs ( r1,r2,r3,jd1,jd2,jd3 );

      fprintf(1,' v2h %15.9f   %15.9f   %15.9f km/s \n\n\n',v2h );
      fprintf(1,' theta %11.6f  theta1  %11.6f  copa %11.6f \n',theta*rad, theta1*rad,copa*rad);
pause;
%
% 
%  this section tests the gibbs routine
%

      fprintf(1,'-------------------- compare gibbs and hgibbs \n');
        constmath;
        constastro;

        r(1)= 1.0;
        r(2)= 0.10;
        r(3)= 0.9760;
        r = r * re;
        magr = mag(r);

        v(1)= 0.3;
        v(2)= 0.10;
        v(3)= 0.2760;
        v = v * velkmps;
        magv = mag(v);

        p    = 1.61*re;
        ecc    = 0.0;
        incl = 34.0/rad;
        omega= 0.0;
        argp = 0.0;
        nu  = 0.0;
        arglat  = 0.0;
        truelon = 0.0;
        lonper  = 0.0;
        [r,v] = coe2rv(p,ecc,incl,omega,argp,nu,arglat,truelon,lonper);

        % use kepler propagation
%        r = [-4550.4256  2220.1946  5090.3547];
%        v = [-4.975006  -5.237109  -2.11811];
%        for i = 0:50
%            [r1,v1,error] =  kepler  ( r,v, i*600.0 );
%            fprintf(1,'%4i x %11.7f  %11.7f  %11.7f  %11.7f  %11.7f  %11.7f \n',i,r1,v1 );
%          end;

        % test gibbs hgibbs accuracies against two body orbits
        i= 1;
        t2= 0.0; % min

        magr = mag( r );
        magv = mag( v );
        rdotv= dot( r,v );

        % -------------  find sme, alpha, and a  ------------------
        sme= ( (magv^2)*0.5  ) - ( mu /magr );
        a= -mu / ( 2.0 *sme );
        period= twopi * sqrt( abs(a)^3.0/mu ); % sec

        fprintf(1,' ktr    dt sec   fraction       ang1          ang2   ');
        fprintf(1,'    coplanar        gibbs      err flg             hgibbs      errflg \n');
        while t2 < period * 0.45
            if t2 < period * 0.001
               t2= t2 + period*0.000005;
              else
                if t2 < period * 0.006
                   t2= t2 + period*0.001;
                  else
                    if t2 < period * 0.015
                        t2= t2 + period*0.005;
                      else
                        if t2 < period * 0.20
                            t2= t2 + period*0.01;
                          else
                            t2= t2 + period*0.05;
                          end;
                      end;
                  end;
               end;

            t1= 0.0;
            t2= t2;  % sec
            t3= 2.0*t2;


            [r1,v1] =  kepler  ( r,v, t1 );
            [r2,v2] =  kepler  ( r,v, t2 );
            [r3,v3] =  kepler  ( r,v, t3 );

            [v2g, theta,theta1,copa, errorg] = gibbs( r1,r2,r3);

            % just needs to be in days - so precision isn't lost
            jd1 = t1/86400.0;
            jd2 = t2/86400.0;
            jd3 = t3/86400.0;
            % use t1, t2, t3 as secs for greater accuracy
            [v2h, theta,theta1,copa, errorh ] =  hgibbs ( r1,r2,r3,jd1,jd2,jd3 );

%v2 = v2 / velkmps;
%v2g = v2g / velkmps;
%v2h = v2h / velkmps;

%            fprintf(1,'v2truth %14.7f  %14.7f  %14.7f  \n',v2(1),v2(2),v2(3) );
%            fprintf(1,'v2g     %14.7f  %14.7f  %14.7f  \n',v2g(1),v2g(2),v2g(3) );
%            fprintf(1,'v2h     %14.7f  %14.7f  %14.7f  \n',v2h(1),v2h(2),v2h(3) );

%            fprintf(1,'angles %14.7f �  %14.7f  %14.7f \n',theta*rad,theta1*rad,copa*rad );

            tempg = v2g - v2;
            temph = v2h - v2;
            magtempg = mag(tempg);
            magtemph = mag(temph);
%            fprintf(1,'dg     %14.7f  %14.7f  %14.7f  \n',tempg(1),tempg(2),tempg(3) );
%            fprintf(1,'dh     %14.7f  %14.7f  %14.7f  \n',temph(1),temph(2),temph(3) );

            fprintf(1,'%3i %11.3f %8.4f  %12.6f � %12.6f %10.5f dg %14.9f %12s dh %14.9f m/s %12s  \n', ...
                      i,t2,t2/period,theta*rad,theta1*rad,copa*rad,magtempg*1000,errorg,magtemph*1000,errorh );
            i = i + 1;
          end;

pause;

% 6 Apr 2004 19:02:28.3860     12709.670465397     9059.657416996      4937.914304833    -2.425788263     1.145863933     4.141391183
% 6 Apr 2004 19:03:28.3860     12562.052651693     9126.923951210      5185.576718081    -2.494671410     1.096292760     4.113797579
% 6 Apr 2004 19:04:28.3860     12410.326214883     9191.205376885      5431.543104340    -2.562738634     1.046363027     4.084858491
% 
% 6 Apr 2004 19:08:28.3860     11763.346541426     9418.074295602      6396.848305457    -2.826407207     0.843386744     3.955843131
% 6 Apr 2004 19:14:28.3860     10678.648538971     9665.704898623      7780.535182624    -3.193778683     0.530988666     3.723761880
       
     % closely spaced
      r1 = [12709.670465397     9059.657416996      4937.914304833];
      r2 = [12562.052651693     9126.923951210      5185.576718081];
      r3 = [12410.326214883     9191.205376885      5431.543104340];

      fprintf(1,'-------------------- problem close spaced circorbit \n');
      fprintf(1,' r1 %15.10f  %15.10f  %15.10f\n',r1 );
      fprintf(1,' r2 %15.10f  %15.10f  %15.10f\n',r2 );
      fprintf(1,' r3 %15.10f  %15.10f  %15.10f\n',r3 );

      [v2g, theta,theta1,copa, errorg] = gibbs( r1,r2,r3);
      fprintf(1,' v2g %15.9f   %15.9f   %15.9f km/s \n',v2g );
      fprintf(1,' theta %11.6f  theta1  %11.6f  copa %11.6f \n',theta*rad, theta1*rad,copa*rad);

      % just needs to be in secs - so precision isn't lost
      jd1 = 0.0;
      jd2 = 60.0/86400.0; %/86400.0;
      jd3 = 120.0/86400.0 ; %/86400.0;

      % use t1, t2, t3 as secs for greater accuracy
      [v2h, theta,theta1,copa, errorh ] =  hgibbs ( r1,r2,r3,jd1,jd2,jd3 );
      fprintf(1,' v2h %15.9f   %15.9f   %15.9f km/s \n',v2h );
      fprintf(1,' theta %11.6f  theta1  %11.6f  copa %11.6f \n',theta*rad, theta1*rad,copa*rad);
   
      % wider spaced 
      r1 = [12709.670465397     9059.657416996      4937.914304833 ];
      r2 = [11763.346541426     9418.074295602      6396.848305457];
      r3 = [10678.648538971     9665.704898623      7780.535182624];
      
      fprintf(1,'-------------------- problem wider spaced circorbit \n');
      fprintf(1,' r1 %15.10f  %15.10f  %15.10f\n',r1 );
      fprintf(1,' r2 %15.10f  %15.10f  %15.10f\n',r2 );
      fprintf(1,' r3 %15.10f  %15.10f  %15.10f\n',r3 );
      [v2g, theta,theta1,copa, errorg] = gibbs( r1,r2,r3);

      fprintf(1,' v2g %15.9f   %15.9f   %15.9f km/s \n',v2g );
      fprintf(1,' theta %11.6f  theta1  %11.6f  copa %11.6f \n',theta*rad, theta1*rad,copa*rad);

      % just needs to be in secs - so precision isn't lost
      jd1 = 0.0; %/86400.0;
      jd2 = 360.0/86400.0; %/86400.0;
      jd3 = 720.0/86400.0; %/86400.0;

      % use t1, t2, t3 as secs for greater accuracy
      [v2h, theta,theta1,copa, errorh ] =  hgibbs ( r1,r2,r3,jd1,jd2,jd3 );
      fprintf(1,' v2h %15.9f   %15.9f   %15.9f km/s \n',v2h );
      fprintf(1,' theta %11.6f  theta1  %11.6f  copa %11.6f \n',theta*rad, theta1*rad,copa*rad);

