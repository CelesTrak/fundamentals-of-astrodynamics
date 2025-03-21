% ------------------------------------------------------------------------------
%
%                              Ex3_9.m
%
%  this file demonstrates example 3-9.
%
%                          companion code for
%             fundamentals of astrodyanmics and applications
%                                 2022
%                            by david vallado
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
% ------------------------------------------------------------------------------

    % -------- hms test
    hr = 15;
    min = 15;
    sec = 53.63;
    fprintf(1,'hr %4i ',hr);
    fprintf(1,'min %4i ',min);
    fprintf(1,'sec %8.6f \n',sec);

    [hms] = hms2rad( hr,min,sec );

    fprintf(1,'hms %11.7f %11.7f \n',hms, hms * 180.0/pi);

    [hr,min,sec] = rad2hms( hms );

    fprintf(1,' hr min sec %4i  %4i  %8.6f \n',hr, min, sec);

