% ------------------------------------------------------------------------------
%
%                              Ex4_1.m
%
%  this file demonstrates example 4-1.
%
%                          companion code for
%             fundamentals of astrodyanmics and applications
%                                 2022
%                            by david vallado
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
% ------------------------------------------------------------------------------

    constmath;
    
    % Helio for Neptune from the Sun
    paral1 = asin(1.0 / 29.664361);
    fprintf(1,'paralax  %11.9f %11.9f arcsec \n', paral1*rad, paral1*rad*3600.0 );
    
    
    % try exact method (neptune to earth)
    au = 149597870.7;
    rnep = [-1757460712.2509    3757470099.1416    1576777174.1537]; %km
   
    %  neptune to sun vector from JPL epehmerides 1994, 5/14
    rearth = [ -1666604612.0985    3868340828.5807    1624846829.1305];  %km
    paral2 = acos( dot(rearth, rnep) / ( mag(rearth)*mag(rnep) ) );
    fprintf(1,'paralax  %11.9f %11.9f arcsec \n', paral2*rad, paral2*rad*3600.0 );
    
    paral1/paral2 * 100
    
    % for Neptune to Earth
    paral1 = asin(6378. /(149597870*29.664361));
    fprintf(1,'paralax  %11.9f %11.9f arcsec \n', paral1*rad, paral1*rad*3600.0 );
    
    % try exact method
    % neptune to Earth
    rnep = [-1757460712.2509    3757470099.1416    1576777174.1537];  % km
    % neptune to site
    rsite = [6378.137  0.0  0.0];
    rearth  = rnep + rsite;
    paral2 = acos( dot(rnep, rearth) / ( mag(rnep)*mag(rearth) ) );
    fprintf(1,'paralax  %11.9f %11.9f arcsec \n', paral2*rad, paral2*rad*3600.0 );
    
    paral1/paral2 * 100
    
    
    
    