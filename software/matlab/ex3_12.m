% ------------------------------------------------------------------------------
%
%                              Ex3_12.m
%
%  this file demonstrates example 3-12.
%
%                          companion code for
%             fundamentals of astrodyanmics and applications
%                                 2022
%                            by david vallado
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
% ------------------------------------------------------------------------------
    
    year = 2001;
    mon = 3;
    day = 18;
    hr  = 12;
    min = 14;
    sec = 0.0;
    [jd, jdfrac] = jday(year, mon, day, hr, min, sec);
    fprintf(1,'input jd %14.7f %11.8f  \n\n', jd, jdfrac );

     % check jdfrac for multiple days
     if abs(jdfrac) >= 1.0 
         jd = jd + floor(jdfrac);
         jdfrac = jdfrac - floor(jdfrac);
     end

     % check for fraction of a day included in the jd
	 dt = jd - floor(jd) - 0.5;
     if (abs(dt) > 0.00000001)
         jd = jd - dt;  
         jdfrac = jdfrac + dt;
     end

     % ----------------- find year and days of the year ---------------
     temp   = jd - 2415019.5; 
     tu     = temp / 365.25;
     year   = 1900 + floor( tu );
     leapyrs= floor( ( year-1901 )*0.25 );
     days   = floor(temp - ((year-1900)*365.0 + leapyrs ));

     % ------------ check for case of beginning of a year -------------
     if days + jdfrac < 1.0
         year   = year - 1;
         leapyrs= floor( ( year-1901 )*0.25 );
         days   = floor(temp - ((year-1900)*365.0 + leapyrs ));
     end

    fprintf(1,'year %6i  \n',year);
    fprintf(1,'mon  %3i  \n',mon);
    fprintf(1,'days  %3i  \n\n',days);
    
    daysf = days + hr/24.0 + min/1440.0 + sec/86400.0;
    fprintf(1,'days  %11.7f  \n',daysf);
    
    utsec = (daysf-days) * 86400.0;
    
    temp   = utsec   / 3600.0;
    hr  = fix( temp  );
    min = fix( (temp - hr)*60.0 );
    sec = (temp - hr - min/60.0 ) * 3600.0;
        
    fprintf(1,' hr min sec %4i  %4i  %8.6f \n',hr, min, sec);
    

