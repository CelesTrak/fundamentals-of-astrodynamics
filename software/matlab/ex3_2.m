    %     -----------------------------------------------------------------
    %
    %                              Ex3_2.m
    %
    %  this file demonstrates example 3-2.
    %
    %                          companion code for
    %             fundamentals of astrodyanmics and applications
    %                                 2013
    %                            by david vallado
    %
    %     (h)               email davallado@gmail.com
    %     (w) 719-573-2600, email dvallado@agi.com
    %
    %     *****************************************************************
    %
    %  current :
    %            16 feb 19  david vallado
    %                         update for new constants
    %  changes :
    %            13 feb 07  david vallado
    %                         original baseline
    %
    %     *****************************************************************

    constastro;
    constmath;

    deg =  -7;
    min =  -54;
    sec = -23.886;
    fprintf(1,'deg %2i min %2i sec %8.6f\n',deg,min,sec );
    [latgd] = dms2rad( deg,min,sec );
    fprintf(1,'dms = %11.7f rad\n',latgd*rad );

    deg = 345;
    min =  35;
    sec =  51.000;
    fprintf(1,'deg %2i min %2i sec %8.6f\n',deg,min,sec );
    [lon] = dms2rad( deg,min,sec );
    fprintf(1,'dms = %11.7f rad\n',lon*rad );

    alt = 0.056;

    % -------------------------  implementation   -----------------
    sinlat      = sin( latgd );
    earthrate(1)= 0.0;
    earthrate(2)= 0.0;
    earthrate(3)= earthrot;

    % ------  find rdel and rk components of site vector  ---------
    cearth= re / sqrt( 1.0 - ( eccearthsqrd*sinlat*sinlat ) )
    rdel  = ( cearth + alt )*cos( latgd );
    rk    = ( (1.0-eccearthsqrd)*cearth + alt )*sinlat;

    (1.0-eccearthsqrd)*cearth

    % ---------------  find site position vector  -----------------
    rs(1) = rdel * cos( lon );
    rs(2) = rdel * sin( lon );
    rs(3) = rk;
    rs = rs';

    rsmag = mag(rs);
    fprintf(1,'site gd %16.9f %16.9f %16.9f   %16.9f \n',rs, rsmag );

    latgc = gd2gc(latgd);

    r(1)= rsmag*cos(latgc)*cos(lon);
    r(2)= rsmag*cos(latgc)*sin(lon);
    r(3)= rsmag*sin(latgc);
    r = r';

    fprintf(1,'site gc %16.9f %16.9f %16.9f  %16.9f  \n',r, latgc*rad );

