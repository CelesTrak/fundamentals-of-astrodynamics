% ------------------------------------------------------------------------------
%
%                              Ex2_4.m
%
%  this file demonstrates example 2-4.
%
%                          companion code for
%             fundamentals of astrodyanmics and applications
%                                 2022
%                            by david vallado
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
% ------------------------------------------------------------------------------

    constastro;

    fprintf(1,'\n-------- kepler  ex 2-4, pg 102 --------- \n' );

    % initial vectors in km and km/s
    ro = [ 1131.340  -2282.343  6672.423];
    vo = [ -5.64305  4.30333  2.42879 ];
    fprintf(1,'input: \n' );
    fprintf(1,'ro %16.8f %16.8f %16.8f km \n',ro );
    fprintf(1,'vo %16.8f %16.8f %16.8f km/s \n',vo );

    % convert 40 minutes to seconds
    dtsec = 40.0*60.0;
    fprintf(1,'dt %16.8f sec \n',dtsec );
    fprintf(1,'intermediate values: \n' );

    [r1,v1] =  kepler ( ro,vo, dtsec );

    % answer in km and km/s
    fprintf(1,'output: \n' );
    fprintf(1,'r1 %16.8f %16.8f %16.8f er \n',r1/re );
    fprintf(1,'r1 %16.8f %16.8f %16.8f km \n',r1 );
    fprintf(1,'v1 %16.8f %16.8f %16.8f er/tu \n',v1/velkmps );
    fprintf(1,'v1 %16.8f %16.8f %16.8f km/s \n\n\n',v1 );


   