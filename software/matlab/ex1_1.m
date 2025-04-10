% ------------------------------------------------------------------------------
%
%                              Ex1_1.m
%
%  this file demonstrates example 1-1.
%
%                          companion code for
%             fundamentals of astrodyanmics and applications
%                                 2022
%                            by david vallado
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
% ------------------------------------------------------------------------------

    fprintf(1,'\n-------- ex 1-1 --------- \n' );

    constmath;
    constastro;
    periodsid = 86164.090524;

    a = (mu*(periodsid/twopi)^2)^(1.0/3.0);

    fprintf(1,'a %16.8f km  %18.10f er \n',a,a/re );

    fprintf(1,'re %16.8f km  mu %18.10f km3/s2 \n', re, mu );
    fprintf(1,'periodsid     %16.10f \n',periodsid );
    fprintf(1,'tuday    %16.14f 1/tuday %18.12f  \n', tuday, 1.0/tuday );
    fprintf(1,'tudaysid %16.14f 1/tuday %18.12f  \n', tudaysid, 1.0/tudaysid );
    fprintf(1,'periodsid alt %16.10f \n',86400/1.002737909350795 );

