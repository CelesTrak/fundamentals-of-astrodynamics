    %
    %  testall
%
% test all the matlab routines
% ths allows comparison withother langugaes (c#, c++, etc)
%
%

%todo - all indicies have to be +1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


%exampleEnumType = struct('zero',0, 'one',1, 'two',2);


function fid = testall

    testnum = -10;

    testnum = input('please enter -10 for all, or test number \n');

    if (testnum == -10)
        optstart = 1;  % 1
        optstop = 80; % 102
    else
        optstart = testnum;
        optstop = testnum;
    end

    fclose('all'); % clear out all the fids except 0, 1, 2
    directory = 'D:\Codes\LIBRARY\matlab\';
    fid = fopen(strcat(directory,'testallm.out'), 'wt');

    for opt = optstart:optstop
        fprintf(fid,'\n\n=================================== Case %i =======================================\n', opt);
        fprintf(1,'\n\n=================================== Case %i =======================================\n', opt);

        switch opt
            case 1
                testvecouter(fid);
            case 2
                testmatadd(fid);
            case 3
                testmatsub(fid);
            case 4
                testmatmult(fid);
            case 5
                testmattransr(fid);
            case 6
                testmattransrx(fid);
            case 7
                testmatinverse(fid);
            case 8
                testdeterminant(fid);
            case 9
                testcholesky(fid);
            case 10
                testposvelcov2pts(fid);
            case 11
                testposcov2pts(fid);
            case 12
                testremakecovpv(fid);
            case 13
                testremakecovp(fid);
            case 14
                testmatequal(fid);
            case 15
                testmatscale(fid);
            case 16
                testunit(fid);
            case 17
                testmag(fid);
            case 18
                testcross(fid);
            case 19
                testdot(fid);
            case 20
                testangle(fid);
            case 21
                testasinh(fid);
            case 22
                testcot(fid);
            case 23
                testacosh(fid);
            case 24
                testaddvec(fid);
            case 25
                testPercentile(fid);
            case 26
                testrot1(fid);
            case 27
                testrot2(fid);
            case 28
                testrot3(fid);
            case 29
                testfactorial(fid);
            case 30
                testcubicspl(fid);
            case 31
                testcubic(fid);
            case 32
                testcubicinterp(fid);
            case 33
                testquadratic(fid);
            case 34
                testconvertMonth(fid);
            case 35
                testjday(fid);
            case 36
                testdays2mdhms(fid);
            case 37
                testinvjday(fid);
            case 38
                % tests eop, spw, and fk5 iau80
                testiau80in(fid);
            case 39
                testfundarg(fid);
            case 40
                testprecess(fid);
            case 41
                testnutation(fid);
            case 42
                testnutationqmod(fid);
            case 43
                testsidereal(fid);
            case 44
                testpolarm(fid);
            case 45
                testgstime(fid);
            case 46
                testlstime(fid);
            case 47
                testhms_sec(fid);
            case 48
                testhms_ut(fid);
            case 49
                testhms_rad(fid);
            case 50
                testdms_rad(fid);
            case 51
                testeci_ecef(fid);
            case 52
                testtod2ecef(fid);
            case 53
                testteme_ecef(fid);
            case 54
                testteme_eci(fid);
            case 55
                testqmod2ecef(fid);
            case 56
                testcsm2efg(fid);
            case 57
                testrv_elatlon(fid);
            case 58
                testrv2radec(fid);
            case 59
                testrv_razel(fid);
            case 60
                testrv_tradec(fid);
            case 61
                testrvsez_razel(fid);
            case 62
                testrv2rsw(fid);
            case 63
                testrv2pqw(fid);
            case 64
                testrv2coe(fid);
            case 65
                testcoe2rv(fid);
            case 66
                testlon2nu(fid);
            case 67
                testnewtonmx(fid);
            case 68
                testnewtonm(fid);
            case 69
                testnewtonnu(fid);
            case 70
                testfindc2c3(fid);
            case 71
                testfindfandg(fid);
            case 72
                testcheckhitearth(fid);
                testcheckhitearthc(fid);
            case 73
                testgibbs(fid);
                testhgibbs(fid);
            case 74
                % in sln directory, testall-Angles.out
                testangles(fid);
                % D\faabook\current\excel\testgeo.out for ch9 plot
                testgeo(fid);
            case 75
                testlambertumins(fid);
                testlambertminT(fid);
            case 76
                testlambhodograph(fid);
            case 77
                testlambertbattin(fid);
            case 78
                testeq2rv(fid);
            case 79
                testrv2eq(fid);
            case 80
                directoryx = 'd\codes\library\matlab\';
                % book example, simple
                testlambertuniv(fid);
                fprintf(fid,'lambert envelope test case results written to %s testall.out ', directoryx);

                % envelope testing
                testAllLamb(fid);
                fprintf(fid,'lambert envelope test case results written to %s testall.out ', directoryx);

            case 81
                testradecgeo2azel(fid);
            case 82
                testijk2ll(fid);
            case 83
                testgd2gc(fid);
            case 84
                testsite(fid);
            case 85
                testsun(fid);
            case 86
                testmoon(fid);
            case 87
                testkepler(fid);
            case 88
                % mean
                testcovct2clmean(fid);
            case 89
                % true
                testcovct2cltrue(fid);
            case 90
                testhill(fid);
            case 91
            case 92
                testcovct2rsw(fid);
                testcovct2ntw(fid);
                testcovcl2eq('truea', fid);
                testcovcl2eq('truen', fid);
                testcovcl2eq('meana', fid);
                testcovcl2eq('meann', fid);
                testcovct2eq('truea', fid);
                testcovct2eq('truen', fid);
                testcovct2eq('meana', fid);
                testcovct2eq('meann', fid);
                testcovct2fl('latlon', fid);
                testcovct2fl('radec', fid);
            case 93
            case 94
            case 95
            case 96
                testcovct2rsw(fid);
            case 97
                testcovct2ntw(fid);
            case 98
                testsunmoonjpl(fid);
            case 99
                testkp2ap(fid);
            case 100
                testproporbit(fid);
                directoryy = 'd\codes\library\matlab\';
                fprintf(fid,'testproporbit (legendre) results written to %s legpoly.out ', directoryy);
                fprintf(fid,'testproporbit (legendre) results written to %s legendreAcc.out ', directoryy);
            case 101
                %testsemianaly(fid);
            case 102
                testazel2radec(fid);
        end

    end % for
    
    fclose('all');
end  % runtests


function testvecouter(fid)

    vec1 = [ 2.3; 4.7; -1.6 ];
    vec2 = [ 0.3; -0.7; 6.0 ];

    mat1 = vec1 * vec2';

    fprintf(fid,'vecout = %15.11f %15.11f %15.11f\n', mat1(1, 1), mat1(1, 2), mat1(1, 3));
    fprintf(fid,'vecout = %15.11f %15.11f %15.11f\n', mat1(2, 1), mat1(2, 2), mat1(2, 3));
    fprintf(fid,'vecout = %15.11f %15.11f %15.11f\n', mat1(3, 1), mat1(3, 2), mat1(3, 3));
end


function testmatadd(fid)
    mat1 = [ 1.0,   2.0,   3.0 ;
        -1.1,   0.5,   2.0 ;
        -2.00,  4.00,  7.0 ];
    mat2 = [ 1.0,  1.4, 1.8 ;
        0.0,  2.6, -0.6 ;
        1.9,  0.1, 7.1  ];

    mat3 = mat1 + mat2;

    fprintf(fid,'matadd = %15.11f %15.11f %15.11f\n', mat3(1, 1), mat3(1, 2), mat3(1, 3));
    fprintf(fid,'matadd = %15.11f %15.11f %15.11f\n', mat3(2, 1), mat3(2, 2), mat3(2, 3));
    fprintf(fid,'matadd = %15.11f %15.11f %15.11f\n', mat3(3, 1), mat3(3, 2), mat3(3, 3));
end


function testmatsub(fid)
    mat1 = [ 1.0,   2.0,   3.0 ;
        -1.1,   0.5,   2.0 ;
        -2.00,  4.00,  7.0 ];
    mat2 = [ 1.0,  1.4, 1.8 ;
        0.0,  2.6, -0.6 ;
        1.9,  0.1, 7.1  ];

    mat3 = mat1 - mat2;

    fprintf(fid,'matsub = %15.11f %15.11f %15.11f\n', mat3(1, 1), mat3(1, 2), mat3(1, 3));
    fprintf(fid,'matsub = %15.11f %15.11f %15.11f\n', mat3(2, 1), mat3(2, 2), mat3(2, 3));
    fprintf(fid,'matsub = %15.11f %15.11f %15.11f\n', mat3(3, 1), mat3(3, 2), mat3(3, 3));
end


function testmatmult(fid)
    mat1 = [ 1.0,   2.0,   3.0 ;
        -1.1,   0.5,   2.0 ;
        -2.00,  4.00,  7.0 ];
    mat2 = [  1.0,  1.4 ;
        0.0,  2.6 ;
        1.9,  0.1  ];
    mat3 = mat1 * mat2;

    fprintf(fid,'matmult = %15.11f %15.11f\n', mat3(1, 1), mat3(1, 2));
    fprintf(fid,'matmult = %15.11f %15.11f\n', mat3(2, 1), mat3(2, 2));
    fprintf(fid,'matmult = %15.11f %15.11f\n', mat3(3, 1), mat3(3, 2));
end


function testmattransr(fid)
    mat1 = [ 1.0,   2.0,   3.0 ;
        -1.1,   0.5,   2.0 ;
        -2.00,  4.00,  7.0 ];
    [mat3] = mat1';

    fprintf(fid,'mattrans = %15.11f %15.11f %15.11f\n', mat3(1, 1), mat3(1, 2), mat3(1, 3));
    fprintf(fid,'mattrans = %15.11f %15.11f %15.11f\n', mat3(2, 1), mat3(2, 2), mat3(2, 3));
    fprintf(fid,'mattrans = %15.11f %15.11f %15.11f\n', mat3(3, 1), mat3(3, 2), mat3(3, 3));
end

function testmattransrx(fid)
    mat1 = [ 1.0,   2.0,   3.0 ;
        -1.1,   0.5,   2.0 ;
        -2.00,  4.00,  7.0 ];

    mat3 = mat1';
end

function testmatinverse(fid)
    % enter by COL!!!!!!!!!!!!!!
    mat1 = [ 3, 5, 6 ;  2, 0, 3  ; 1, 2, 8  ];
    matinv = inv(mat1);

    fprintf(fid,'matinv = %15.11f %15.11f %15.11f\n', matinv(1, 1), matinv(1, 2), matinv(1, 3));
    fprintf(fid,'matinv = %15.11f %15.11f %15.11f\n', matinv(2, 1), matinv(2, 2), matinv(2, 3));
    fprintf(fid,'matinv = %15.11f %15.11f %15.11f\n', matinv(3, 1), matinv(3, 2), matinv(3, 3));

    %Results: test before
    % 0.1016949    0.4745763 - 0.2542373
    % 0.2203390 - 0.3050847 - 0.0508475
    %- 0.0677966    0.0169492    0.1694915

    mat1 = [ 1, 3, 3  ; 1, 4, 3 ;  1, 3, 4  ];
    matinv = inv(mat1);

    fprintf(fid,'matinv = %15.11f %15.11f %15.11f\n', matinv(1, 1), matinv(1, 2), matinv(1, 3));
    fprintf(fid,'matinv = %15.11f %15.11f %15.11f\n', matinv(2, 1), matinv(2, 2), matinv(2, 3));
    fprintf(fid,'matinv = %15.11f %15.11f %15.11f\n', matinv(3, 1), matinv(3, 2), matinv(3, 3));


    ata = [ 264603537.493561, 206266447.729262, 274546062925.826, -282848493891885, 362835957483807, -4.3758299682612E+17 ;...
        206266447.729262, 160790924.64848, 214016946538.904, -220488942186083, 282841585443473, -3.41109159752805E+17 ;...
        274546062925.826, 214016946538.904, 284862180536440, -2.93476576836794E+17, 3.76469583735348E+17, -4.54025256502168E+20;...
        -282848493891885, -220488942186083, -2.93476576836794E+17, 3.02351477439543E+20, -3.87854240635812E+20, 4.67755241586584E+23;...
        362835957483807, 282841585443473, 3.76469583735348E+17, -3.87854240635812E+20, 4.97536553328938E+20, -6.00032966815125E+23;...
        -4.3758299682612E+17, -3.41109159752805E+17, -4.54025256502168E+20, 4.67755241586584E+23, -6.00032966815125E+23, 7.23644441510866E+26 ];

    matinv = inv(ata);

    fprintf(fid,'matinv = %15.11f %15.11f %15.11f %15.11f %15.11f %15.11f\n', matinv(1, 1), matinv(1, 2), matinv(1, 3), matinv(1, 4), matinv(1, 5), matinv(1, 6));
    fprintf(fid,'matinv = %15.11f %15.11f %15.11f %15.11f %15.11f %15.11f\n', matinv(2, 1), matinv(2, 2), matinv(2, 3), matinv(2, 4), matinv(2, 5), matinv(2, 6));
    fprintf(fid,'matinv = %15.11f %15.11f %15.11f %15.11f %15.11f %15.11f\n', matinv(3, 1), matinv(3, 2), matinv(3, 3), matinv(3, 4), matinv(3, 5), matinv(3, 6));
    fprintf(fid,'matinv = %15.11f %15.11f %15.11f %15.11f %15.11f %15.11f\n', matinv(4, 1), matinv(4, 2), matinv(4, 3), matinv(4, 4), matinv(4, 5), matinv(4, 6));
    fprintf(fid,'matinv = %15.11f %15.11f %15.11f %15.11f %15.11f %15.11f\n', matinv(5, 1), matinv(5, 2), matinv(5, 3), matinv(5, 4), matinv(5, 5), matinv(5, 6));
    fprintf(fid,'matinv = %15.11f %15.11f %15.11f %15.11f %15.11f %15.11f\n', matinv(6, 1), matinv(6, 2), matinv(6, 3), matinv(6, 4), matinv(6, 5), matinv(6, 6));

    mat3 = ata * matinv;

    fprintf(fid,'mat3 = %15.11f %15.11f %15.11f %15.11f %15.11f %15.11f\n', mat3(1, 1), mat3(1, 2), mat3(1, 3), mat3(1, 4), mat3(1, 5), mat3(1, 6));
    fprintf(fid,'mat3 = %15.11f %15.11f %15.11f %15.11f %15.11f %15.11f\n', mat3(2, 1), mat3(2, 2), mat3(2, 3), mat3(2, 4), mat3(2, 5), mat3(2, 6));
    fprintf(fid,'mat3 = %15.11f %15.11f %15.11f %15.11f %15.11f %15.11f\n', mat3(3, 1), mat3(3, 2), mat3(3, 3), mat3(3, 4), mat3(3, 5), mat3(3, 6));
    fprintf(fid,'mat3 = %15.11f %15.11f %15.11f %15.11f %15.11f %15.11f\n', mat3(4, 1), mat3(4, 2), mat3(4, 3), mat3(4, 4), mat3(4, 5), mat3(4, 6));
    fprintf(fid,'mat3 = %15.11f %15.11f %15.11f %15.11f %15.11f %15.11f\n', mat3(5, 1), mat3(5, 2), mat3(5, 3), mat3(5, 4), mat3(5, 5), mat3(5, 6));
    fprintf(fid,'mat3 = %15.11f %15.11f %15.11f %15.11f %15.11f %15.11f\n', mat3(6, 1), mat3(6, 2), mat3(6, 3), mat3(6, 4), mat3(6, 5), mat3(6, 6));
end


function testdeterminant(fid)
    order = 3;

    mat1(1, 1) = 6.0;
    mat1(1, 2) = 1.0;
    mat1(1, 3) = 1.0;
    mat1(2, 1) = 4.0;
    mat1(2, 2) = -2.0;
    mat1(2, 3) = 5.0;
    mat1(3, 1) = 2.0;
    mat1(3, 2) = 8.0;
    mat1(3, 3) = 7.0;

    deter = det(mat1);

    fprintf(fid,'det = %15.11f ansr -306/n', deter);
end

function testcholesky(fid)

    a = [1 0 1; 0 2 0; 1 0 3];
  
    matcho = chol(a);

    fprintf(fid,'matcho = %15.11f %15.11f %15.11f\n', matcho(1, 1), matcho(1, 2), matcho(1, 3));
    fprintf(fid,'matcho = %15.11f %15.11f %15.11f\n', matcho(2, 1), matcho(2, 2), matcho(2, 3));
    fprintf(fid,'matcho = %15.11f %15.11f %15.11f\n', matcho(3, 1), matcho(3, 2), matcho(3, 3));
end

function testposvelcov2pts(fid)
    reci = [5102.50895791863; 6123.01140072233; 6378.13692818659];
    veci = [-4.74322014817; 0.79053648924; 5.53375572723];

    cov2 = zeros(6,6);

    cov2(1, 1) = 12559.93762571587;
    cov2(1, 2) =  12101.56371305036;
    cov2(2, 1) = cov2(1, 2);
    cov2(1, 3) =  -440.3145384949657;
    cov2(3, 1) = cov2(1, 3);
    cov2(1, 4) = -0.8507401236198346;
    cov2(4, 1) = cov2(1, 4);
    cov2(1, 5) =  0.9383675791981778;
    cov2(5, 1) = cov2(1, 5);
    cov2(1, 6) =  -0.0318596430999798;
    cov2(6, 1) = cov2(1, 6);
    cov2(2, 2) = 12017.77368889201;
    cov2(2, 3) =  270.3798093532698;
    cov2(3, 2) = cov2(2, 3);
    cov2(2, 4) =  -0.8239662300032132;
    cov2(4, 2) = cov2(2, 4);
    cov2(2, 5) =  0.9321640899868708;
    cov2(5, 2) = cov2(2, 5);
    cov2(2, 6) =  -0.001327326827629336;
    cov2(6, 2) = cov2(2, 6);
    cov2(3, 3) = 4818.009967057008;
    cov2(3, 4) =  0.02033418761460195;
    cov2(4, 3) = cov2(3, 4) ;
    cov2(3, 5) =  0.03077663516695039;
    cov2(5, 3) = cov2(3, 5);
    cov2(3, 6) =  0.1977541628188323;
    cov2(6, 3) = cov2(3, 6);
    cov2(4, 4) = 5.774758755889862e-005;
    cov2(4, 5)  = -6.396031584925255e-005;
    cov2(5, 4) = cov2(4, 5);
    cov2(4, 6) =  1.079960679599204e-006;
    cov2(6, 4) = cov2(4, 6);
    cov2(5, 5) = 7.24599391355188e-005;
    cov2(5, 6) = 1.03146660433274e-006;
    cov2(6, 5) = cov2(5, 6);
    cov2(6, 6) = 1.870413627417302e-005;

    [sigmapts] = posvelcov2pts(reci, veci, cov2);

    fprintf(fid,'sigmapts = %15.11f %15.11f %15.11f %15.11f %15.11f %15.11f\n', sigmapts(1, 1), sigmapts(1, 2), sigmapts(1, 3), sigmapts(1, 4), sigmapts(1, 5), sigmapts(1, 6));
    fprintf(fid,'sigmapts = %15.11f %15.11f %15.11f %15.11f %15.11f %15.11f\n', sigmapts(2, 1), sigmapts(2, 2), sigmapts(2, 3), sigmapts(2, 4), sigmapts(2, 5), sigmapts(2, 6));
    fprintf(fid,'sigmapts = %15.11f %15.11f %15.11f %15.11f %15.11f %15.11f\n', sigmapts(3, 1), sigmapts(3, 2), sigmapts(3, 3), sigmapts(3, 4), sigmapts(3, 5), sigmapts(3, 6));
    fprintf(fid,'sigmapts = %15.11f %15.11f %15.11f %15.11f %15.11f %15.11f\n', sigmapts(4, 1), sigmapts(4, 2), sigmapts(4, 3), sigmapts(4, 4), sigmapts(4, 5), sigmapts(4, 6));
    fprintf(fid,'sigmapts = %15.11f %15.11f %15.11f %15.11f %15.11f %15.11f\n', sigmapts(5, 1), sigmapts(5, 2), sigmapts(5, 3), sigmapts(5, 4), sigmapts(5, 5), sigmapts(5, 6));
    fprintf(fid,'sigmapts = %15.11f %15.11f %15.11f %15.11f %15.11f %15.11f\n', sigmapts(6, 1), sigmapts(6, 2), sigmapts(6, 3), sigmapts(6, 4), sigmapts(6, 5), sigmapts(6, 6));

end

function testposcov2pts(fid)
    r1 = [5102.50895791863; 6123.01140072233; 6378.13692818659];
    cov2 = zeros(3,3);

    cov2(1, 1) = 12559.93762571587;
    cov2(1, 2) =  12101.56371305036;
    cov2(2, 1) = cov2(1, 2);
    cov2(1, 3) =  -440.3145384949657;
    cov2(3, 1) = cov2(1, 3);
    cov2(2, 2) = 12017.77368889201;
    cov2(2, 3) =  270.3798093532698;
    cov2(3, 2) = cov2(2, 3);
    cov2(3, 3) = 4818.009967057008;
 
    % form sigmapts pos/vel
    [sigmapts] = poscov2pts(r1', cov2);
    fprintf(fid,'sigmapts = %15.11f %15.11f %15.11f\n', sigmapts(1, 1), sigmapts(1, 2), sigmapts(1, 3));
    fprintf(fid,'sigmapts = %15.11f %15.11f %15.11f\n', sigmapts(2, 1), sigmapts(2, 2), sigmapts(2, 3));
    fprintf(fid,'sigmapts = %15.11f %15.11f %15.11f\n', sigmapts(3, 1), sigmapts(3, 2), sigmapts(3, 3));

    % reassemble covariance at each step and write out
    % [yu, covout] = remakecovpv(sigmapts);
    % [sigmapts] = poscov2pts(reci, cov);
    %
    % fprintf(fid,'sigmapts = %15.11f %15.11f %15.11f %15.11f %15.11f %15.11f\n', sigmapts(1, 1), sigmapts(1, 2), sigmapts(1, 3), sigmapts(1, 4), sigmapts(1, 5), sigmapts(1, 6));
    % fprintf(fid,'sigmapts = %15.11f %15.11f %15.11f %15.11f %15.11f %15.11f\n', sigmapts(2, 1), sigmapts(2, 2), sigmapts(2, 3), sigmapts(2, 4), sigmapts(2, 5), sigmapts(2, 6));
    % fprintf(fid,'sigmapts = %15.11f %15.11f %15.11f %15.11f %15.11f %15.11f\n', sigmapts(3, 1), sigmapts(3, 2), sigmapts(3, 3), sigmapts(3, 4), sigmapts(3, 5), sigmapts(3, 6));
end

function testremakecovp(fid)
    sigmapts = [5377.0260353 4827.9918806 5102.5089579 5102.5089579 5102.5089579 5102.5089579; ...
        6123.0114007 6123.0114007 6391.5382004 5854.4846011 6123.0114007 6123.0114007; ...
        6378.1369282 6378.1369282 6378.1369282 6378.1369282 6548.1606318 6208.1132245];

    % [yu, sigmapts] = remakecovpv(sigmapts);
    %
    % fprintf(fid,'cov = %15.11f %15.11f %15.11f\n', cov(1, 1), cov(1, 2), cov(1, 3));
    % fprintf(fid,'cov = %15.11f %15.11f %15.11f\n', cov(2, 1), cov(2, 2), cov(2, 3));
    % fprintf(fid,'cov = %15.11f %15.11f %15.11f\n', cov(3, 1), cov(3, 2), cov(3, 3));
end


function testremakecovpv(fid)
    sigmapts = [5377.0260353 4827.9918806 5102.5089579 5102.5089579 5102.5089579 5102.5089579; ...
        6123.0114007 6123.0114007 6391.5382004 5854.4846011 6123.0114007 6123.0114007; ...
        6378.1369282 6378.1369282 6378.1369282 6378.1369282 6548.1606318 6208.1132245];

    % [yu, sigmapts] = remakecovpv(sigmapts);
    %
    % fprintf(fid,'cov = %15.11f %15.11f %15.11f\n', cov(1, 1), cov(1, 2), cov(1, 3));
    % fprintf(fid,'cov = %15.11f %15.11f %15.11f\n', cov(2, 1), cov(2, 2), cov(2, 3));
    % fprintf(fid,'cov = %15.11f %15.11f %15.11f\n', cov(3, 1), cov(3, 2), cov(3, 3));
end


function testmatequal(fid)
    mat1 = [ 1, 3, 3  ; 1, 4, 3 ;  1, 3, 4  ];

    mat3 = mat1;

    fprintf(fid,'matequal = %15.11f %15.11f %15.11f\n', mat3(1, 1), mat3(1, 2), mat3(1, 3));
    fprintf(fid,'matequal = %15.11f %15.11f %15.11f\n', mat3(2, 1), mat3(2, 2), mat3(2, 3));
    fprintf(fid,'matequal = %15.11f %15.11f %15.11f\n', mat3(3, 1), mat3(3, 2), mat3(3, 3));
end

function testmatscale(fid)
    scale = 1.364;
    mat1 = [ 1, 3, 3  ; 1, 4, 3 ;  1, 3, 4  ];

    mat3 = mat1 * scale;

    fprintf(fid,'matscale = %15.11f %15.11f %15.11f\n', mat3(1, 1), mat3(1, 2), mat3(1, 3));
    fprintf(fid,'matscale = %15.11f %15.11f %15.11f\n', mat3(2, 1), mat3(2, 2), mat3(2, 3));
    fprintf(fid,'matscale = %15.11f %15.11f %15.11f\n', mat3(3, 1), mat3(3, 2), mat3(3, 3));
end


function testunit(fid)
    vec1 = [ 2.3; 4.7; -1.6 ];

    vec2 = unit(vec1);

    fprintf(fid,'norm = %15.11f %15.11f %15.11f\n', vec2(1), vec2(2), vec2(3));
end


function testmag(fid)
    x = [ 1.0; 2.0; 5.0 ];

    magx = mag(x);

    fprintf(fid,'mag = %15.11f\n', magx);
end


function testcross(fid)
    vec1 = [ 1.0; 2.0; 5.0 ];
    vec2 = [ 2.3; 4.7; -1.6 ];

    [outvec] = cross(vec1, vec2);

    fprintf(fid,'cross = %15.11f %15.11f %15.11f\n', outvec(1), outvec(2), outvec(3));
end


function testdot(fid)
    x = [ 1; 2; 5 ];
    y = [ 2.3; 4.7; -1.6 ];

    dotp = dot(x, y);

    fprintf(fid,'x %15.11f %15.11f %15.11f\n',  x(1), x(2), x(3));
    fprintf(fid,'y %15.11f %15.11f %15.11f\n', y(1), y(2), y(3));

    fprintf(fid,'dot = %15.11f\n', dotp);
end


function testangle(fid)
    vec1 = [ 1; 2; 5 ];
    vec2 = [ 2.3; 4.7; -1.6 ];

    ang = angl(vec1, vec2);

    fprintf(fid,'angle = %15.11f\n', ang);
end


function testasinh(fid)
    xval = 1.45;

    ansr = asinh(xval);

    fprintf(fid,'asinh = %15.11f\n', ansr);
end


function testcot(fid)
    xval = 0.47238734;

    ansr = cot(xval);

    fprintf(fid,'cot = %15.11f\n',ansr);
end


function testacosh(fid)
    xval = 1.43;

    ansr = acosh(xval);

    fprintf(fid,'acosh = %15.11f\n', ansr);
end


function testaddvec(fid)
    vec1 = [ 1; 2; 5 ];
    vec2 = [ 2.3; 4.7; -5.6 ];
    a1 = 1.0;
    a2 = 2.0;

    vec3 = (a1* vec1) + (a2 * vec2);

    fprintf(fid,'vec1 %15.11f %15.11f  %15.11f\n', vec1(1), vec1(2), vec1(3));
    fprintf(fid,'vec2 %15.11f %15.11f  %15.11f\n', vec2(1), vec2(2), vec2(3));

    fprintf(fid,'addvec = %15.11f %15.11f  %15.11f\n', vec3(1), vec3(2), vec3(3));
end


function testPercentile(fid)
    excelPercentile = 0.3;
    arrSize = 7;
    sequence(1) = 45.3;
    sequence(2) = 5.63;
    sequence(3) = 5.13;
    sequence(4) = 345.3;
    sequence(5) = 45.3;
    sequence(6) = 3445.3;
    sequence(7) = 0.03;

    ansr = prctile(sequence, excelPercentile, arrSize);

    fprintf(fid,'percentile = %15.11f\n', ansr);
end


function testrot1(fid)
    rad = 180.0 / pi;
    vec(1) = 23.4;
    vec(2) = 6723.4;
    vec(3) = -2.4;
    xval = 225.0 / rad;

    vec3 = rot1(vec, xval);

    fprintf(fid,'testrot1 = %15.11f  %15.11f  %15.11f\n', vec3(1), vec3(2), vec3(3));
end


function testrot2(fid)
    rad = 180.0 / pi;
    vec(1) = 23.4;
    vec(2) = 6723.4;
    vec(3) = -2.4;
    xval = 23.4 / rad;
    vec3 = rot2(vec, xval);

    fprintf(fid,'testrot2 = %15.11f  %15.11f  %15.11f\n', vec3(1), vec3(2), vec3(3));
end


function testrot3(fid)
    rad = 180.0 / pi;
    vec(1) = 23.4;
    vec(2) = 6723.4;
    vec(3) = -2.4;
    xval = 323.4 / rad;
    vec3 = rot3(vec, xval);

    fprintf(fid,'testrot3 = %15.11f  %15.11f  %15.11f\n', vec3(1), vec3(2), vec3(3));
end


function testfactorial(fid)
    n = 4;

    ansr = factorial(n);

    fprintf(fid,'factorial = %f\n', ansr);
end


function testcubicspl(fid)
    p1 = 1.0;
    p2 = 3.5;
    p3 = 5.6;
    p4 = 32.0;

    [acu0, acu1, acu2, acu3] = cubicspl(p1, p2, p3, p4);

    fprintf(fid,'cubicspl = %15.11f  %15.11f  %15.11f  %15.11f\n', acu0, acu1, acu2, acu3);
end


function testcubic(fid)
    a3 = 1.7;
    b2 = 3.5;
    c1 = 5.6;
    d0 = 32.0;
    opt = 'I';  % all roots, unique, real

    [r1r, r1i, r2r, r2i, r3r, r3i] = cubic(a3, b2, c1, d0, opt);

    fprintf(fid,'cubic = %15.11f  %15.11f  %15.11f  %15.11f  %15.11f  %15.11f\n', r1r, r1i, r2r, r2i, r3r, r3i);
end


function testcubicinterp(fid)
    p1a = 1.7;
    p1b = 3.5;
    p1c = 5.6;
    p1d = 11.7;
    p2a = 21.7;
    p2b = 35.5;
    p2c = 57.6;
    p2d = 181.7;
    valuein = 4.0;

    ansr = cubicinterp(p1a, p1b, p1c, p1d, p2a, p2b, p2c, p2d, valuein);

    fprintf(fid,'cubicint = %15.11f\n', ansr);
end

function testquadratic(fid)
    a = 1.7;
    b = 3.5;
    c = 5.6;
    opt = 'I';  % all roots, unique, real

    [r1r, r1i, r2r, r2i] = quadric(a, b, c, opt);

    fprintf(fid,'quad = %15.11f  %15.11f  %15.11f  %15.11f\n', r1r, r1i, r2r, r2i);
end

function testconvertMonth(fid)
    monstr = 'Jan';

    mon = getintmo(monstr); 
end

function testjday(fid)
    year = 2020;
    mon = 12;
    day = 15;
    hr = 16;
    minute = 58;
    second = 50.208;
    [jd, jdfrac] = jday(year, mon, day, hr, minute, second);

    fprintf(fid,'jd %15.11f  %15.11f\n', jd, jdfrac);
    % alt tests
    [year,mon,day,hr,minute,secs] = invjday ( 2450382.5, jdfrac );
    fprintf(fid,'year %5i   mon %4i day %3i %3i:%2i:%8.6f\n',year, mon, day, hr, minute, secs);

    [year,mon,day,hr,minute,secs] = invjday ( 2450382.5, jdfrac - 0.2 );
    fprintf(fid,'year %5i   mon %4i day %3i %3i:%2i:%8.6f\n',year, mon, day, hr, minute, secs);

    [year,mon,day,hr,minute,secs] = invjday ( 2450382.5 + 1, jdfrac + 1.5 );
    fprintf(fid,'year %5i   mon %4i day %3i %3i:%2i:%8.6f\n',year, mon, day, hr, minute, secs);

    [year,mon,day,hr,minute,secs] = invjday ( 2450382.5, -0.5 );
    fprintf(fid,'year %5i   mon %4i day %3i %3i:%2i:%8.6f\n',year, mon, day, hr, minute, secs);

    [year,mon,day,hr,minute,secs] = invjday ( 2450382.5, +0.5 );
    fprintf(fid,'year %5i   mon %4i day %3i %3i:%2i:%8.6f\n',year, mon, day, hr, minute, secs);

end


function testdays2mdhms(fid)
    year = 2020;
    days = 237.456982367;
    [mon, day, hr, minute, second] = days2mdh(year, days);

    fprintf(fid,'year %5i   days %14.9f\n',year, days);
    fprintf(fid,'year %5i   mon %4i day %3i %3i:%2i:%8.6f\n',year, mon, day, hr, minute, second);
end


function testinvjday(fid)
    jd = 2449877.0;
    jdFrac = 0.3458762;

    [year, mon, day, hr, minute, second] = invjday(jd, jdFrac);

    fprintf(fid,'year %5i   mon %4i day %3i %3i:%2i:%8.6f\n',year, mon, day, hr, minute, second);

    % some stressing cases
    for i = 1:11
        if (i == 1)
            jd = 2457884.5;
            jdF = 0.160911640046296;
        end
        if (i == 2)
            jd = 2457884.5;
            jdF = -0.160911640046296;
        end
        if (i == 3)
            jd = 2457884.5;
            jdF = 0.660911640046296;
        end
        if (i == 4)
            jd = 2457884.5;
            jdF = -0.660911640046296;
        end
        if (i == 5)
            jd = 2457884.5;
            jdF = 2.160911640046296;
        end
        if (i == 6)
            jd = 2457884.660911640046296;
            jdF = 0.0;
        end
        if (i == 7)
            jd = 2457884.0;
            jdF = 2.160911640046296;
        end
        if (i == 8)
            jd = 2457884.5;
            jdF = 0.0;
        end
        if (i == 9)
            jd = 2457884.5;
            jdF = 0.5;
        end
        if (i == 10)
            jd = 2457884.5;
            jdF = 1.0;
        end
        if (i == 11)
            jd = 2457884.3;
            jdF = 1.0;
        end

        jdb = floor(jd + jdF) + 0.5;
        mfme = (jd + jdF - jdb) * 1440.0;
        if (mfme < 0.0)
            mfme = 1440.0 + mfme;
        end
        [year,mon,day,hr,minute,second] = invjday ( jd, jdF );
        if (abs(jdF) >= 1.0)
            jd = jd + floor(jdF);
            jdF = jdF - floor(jdF);
        end
        dt = jd - floor(jd) - 0.5;
        if (abs(dt) > 0.00000001)
            jd = jd - dt;
            jdF = jdF + dt;
        end
        % this gets it even to the day
        if (jdF < 0.0)
            jd = jd - 1.0;
            jdF = 1.0 + jdF;
        end
        fprintf(fid,'%2i %4i %2i %2i %2i:%2i:%7.4f   %9.4f %9.4f %12.4f %8.5f %5.4f\n',...
            i, year, mon, day, hr, minute, second, mfme, hr*60.0+minute+second/60.0, jd, jdF, dt);
    end  % through stressing cases

    % pause;
    % 
    % [jdo, jdfraco] = jday(2017, 8, 23, 12, 15, 16);
    % 
    % [jdo, jdfraco] = jday(2017, 12, 31, 12, 15, 16);
    % fprintf(fid,'%11.6f  %11.6f\n',jdo,jdfraco);
    % [year,mon,day,hr,minute,second] = invjday ( jdo + jdfraco, 0.0 );
    % fprintf(fid,'%4i  %3i  %3i  %2i:%2i:%6.4f\n', year,mon,day,hr,minute,second);
    % [year,mon,day,hr,minute,second] = invjday ( jdo, jdfraco );
    % fprintf(fid,'%4i  %3i  %3i  %2i:%2i:%6.4f\n', year,mon,day,hr,minute,second);
    % 
    % for i = -50:50
    %     jd = jdo + i/10.0;
    %     jdfrac = jdfraco;
    %     [year,mon,day,hr,minute,second] = invjday ( jd + jdfrac, 0.0 );
    %     fprintf(fid,'%4i  %3i  %3i  %2i:%2i:%6.4f\n', year,mon,day,hr,minute,second);
    %     [year,mon,day,hr,minute,second] = invjday ( jd, jdfrac );
    %     dt = jd - floor(jd);
    %     fprintf(fid,'%4i  %3i  %3i  %2i:%2i:%6.4f  \n\n', year,mon,day,hr,minute,second);
    % end
    % 
    % fprintf(fid,'end first half\n');
    % 
    % for i = -50:50
    %     jd = jdo;
    %     jdfrac = jdfraco + i/10.0;
    %     [year,mon,day,hr,minute,second] = invjday ( jd + jdfrac, 0.0 );
    %     fprintf(fid,'%4i  %3i  %3i  %2i:%2i:%6.4f\n', year,mon,day,hr,minute,second);
    %     [year,mon,day,hr,minute,second] = invjday ( jd, jdfrac );
    %     fprintf(fid,'%4i  %3i  %3i  %2i:%2i:%6.4f  \n\n', year,mon,day,hr,minute,second);
    % end
end   % testinvjday


% tests eop, spw, and fk5 iau80
function testiau80in(fid)
    year = 2017;
    mon = 4;
    day = 6;
    hr = 0;
    minute = 0;
    second = 0.0;

    fileLoc = 'D:\Codes\LIBRARY\DataLib\';
    [iau80arr] = iau80in(fileLoc);

    eopFileName = 'D:\Codes\LIBRARY\DataLib\EOP-All-v1.1_2018-01-04.txt';
    [eoparr] = readeop(eopFileName);

    fprintf(fid,'EOP tests  mfme    dut1  dat    lod           xp                      yp               ddpsi                   ddeps               ddx                 ddy\n');
    for i = 0: 90

        [jd, jdFrac] = jday(year, mon, day, hr + i, minute, second);
        [dut1, dat, lod, xp, yp, ddpsi,ddeps, ddx, ddy] = findeopparam(jd, jdFrac, 's', eoparr);
        [y, m, d, h, mm, ss] = invjday(jd, jdFrac);
        fprintf(fid,' %i  %i  %i  %i  %15.11f  %3i  %15.11f  %15.11g  %15.11g  %15.11g  %15.11g  %15.11g  %15.11g\n', ...
            y, m, d, h * 60 + mm, dut1, dat, lod, xp, yp, ddpsi, ddeps, ddx, ddy);
    end

    spwFileName = 'D:\Codes\LIBRARY\DataLib\SpaceWeather-All-v1.2_2018-01-04.txt';
    [spwarr, updDate] = readspw(spwFileName);
    fprintf(fid,'SPW tests  mfme f107 f107bar ap apavg  kp sumkp aparr[]\n');
    for i = 0: 45

        [jd, jdF] = jday(year, mon, day, hr + i, minute, second);
        % adj obs, last ctr, act con
        [f107, f107bar, ap, avgap, aparr, kp, sumkp, kparr] = findspwparam( jd, jdF, 's', 'a', 'l', 'a', spwarr);
        [y, m, d, h, mm, ss] = invjday(jd, jdF);
        fprintf(fid,'sala %i  %i  %i  %i  %15.11f  %15.11f  %15.11f  %15.11f  %15.11f  %15.11f  %15.11f %15.11f  %15.11f\n', ...
            y, m, d, h * 60 + mm, f107, f107bar, ap, avgap, kp, sumkp, aparr(1), aparr(2), aparr(3));
       
       [f107, f107bar, ap, avgap, aparr, kp, sumkp, kparr] = findspwparam( jd, jdF, 's', 'o', 'c', 'a', spwarr);
        [y, m, d, h, mm, ss] = invjday(jd, jdF);
        fprintf(fid,'soca %i  %i  %i  %i  %15.11f  %15.11f  %15.11f  %15.11f  %15.11f  %15.11f  %15.11f %15.11f  %15.11f\n', ...
            y, m, d, h * 60 + mm, f107, f107bar, ap, avgap, kp, sumkp, aparr(1), aparr(2), aparr(3));
    end
end  % testiau80in

function testfundarg(fid)
    opt = '80';
    ttt = 0.042623631889;

    [fArgs] = fundarg(ttt, opt);
    fprintf(fid,'fundarg = %15.11f %15.11f  %15.11f %15.11f  %15.11f %15.11f  %15.11f %15.11f  %15.11f %15.11f  %15.11f %15.11f  %15.11f %15.11f\n', ...
        fArgs(1), fArgs(2), fArgs(3), fArgs(4), fArgs(5), fArgs(6), fArgs(7), fArgs(8),...
        fArgs(9), fArgs(10), fArgs(11), fArgs(12), fArgs(13), fArgs(14));

    % do in deg
    for i = 1: 14
        fArgs(i) = fArgs(i) * 180.0 / pi;
    end

    fprintf(fid,'fundarg = %15.11f %15.11f  %15.11f %15.11f  %15.11f %15.11f  %15.11f %15.11f  %15.11f %15.11f  %15.11f %15.11f  %15.11f %15.11f\n', ...
        fArgs(1), fArgs(2), fArgs(3), fArgs(4), fArgs(5), fArgs(6), fArgs(7), fArgs(8),...
        fArgs(9), fArgs(10), fArgs(11), fArgs(12), fArgs(13), fArgs(14));


    [fArgs] = fundarg(ttt, '06');
    fprintf(fid,'fundarg = %15.11f %15.11f  %15.11f %15.11f  %15.11f %15.11f  %15.11f %15.11f  %15.11f %15.11f  %15.11f %15.11f  %15.11f %15.11f\n', ...
        fArgs(1), fArgs(2), fArgs(3), fArgs(4), fArgs(5), fArgs(6), fArgs(7), fArgs(8),...
        fArgs(9), fArgs(10), fArgs(11), fArgs(12), fArgs(13), fArgs(14));

    % do in deg
    for i = 1: 14
        fArgs(i) = fArgs(i) * 180.0 / pi;
    end

    fprintf(fid,'fundarg = %15.11f %15.11f  %15.11f %15.11f  %15.11f %15.11f  %15.11f %15.11f  %15.11f %15.11f  %15.11f %15.11f  %15.11f %15.11f\n', ...
        fArgs(1), fArgs(2), fArgs(3), fArgs(4), fArgs(5), fArgs(6), fArgs(7), fArgs(8),...
        fArgs(9), fArgs(10), fArgs(11), fArgs(12), fArgs(13), fArgs(14));
end  % testfundarg


function testprecess(fid)
    opt = '80';
    ttt = 0.042623631889;
    % ttt = 0.04262362174880504;
    [prec, psia, wa, epsa, chia] = precess(ttt, opt);

    fprintf(fid,'prec = %15.11f  %15.11f  %15.11f\n', prec(1, 1), prec(1, 2), prec(1, 3));
    fprintf(fid,'prec = %15.11f  %15.11f  %15.11f\n', prec(2, 1), prec(2, 2), prec(2, 3));
    fprintf(fid,'prec = %15.11f  %15.11f  %15.11f\n', prec(3, 1), prec(3, 2), prec(3, 3));

    [psia, wa, epsa, chia] = precess(ttt, '06');

    fprintf(fid,'prec06 = %15.11f  %15.11f  %15.11f\n', prec(1, 1), prec(1, 2), prec(1, 3));
    fprintf(fid,'prec06 = %15.11f  %15.11f  %15.11f\n', prec(2, 1), prec(2, 2), prec(2, 3));
    fprintf(fid,'prec06 = %15.11f  %15.11f  %15.11f\n', prec(3, 1), prec(3, 2), prec(3, 3));
end


function testnutation(fid)
    opt = '80';
    fileLoc = 'D:\Codes\LIBRARY\DataLib\';
    [iau80arr] = iau80in(fileLoc);
    fileLoc = 'D:\Codes\LIBRARY\DataLib\';
    [iau06arr] = iau06in(fileLoc);

    ttt = 0.042623631889;
    ddpsi = -0.052195;
    ddeps = -0.003875;

    [fArgs] = fundarg(ttt, opt);

    [deltapsi, trueeps, meaneps, nut] = nutation  (ttt, ddpsi, ddeps, iau80arr, fArgs);

    fprintf(fid,'nut = %15.11f  %15.11f  %15.11f\n', nut(1, 1), nut(1, 2), nut(1, 3));
    fprintf(fid,'nut = %15.11f  %15.11f  %15.11f\n', nut(2, 1), nut(2, 2), nut(2, 3));
    fprintf(fid,'nut = %15.11f  %15.11f  %15.11f\n', nut(3, 1), nut(3, 2), nut(3, 3));

    [fArgs] = fundarg(ttt, '06');
    % nut00 = precnutbias00a(ttt, ddpsi, ddeps, iau06arr, '06', fArgs);
    % fprintf(fid,'nut06 c= %15.11f  %15.11f  %15.11f\n', nut(1, 1), nut(1, 2), nut(1, 3));
    % fprintf(fid,'nut06 c= %15.11f  %15.11f  %15.11f\n', nut(2, 1), nut(2, 2), nut(2, 3));
    % fprintf(fid,'nut06 c= %15.11f  %15.11f  %15.11f\n', nut(3, 1), nut(3, 2), nut(3, 3));

    %fundarg(ttt, AstroLib.EOpt.e00a, out fArgs);
    %nut00 = nutation00a(ttt, ddpsi, ddeps, iau06arr, AstroLib.EOpt.e00a);
    %fprintf(fid,'nut06 a= ' + nut(1, 1), nut(1, 2), nut(1, 3));
    %fprintf(fid,'nut06 a= ' + nut(2, 1), nut(2, 2), nut(2, 3));
    %fprintf(fid,'nut06 a= ' + nut(3, 1), nut(3, 2), nut(3, 3));
end


function testnutationqmod(fid)
    opt = '80';
    string fileLoc;
    fileLoc = 'D:\Codes\LIBRARY\DataLib\';
    [iau80arr] = iau80in(fileLoc);
    ttt = 0.042623631889;
    % ttt = 0.04262362174880504;

    [fArgs] = fundarg(ttt, opt);

  %  [opt, fArgs, deltapsi, deltaeps, meaneps, nutq] = nutationqmod(ttt, iau80arr);
end


function testsidereal(fid)
    eqeterms = 2;
    opt = '80';
    jdut1 = 2453101.82740678310;
    ttt = 0.042623631889;
    %ttt = 0.04262362174880504;
    lod = 0.001556;

    fileLoc = 'D:\Codes\LIBRARY\DataLib\';
    [iau80arr] = iau80in(fileLoc);

    ddpsi = -0.052195;
    ddeps = -0.003875;

    [fArgs] = fundarg(ttt, opt);
    [deltapsi, trueeps, meaneps, nut] = nutation(ttt, ddpsi, ddeps, iau80arr, fArgs);
    [st, stdot]  = sidereal (jdut1, deltapsi, meaneps, fArgs(5), lod, eqeterms, '80' );
    fprintf(fid,'st = %15.11f  %15.11f  %15.11f\n', st(1, 1), st(1, 2), st(1, 3));
    fprintf(fid,'st = %15.11f  %15.11f  %15.11f\n', st(2, 1), st(2, 2), st(2, 3));
    fprintf(fid,'st = %15.11f  %15.11f  %15.11f\n', st(3, 1), st(3, 2), st(3, 3));


    [st, stdot]  = sidereal (jdut1, deltapsi, meaneps, fArgs(5), lod, eqeterms, '06' );
    fprintf(fid,'st00 = %15.11f  %15.11f  %15.11f\n', st(1, 1), st(1, 2), st(1, 3));
    fprintf(fid,'st00 = %15.11f  %15.11f  %15.11f\n', st(2, 1), st(2, 2), st(2, 3));
    fprintf(fid,'st00 = %15.11f  %15.11f  %15.11f\n', st(3, 1), st(3, 2), st(3, 3));
end


function testpolarm(fid)
    opt = '80';

    ttt = 0.042623631889;
    %ttt = 0.04262363188899416;
    xp = 0.0;
    yp = 0.0;

    pm = polarm(xp, yp, ttt, opt);

    fprintf(fid,'pm = %15.11f  %15.11f  %15.11f\n', pm(1, 1), pm(1, 2), pm(1, 3));
    fprintf(fid,'pm = %15.11f  %15.11f  %15.11f\n', pm(2, 1), pm(2, 2), pm(2, 3));
    fprintf(fid,'pm = %15.11f  %15.11f  %15.11f\n', pm(3, 1), pm(3, 2), pm(3, 3));

    pm = polarm(xp, yp, ttt, '06');

    fprintf(fid,'pm06 = %15.11f  %15.11f  %15.11f\n', pm(1, 1), pm(1, 2), pm(1, 3));
    fprintf(fid,'pm06 = %15.11f  %15.11f  %15.11f\n', pm(2, 1), pm(2, 2), pm(2, 3));
    fprintf(fid,'pm06 = %15.11f  %15.11f  %15.11f\n', pm(3, 1), pm(3, 2), pm(3, 3));
end


function testgstime(fid)
    jdut1 = 2453101.82740678310;

    gst = gstime(jdut1);

    fprintf(fid,'gst = %15.11f  %15.11f\n', gst, (gst * 180.0 / pi));
end


function testlstime(fid)
    rad = 180.0 / pi;
 
    lon = -104.0 / rad;
    jdut1 = 2453101.82740678310;

    [lst, gst] = lstime(lon, jdut1);

    fprintf(fid,'lst = %15.11f  %15.11f\n', lst, (lst * 180.0 / pi));
end


function testhms_sec(fid)
    hr = 12;
    minute = 34;
    second = 56.233;
    fprintf(fid,' hr minute second %4i  %4i  %8.6f\n',hr, minute, second);

    [utsec ] = hms2sec( hr,minute,second );
    fprintf(fid,'utsec = %15.11f\n', utsec);

    [hr,minute,second ] = sec2hms( utsec);
    fprintf(fid,' hr minute second %4i  %4i  %8.6f\n',hr, minute, second);
end


function testhms_ut(fid)
    hr = 13;
    minute = 22;
    second = 45.98;

    fprintf(fid,' hr minute second %4i  %4i  %8.6f\n',hr, minute, second);

    [ut] = hms2ut(hr, minute, second);
    fprintf(fid,'ut = %15.11f\n', ut);

    [hr,minute,second] = ut2hms( ut );
    fprintf(fid,' hr minute second %4i  %4i  %8.6f\n',hr, minute, second);

end


function testhms_rad(fid)
    hr = 15;
    minute = 15;
    second = 53.63;

    fprintf(fid,' hr minute second %4d  %4d  %8.6f\n',hr, minute, second);
    
    [hms] = hms2rad( hr,minute,second );
    fprintf(fid,'hms %15.11f \n',hms);

    [hr,minute,second] = rad2hms( hms );
    fprintf(fid,' hr minute second %4d  %4d  %8.6f\n',hr, minute, second);
end


function testdms_rad(fid)
    deg = -35;
    minute = -15;
    second = -53.63;

    fprintf(fid,' deg minute second %4d  %4d  %8.6f\n', deg, minute, second);

    [dms] = dms2rad( deg, minute, second );
    fprintf(fid,'dms = %15.11f\n', dms);

    [deg, minute, second] = rad2dms( dms );
    fprintf(fid,' deg minute second %4d  %4d  %8.6f\n', deg, minute, second);
end


function testeci_ecef(fid)
    conv = pi / (180.0 * 3600.0);  % arcsec to rad

    recef = [ -1033.4793830; 7901.2952754; 6380.3565958 ];
    vecef = [ -3.225636520; -2.872451450; 5.531924446 ];
    aecef = [  0.0; 0.0; 0.0 ];
    year = 2004;
    mon = 4;
    day = 6;
    hr = 7;
    minute = 51;
    second = 28.386009;

    % sofa example -------------
    %year = 2007;
    %mon = 4;
    %day = 5;
    %hr = 12;
    %minute = 0;
    %second = 0.0;

    [jd, jdFrac] = jday(year, mon, day, hr, minute, second);

    dut1 = -0.4399619;      % second
    dat = 32;               % second
    xp = -0.140682 * conv;  % ' to rad
    yp = 0.333309 * conv;
    lod = 0.0015563;
    ddpsi = -0.052195 * conv;  % ' to rad
    ddeps = -0.003875 * conv;
    ddx = -0.000205 * conv;    % ' to rad
    ddy = -0.000136 * conv;
    eqeterms = 2;

    % sofa example -------------
    %dut1 = -0.072073685;      % second
    %dat = 33;                % second
    %xp = 0.0349282 * conv;  % ' to rad
    %yp = 0.4833163 * conv;
    %lod = 0.0;
    %ddpsi = -0.0550655 * conv;  % ' to rad
    %ddeps = -0.006358 * conv;
    %ddx = 0.000175 * conv;    % ' to rad
    %ddy = -0.0002259 * conv;

    jdtt = jd;
    jdftt = jdFrac + (dat + 32.184) / 86400.0;
    ttt = (jdtt + jdftt - 2451545.0) / 36525.0;
    fprintf(fid,'ttt wo base (use this) %15.11f\n', ttt);
    jdut1 = jd + jdFrac + dut1 / 86400.0;

    fprintf(fid,'ITRF          IAU-76/FK5  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', recef(1), recef(2), recef(3), vecef(1), vecef(2), vecef(3));

    fileLoc = 'D:\Codes\LIBRARY\DataLib\';
    [iau80arr] = iau80in(fileLoc);

    fileLoc = 'D:\Codes\LIBRARY\DataLib\';
    [iau06arr] = iau06in(fileLoc);

    % test creating xys file
    %fileLoc = 'D:\Codes\LIBRARY\DataLib\';
    % done, works in c#. :-)
    %createXYS(fileLoc, iau06arr, fArgs);

    % now read it in
    fileLoc = 'D:\Codes\LIBRARY\DataLib\';
    [xys06table] = readxys(fileLoc);

    % now test it for interpolation
    %jdtt = jd + jdFrac + (dat + 32.184) / 86400.0;
    [fArgs06] = fundarg(ttt, '06');
    [x ,y, s] = iau06xysS (iau06arr, fArgs06, ttt );
    fprintf(fid,'iau06xysS wo   x   %11.9f  y  %11.9g  s %11.9g\n',x, y, s);
    fprintf(fid,'iau06xysS wo   x   %11.9f  y  %11.9g  s %11.9g arsec\n',x/conv, y/conv, s/conv);
    x = x + ddx;
    y = y + ddy;
    fprintf(fid,'iau06xysS      x   %11.9f  y  %11.9g  s %11.9g\n',x, y, s);
    fprintf(fid,'iau06xysS      x   %11.9f  y  %11.9g  s %11.9g arcsec\n',x/conv, y/conv, s/conv);
    [x, y, s] = findxysparam(jdtt + jdftt, 0.0, 'n', xys06table);
    fprintf(fid,'findxysparam n x   %11.9f  y  %11.9g  s  %11.9g\n', x, y, s);
    [x, y, s] = findxysparam(jdtt + jdftt, 0.0, 'l', xys06table);
    fprintf(fid,'findxysparam l x   %11.9f  y  %11.9g  s  %11.9g\n', x, y, s);
    [x, y, s] = findxysparam(jdtt + jdftt, 0.0, 's', xys06table);
    fprintf(fid,'findxysparam s x   %11.9f  y  %11.9g  s  %11.9g\n', x, y, s);
    [x ,y, s] = iau06xysS (iau06arr, fArgs06, ttt );

    [recii, vecii, aecii] = ecef2eci(recef, vecef, aecef, iau80arr, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps );
    fprintf(fid,'GCRF          IAU-76/FK5  %15.11f  %15.11f  %15.11f  %15.11f  %15.11f  %15.11f\n', recii(1), recii(2), recii(3), vecii(1), vecii(2), vecii(3));

    [reci veci, aeci] = ecef2eci06(recef, vecef, aecef,   iau06arr, xys06table, ttt, jdut1, lod, xp, yp, ddx, ddy, '06x');
    fprintf(fid,'GCRF          IAU-2006 CIO  %15.11f  %15.11f  %15.11f  %15.11f  %15.11f  %15.11f\n', reci(1), reci(2), reci(3), veci(1), veci(2), veci(3));

    % try backwards
    [recef,vecef,aecef] = eci2ecef  ( recii,vecii,aecii,iau80arr, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps );
    fprintf(fid,'ITRF rev       IAU-76/FK5  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', recef(1), recef(2), recef(3), vecef(1), vecef(2), vecef(3));

    % these are not correct
    [recef, vecef, aecef] = eci2ecef06(reci, veci, aeci, iau06arr, xys06table, ttt,jdut1,lod,xp,yp, ddx, ddy,'06x' );
    fprintf(fid,'ITRF rev       IAU-2006 CIO  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', recef(1), recef(2), recef(3), vecef(1), vecef(2), vecef(3));


    % writeout data for table interpolation
    year = 1980;
    mon = 1;
    day = 1;
    hr = 0;
    minute = 0;
    second = 0.0;
    interp = 'x';  % full series

    fileLoc = 'D:\Codes\LIBRARY\DataLib\';
    [iau06arr] = iau06in(fileLoc);

    % read interpolated one
    %EOPSPWLibr.initEOPArrayP(ref EOPSPWLibr.eopdataP);

    % read existing data - this does not find x, y, s!
    eopFileName = 'D:\Codes\LIBRARY\DataLib\EOP-All-v1.1_2018-01-04.txt';
    [eoparr] = readeop(eopFileName);

    % now find table of CIO values
    fileLoc = 'D:\Codes\LIBRARY\DataLib\';
    [xys06table] = readxys(fileLoc);

    % rad to '
     convrt = (180.0 * 3600.0) / pi;
    fprintf(fid,'CIO tests      x          y          s          ddpsi        ddeps     ddx       ddy      mjd\n');
    for i = 0:13 % 14500

        [jd, jdFrac] = jday(year, mon, day + i, hr, minute, second);
        [dut1, dat, lod, xp, yp, ddpsi, ddeps, ddx, ddy] = findeopparam(jd, jdFrac, 's', eoparr);
        jdtt = jd;
        jdftt = jdFrac + (dat + 32.184) / 86400.0;
        ttt = (jdtt + jdftt - 2451545.0) / 36525.0;

        [fArgs] = fundarg(ttt, '80');

        % ddpsi = 0.0;
        % ddeps = 0.0;
        deltapsi= 0.0;
        deltaeps= 0.0;
         for ii= 106:-1: 1
            tempval = iau80arr.iar80(ii, 1) * fArgs(1) + iau80arr.iar80(ii, 2) * fArgs(2) + iau80arr.iar80(ii, 3) * fArgs(3) + ...
                      iau80arr.iar80(ii, 4) * fArgs(4) + iau80arr.iar80(ii, 5) * fArgs(5);
            deltapsi = deltapsi + (iau80arr.rar80(ii, 1) + iau80arr.rar80(ii, 2) * ttt) * sin(tempval);
            deltaeps = deltaeps + (iau80arr.rar80(ii, 3) + iau80arr.rar80(ii, 4) * ttt) * cos(tempval);
        end

        % --------------- find nutation parameters --------------------
        deltapsi = (deltapsi + ddpsi); % (2.0 * pi);
        deltaeps = (deltaeps + ddeps); % (2.0 * pi);
        deltapsi = deltapsi * convrt;
        deltaeps = deltaeps * convrt;

        % CIO parameters
        [fArgs06] = fundarg(ttt, '06');
        % ddx = 0.0;
        % ddy = 0.0;
        [x, y, s, pn] = iau06xys (iau06arr, fArgs06, xys06table, ttt, ddx, ddy, '06x');
        x = x * convrt;
        y = y * convrt;
        s = s * convrt;

        fprintf(fid,'%15.11f %15.11f %15.11f %15.11f %15.11f %15.11f %15.11f %11.0f\n', x, y, s, deltapsi, deltaeps, ddx, ddy, (jd + jdFrac - 2400000.5));
    end
end


function testtod2ecef(fid)
    conv = pi / (180.0 * 3600.0);

    recef = [ -1033.4793830; 7901.2952754; 6380.3565958 ];
    vecef = [ -3.225636520; -2.872451450; 5.531924446 ];
    aecef = [0.00001; 0.00001; 0.00001];
    year = 2004;
    mon = 4;
    day = 6;
    hr = 7;
    minute = 51;
    second = 28.386009;
    [jd, jdFrac] = jday(year, mon, day, hr, minute, second);

    dut1 = -0.4399619;      % second
    dat = 32;               % second
    xp = -0.140682 * conv;  % ' to rad
    yp = 0.333309 * conv;
    lod = 0.0015563;
    ddpsi = -0.052195 * conv;  % ' to rad
    ddeps = -0.003875 * conv;
    ddx = -0.000205 * conv;    % ' to rad
    ddy = -0.000136 * conv;
    eqeterms = 2;

    reci = [ 0.0; 0.0; 0.0 ];
    veci = [ 0.0; 0.0; 0.0 ];

    fileLoc = 'D:\Codes\LIBRARY\DataLib\';
    [iau80arr] = iau80in(fileLoc);
    fileLoc = 'D:\Codes\LIBRARY\DataLib\';
    [iau06arr] = iau06in(fileLoc);
    jdtt = jd;
    jdftt = jdFrac + (dat + 32.184) / 86400.0;
    ttt = (jdtt + jdftt - 2451545.0) / 36525.0;

    % now read it in
    fileLoc = 'D:\Codes\LIBRARY\DataLib\';
    [xys06table] = readxys(fileLoc);

    jdut1 = jd + jdFrac + dut1 / 86400.0;

    fprintf(fid,'ITRF start    IAU-76/FK5  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', recef(1), recef(2), recef(3), vecef(1), vecef(2), vecef(3));

    % PEF
    [rpef, vpef, apef] = ecef2pef ( recef, vecef, aecef, '80', ttt, xp, yp );
    fprintf(fid,'PEF           IAU-76/FK5  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', rpef(1), rpef(2), rpef(3), vpef(1), vpef(2), vpef(3));
    [recii, vecii, aecii] = pef2ecef(rpef, vpef, apef, '80', ttt, xp, yp);
    fprintf(fid,'ITRF  rev     IAU-76/FK5  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', recii(1), recii(2), recii(3), vecii(1), vecii(2), vecii(3));

    [rtirs, vtirs, atirs] = ecef2tirs(recef, vecef, aecef, ttt, xp, yp );
    fprintf(fid,'TIRS          IAU-2006 CIO  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', rtirs(1), rtirs(2), rtirs(3), vtirs(1), vtirs(2), vtirs(3));
    [recefi, vecefi, aecefi] = tirs2ecef(rtirs, vtirs, atirs, ttt, xp, yp );
    fprintf(fid,'ECI rev       IAU-2006 CIO  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', recefi(1), recefi(2), recefi(3), vecefi(1), vecefi(2), vecefi(3));

    % TOD
    [rtod, vtod, atod] = ecef2tod(recef, vecef, aecef, iau80arr, ttt, jdut1, lod, xp, yp, eqeterms, 0.0, 0.0);
    fprintf(fid,'TOD wo corr   IAU-76/FK5  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', rtod(1), rtod(2), rtod(3), vtod(1), vtod(2), vtod(3));
    [rtod, vtod, atod] = ecef2tod(recef, vecef, aecef, iau80arr, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps );
    fprintf(fid,'TOD w corr    IAU-76/FK5  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', rtod(1), rtod(2), rtod(3), vtod(1), vtod(2), vtod(3));
    [recefi,vecefi,aecefi] = tod2ecef ( rtod, vtod, atod, iau80arr, jdut1, ttt, lod, xp, yp, eqeterms, ddpsi, ddeps );
    fprintf(fid,'ITRFi         IAU-76/FK5  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', recefi(1), recefi(2), recefi(3), vecefi(1), vecefi(2), vecefi(3));

    [rcirs, vcirs, acirs] = ecef2cirs(recef, vecef, aecef, iau06arr, ttt, jdut1, lod, xp, yp, ddx, ddy, '06x' );
    fprintf(fid,'CIRS          IAU-2006 CIO  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', rcirs(1), rcirs(2), rcirs(3), vcirs(1), vcirs(2), vcirs(3));
    [recefi,vecefi,aecefi] = cirs2ecef(rcirs, vcirs, acirs, iau06arr, ttt, jdut1, lod, xp, yp, ddx, ddy, '06x' );
    fprintf(fid,'ITRF rev      IAU-2006 CIO  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', recefi(1), recefi(2), recefi(3), vecefi(1), vecefi(2), vecefi(3));

    % MOD
    [rmod, vmod, amod] = ecef2mod( recef, vecef, aecef, iau80arr, ttt, jdut1, lod, xp, yp, eqeterms, 0.0, 0.0);
    fprintf(fid,'MOD wo corr   IAU-76/FK5  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', rmod(1), rmod(2), rmod(3), vmod(1), vmod(2), vmod(3));
    [rmod, vmod, amod] = ecef2mod( recef, vecef, aecef, iau80arr, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps );
    fprintf(fid,'MOD  w corr   IAU-76/FK5  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', rmod(1), rmod(2), rmod(3), vmod(1), vmod(2), vmod(3));
    [recefi,vecefi,aecefi] = mod2ecef( rmod, vmod, amod, iau80arr, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps );
    fprintf(fid,'ITRF  rev     IAU-76/FK5  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', recefi(1), recefi(2), recefi(3), vecefi(1), vecefi(2), vecefi(3));


    % J2000
    [recii, vecii, aecii] = ecef2eci(recef, vecef, aecef, iau80arr, ttt, jdut1, lod, xp, yp, eqeterms, 0.0, 0.0 );
    fprintf(fid,'J2000 wo corr IAU-76/FK5  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', recii(1), recii(2), recii(3), vecii(1), vecii(2), vecii(3));

    % GCRF
    [reci, veci, aeci] = ecef2eci(recef, vecef, aecef, iau80arr, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps );
    fprintf(fid,'GCRF w corr   IAU-76/FK5  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', reci(1), reci(2), reci(3), veci(1), veci(2), veci(3));
    [recefi, vecefi, aecefi] = eci2ecef(reci, veci, aeci, iau80arr, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps );
    fprintf(fid,'ITRF rev      IAU-76/FK5  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', recefi(1), recefi(2), recefi(3), vecefi(1), vecefi(2), vecefi(3));

    [recii, vecii, aecii] = ecef2eci06(recef, vecef, aecef,   iau06arr, xys06table, ttt, jdut1, lod, xp, yp, ddx, ddy, '06x');
    fprintf(fid,'GCRF          IAU-2006 CIO  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', recii(1), recii(2), recii(3), vecii(1), vecii(2), vecii(3));
    [recefi,vecefi,aecefi] = eci2ecef06  ( recii,vecii,aecii, iau06arr, xys06table, ttt, jdut1, lod, xp, yp, ddx, ddy, '06x' );    
    fprintf(fid,'ITRF rev      IAU-2006 CIO  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', recefi(1), recefi(2), recefi(3), vecefi(1), vecefi(2), vecefi(3));

    % sofa
    fprintf(fid,'SOFA ECI CIO  5102.508959486507   6123.011392959787   6378.136934384333\n');
    fprintf(fid,'SOFA ECI 06a  5102.508965811828   6123.011397147659   6378.136925303720\n');
    fprintf(fid,'SOFA ECI 00a  5102.508965732931   6123.011397847143   6378.136924695331\n');


    % now reverses from eci
    fprintf(fid,'GCRF wco STARTIAU-76/FK5 %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', reci(1), reci(2), reci(3), veci(1), veci(2), veci(3));

    % PEF
    [rpef, vpef, apef] = eci2pef(reci, veci, aeci, iau80arr, ttt, jdut1, lod, eqeterms, ddpsi, ddeps);
    fprintf(fid,'PEF           IAU-76/FK5  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', rpef(1), rpef(2), rpef(3), vpef(1), vpef(2), vpef(3));
    [recii, vecii, aecii] = pef2eci(rpef, vpef, apef, iau80arr, ttt, jdut1, lod, eqeterms, ddpsi, ddeps);
    fprintf(fid,'ECI rev       IAU-76/FK5  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', recii(1), recii(2), recii(3), vecii(1), vecii(2), vecii(3));

    [rtirs, vtirs, atirs] = eci2tirs(recii, vecii, aecii, iau06arr, xys06table, ttt, jdut1, lod, ddx, ddy, '06x' );
    fprintf(fid,'TIRS          IAU-2006 CIO  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', rtirs(1), rtirs(2), rtirs(3), vtirs(1), vtirs(2), vtirs(3));
    [recii, vecii, aecii] = tirs2eci(rtirs, vtirs, atirs, iau06arr, xys06table, ttt, jdut1, lod, ddx, ddy, '06x' );
    fprintf(fid,'ECI rev       IAU-2006 CIO  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', recii(1), recii(2), recii(3), vecii(1), vecii(2), vecii(3));

    % TOD
    [rtod, vtod, atod] = eci2tod(reci, veci, aeci, iau80arr, ttt, ddpsi, ddeps);
    fprintf(fid,'TOD           IAU-76/FK5  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', rtod(1), rtod(2), rtod(3), vtod(1), vtod(2), vtod(3));
    [recii, vecii, aecii] = tod2eci(rtod, vtod, atod, iau80arr, ttt, ddpsi, ddeps);
    fprintf(fid,'ECI rev       IAU-76/FK5  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', recii(1), recii(2), recii(3), vecii(1), vecii(2), vecii(3));

    [rcirs, vcirs, acirs] = eci2cirs( recii, vecii, aecii, iau06arr, xys06table, ttt, ddx, ddy, '06x' );
    fprintf(fid,'CIRS          IAU-2006 CIO  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', rcirs(1), rcirs(2), rcirs(3), vcirs(1), vcirs(2), vcirs(3));
    [recii, vecii, aecii] = cirs2eci(rcirs, vcirs, acirs, iau06arr,  xys06table, ttt, ddx, ddy, '06x');
    fprintf(fid,'ECI rev       IAU-2006 CIO  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', recii(1), recii(2), recii(3), vecii(1), vecii(2), vecii(3));

    % MOD
    [rmod, vmod, amod] = eci2mod(reci, veci, aeci, ttt);
    fprintf(fid,'MOD           IAU-76/FK5  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', rmod(1), rmod(2), rmod(3), vmod(1), vmod(2), vmod(3));
    [recii, vecii, aecii] = mod2eci(rmod, vmod, amod, ttt);
    fprintf(fid,'ECI rev       IAU-76/FK5  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', recii(1), recii(2), recii(3), vecii(1), vecii(2), vecii(3));
end


function testteme_ecef(fid)
    conv = pi / (180.0 * 3600.0);

    recef = [ -1033.4793830; 7901.2952754; 6380.3565958 ];
    vecef = [ -3.225636520; -2.872451450; 5.531924446 ];
    year = 2004;
    mon = 4;
    day = 6;
    hr = 7;
    minute = 51;
    second = 28.386009;
    [jd, jdFrac] = jday(year, mon, day, hr, minute, second);

    dut1 = -0.4399619;      % second
    dat = 32;               % second
    xp = -0.140682 * conv;  % ' to rad
    yp = 0.333309 * conv;
    lod = 0.0015563;
    ddpsi = -0.052195 * conv;  % ' to rad
    ddeps = -0.003875 * conv;
    ddx = -0.000205 * conv;    % ' to rad
    ddy = -0.000136 * conv;
    eqeterms = 2;
    aecef = [ 0.0; 0.0; 0.0 ];

    % note you have to use tdb for time of ineterst AND j2000 (when dat = 32)
    ttt = (jd + jdFrac + (dat + 32.184) / 86400.0 - 2451545.0 - (32 + 32.184) / 86400.0) / 36525.0;
    jdut1 = jd + jdFrac + dut1 / 86400.0;

    fprintf(fid,'ITRF          IAU-76/FK5  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', recef(1), recef(2), recef(3), vecef(1), vecef(2), vecef(3));

    [rteme, vteme, ateme] = ecef2teme(recef, vecef, aecef, ttt, jdut1, lod, xp, yp,    eqeterms);
    fprintf(fid,'TEME          IAU-76/FK5  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', rteme(1), rteme(2), rteme(3), vteme(1), vteme(2), vteme(3));

    vecef = [ 0.0; 0.0; 0.0 ];
    [recef, vecef, aecef] = teme2ecef(rteme, vteme, ateme, ttt, jdut1, lod, xp, yp,    eqeterms);
    fprintf(fid,'ITRF          IAU-76/FK5  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', recef(1), recef(2), recef(3), vecef(1), vecef(2), vecef(3));

end


function testteme_eci(fid)
    conv = pi / (180.0 * 3600.0);

    reci = [ 5102.5089579; 6123.0114007; 6378.1369282 ];
    veci = [ -4.743220157; 0.790536497; 5.533755727 ];
    aeci = [0.0; 0.0; 0.0];
    year = 2004;
    mon = 4;
    day = 6;
    hr = 7;
    minute = 51;
    second = 28.386009;
    [jd jdFrac] = jday(year, mon, day, hr, minute, second);

    dut1 = -0.4399619;      % second
    dat = 32;               % second
    xp = -0.140682 * conv;  % ' to rad
    yp = 0.333309 * conv;
    ddpsi = -0.052195 * conv;  % ' to rad
    ddeps = -0.003875 * conv;
    ddx = -0.000205 * conv;    % ' to rad
    ddy = -0.000136 * conv;

    fileLoc = 'D:\Codes\LIBRARY\DataLib\';
    [iau80arr] = iau80in(fileLoc);

    % note you have to use tdb for time of ineterst AND j2000 (when dat = 32)
    ttt = (jd + jdFrac + (dat + 32.184) / 86400.0 - 2451545.0 - (32 + 32.184) / 86400.0) / 36525.0;
    jdut1 = jd + jdFrac + dut1 / 86400.0;

    fprintf(fid,'GCRF          IAU-76/FK5  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', reci(1), reci(2), reci(3), veci(1), veci(2), veci(3));
    [rteme, vteme, ateme] = eci2teme  ( reci, veci, aeci, iau80arr, ttt, ddpsi, ddeps);
    fprintf(fid,'TEME          IAU-76/FK5  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', rteme(1), rteme(2), rteme(3), vteme(1), vteme(2), vteme(3));

    [reci, veci, aeci] = teme2eci(rteme, vteme, ateme, iau80arr, ttt, ddpsi, ddeps);
    fprintf(fid,'GCRF          IAU-76/FK5  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', reci(1), reci(2), reci(3), veci(1), veci(2), veci(3));
end


function testqmod2ecef(fid)
    ttt = 0.042623631889;
    jdutc = 2453101.82740678310;

    [fArgs] = fundarg(ttt, '80');

    rqmod = [ 5102.5089579; 6123.0114007; 6378.1369282 ];
    vqmod = [ -4.743220157; 0.790536497; 5.533755727 ];

    fileLoc = 'D:\Codes\LIBRARY\DataLib\';
    [iau80arr] = iau80in(fileLoc);

   % [recef, vecef] = qmod2ecef(rqmod, vqmod, iau80arr, ttt, jdutc);
end

function testcsm2efg(fid)
    xp = 0.0;
    yp = 0.0;
    lod = 0.0;
    jdut1 = 2453101.82740678310;
    ttt = 0.042623631889;
    ddpsi = -0.052195;
    ddeps = -0.003875;
    eqeterms = 2;

    % [r1ecef, v1ecef, r2ecef, v2ecef] = csm2efg(r1pef, v1pef, r2ric, v2ric, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps);
    % 
    % fprintf(fid,'csm2efg  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', r1ecef(1), r1ecef(2), r1ecef(3), r2ecef(1), r2ecef(2), r2ecef(3));
end

function testrv_elatlon(fid)
    rad = 180.0 / pi;
    rr = 12756.00;
    ecllat = 60.04570;
    ecllon = 256.004567345;
    drr = 0.045670;
    decllat = 6.798614;
    decllon = 0.00768;

    [rijk, vijk] = elatlon2rv(rr, ecllat, ecllon, drr, decllat, decllon);
    fprintf(fid,'rv ecllat  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', rijk(1), rijk(2), rijk(3), vijk(1), vijk(2), vijk(3));

    [rr, ecllat, ecllon, drr, decllat, decllon] = rv2elatlon(rijk, vijk);
    fprintf(fid,'ecl lat lon  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', rr, ecllat * rad, ecllon * rad, drr, decllat, decllon);



    % additional tests
    rad    = 180.0 / pi;
    twopi = 2.0 * pi;
    % 1" to rad
    convrt = pi / (3600.0*180.0);

    latgd = 20.7071/rad;
    lon = -156.257/rad;
    alt = 3.073;  % km
    [rsecef, vsecef] = site ( latgd,lon,alt );
    fprintf(fid,'\n\nrsecef = -5466.080829    -2404.282897    2242.177454\n');
    fprintf(fid,'rsecef  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', rsecef(1), rsecef(2), rsecef(3), vsecef(1), vsecef(2), vsecef(3));

    year = 2021;
    mon  =  10;
    day  =  12;
    hr   =    4;
    minute  =   10;
    second  =   0.00;
    [jd, jdfrac] =  jday( year, mon, day, hr, minute, second );

    % 2021 10 11 59498  0.204294  0.265940 -0.1059506 -0.0002003 -0.116845 -0.008499  0.000164 -0.000071  37
    % 2021 10 12 59499  0.202929  0.265318 -0.1056985 -0.0003003 -0.116696 -0.008264  0.000206 -0.000037  37
    % 2021 10 13 59500  0.201347  0.265115 -0.1053769 -0.0003296 -0.116533 -0.008163  0.000249 -0.000003  37
    % 2021 10 14 59501  0.199531  0.264474 -0.1050561 -0.0002785 -0.116392 -0.008233  0.000253 -0.000031  37

    dut1 =  -0.1056985;  % s
    dat  = 37;
    xp   =  0.202929* convrt;  % " to rad
    yp   =  0.265318 * convrt;
    lod  =  -0.0003003;
    timezone = 0;
    ddpsi =  -0.116696* convrt;
    ddeps =  -0.008264* convrt;
    dx  = 0.000206* convrt;
    dy  = -0.000037* convrt;
    order =  106;
    terms = 2;
    % , tcg, jdtcg,jdtcgfrac, tcb, jdtcb,jdtcbfrac
    [ ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb, jdtdb, jdtdbfrac ] ...
        = convtime ( year, mon, day, hr, minute, second, timezone, dut1, dat );
    fprintf(fid,'ut1 %8.6f tut1 %16.12f jdut1 %18.11f\n',ut1,tut1,jdut1 );
    fprintf(fid,'utc %8.6f\n',utc );
    fprintf(fid,'tai %8.6f\n',tai );
    fprintf(fid,'tt  %8.6f ttt  %16.12f jdtt  %18.11f\n',tt,ttt,jdtt );
    fprintf(fid,'tdb %8.6f ttdb %16.12f jdtdb %18.11f\n',tdb,ttdb,jdtdb );

    [lst, gst] = lstime(lon, jd + jdfrac);
    fprintf(fid,'lst %15.11f gst %15.11f\n',lst*rad, gst*rad );

    %     rho = 100000.0;
    %     az = 40.0/rad;
    %     el = 20.0/rad;
    %
    %     % horizontal
    %     [r, v] = razel2rv ( rho,az,el,0.0, 0.0, 0.0,latgd,lon,alt, ttt,jdut1,lod,xp,yp,terms,ddpsi,ddeps );
    fileLoc = 'D:\Codes\LIBRARY\DataLib\';
    [iau80arr] = iau80in(fileLoc);
    
    reci = [2919.71566515;    -6559.47300411;      276.48177946];
    veci = [-1.168779561;    -0.198323404;    7.352918872];
    aeci = [0.0; 0.0; 0.0];
    fprintf(fid,'reci    %14.7f %14.7f %14.7f',reci );
    fprintf(fid,' v %14.9f %14.9f %14.9f\n',veci );

   % [rho, az, el, drho, daz, del] = rv2razel(recef, vecef, latgd, lon, alt);
   %  if az < 0.0
   %      az = az + twopi;
   %  end
   %  fprintf(fid,'rvraz   %14.7f %14.7f %14.7f',rho, az*rad, el*rad );
   %  fprintf(fid,' %14.7f %14.12f %14.12f\n',drho, daz*rad, del*rad );
   %  fprintf(fid,'STK 12 Oct 2021 04:09:07.155          159.523              5.000    2788.517174\n');
   %  fprintf(fid,'STK 12 Oct 2021 04:10:07.000          158.339              9.710    2393.490995\n');
   %  %  rtascdecl report
   %  %             147.238               57.942    12 Oct 2021 04:09:00.000
   %  %             144.255               53.546    12 Oct 2021 04:10:00.000
   % 
   %  %rtascdav report
   %  %12 Oct 2021 04:09:07.154735    2788.517174       159.522583           5.000332           326.863806           -57.466141
   %  %12 Oct 2021 04:10:07.000000    2393.490995       158.338742           9.710252           323.927386           -52.964440
   % 
   %  % geocentric
   %  [rr,rtasc,decl,drr,drtasc,ddecl] = rv2radec( reci, veci );
   %  fprintf(fid,'            rho km       rtasc deg     decl deg      drho km/s      drtasc deg/s   ddecl deg/s\n' );
   %  if rtasc < 0.0
   %      rtasc = rtasc + twopi;
   %  end
   %  fprintf(fid,'radec  %14.7f %14.7f %14.7f',rr,rtasc*rad,decl*rad );
   %  fprintf(fid,' %14.7f %14.12f %14.12f\n',drr,drtasc*rad,ddecl*rad );
   % 
   %  [reci,veci] = radec2rv(rr,rtasc,decl,drr,drtasc,ddecl);
   %  fprintf(fid,'reci    %14.7f %14.7f %14.7f',reci );
   %  fprintf(fid,' v %14.9f %14.9f %14.9f\n',veci );

    % topocentric
    % ddpsi = 0.0;
    % ddeps = 0.0;
    % [rho, trtasc, tdecl, drho, dtrtasc, dtdecl] = rv2tradec ( recef, vecef, rsecef, vsecef );
    % fprintf(fid,'           trho km      trtasc deg    tdecl deg     tdrho km/s     tdrtasc deg/s  tddecl deg/s\n' );
    % if trtasc < 0.0
    %     trtasc = trtasc + twopi;
    % end
    % fprintf(fid,'tradec  %14.7f %14.7f %14.7f',trr,trtasc*rad,tdecl*rad );
    % fprintf(fid,' %14.7f %14.12f %14.12f\n',tdrr,tdrtasc*rad,tddecl*rad );
    % fprintf(fid,'STK 12 Oct 2021 04:09:07.155              326.864              -57.466\n');
    % fprintf(fid,'STK 12 Oct 2021 04:10:07.000              323.927              -52.964\n');
    % 

    % satellite ecef
    % 12 Oct 2021 04:10:00.0000    -6166.715106346013    -3676.915653933220      282.460443869915
    %                              -0.605685283644995     1.601875853811738     7.350462613646496
    % satellite eci vector
    % 12 Oct 2021 04:10:00.0000     2919.71566515    -6559.47300411      276.48177946
    %                               -1.168779561    -0.198323404     7.352918872


end


function testrv2radec(fid)
    r = [ -605.79221660; -5870.22951108; 3493.05319896 ];
    v = [ -1.56825429; -3.70234891; -6.47948395 ];
    rr = 0.0;
    rtasc = 0.0;
    decl = 0.0;
    drr = 0.0;
    drtasc = 0.0;
    ddecl = 0.0;

    [rr, rtasc, decl, drr, drtasc, ddecl] = rv2radec(r, v);
    fprintf(fid,'rv radec  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', r(1), r(2), r(3), v(1), v(2), v(3));
    fprintf(fid,'rhosez  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', rr, rtasc, decl, drr, drtasc, ddecl);

    [r, v] = radec2rv(rr, rtasc, decl, drr, drtasc, ddecl);
    fprintf(fid,'rv radec  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', r(1), r(2), r(3), v(1), v(2), v(3));
end


function testrv_razel(fid)
    rad = 180.0 / pi;
    recef = [ -605.79221660; -5870.22951108; 3493.05319896 ];
    vecef = [ -1.56825429; -3.70234891; -6.47948395 ];
    %rsecef = [ -1605.79221660; -570.22951108; 193.05319896 ];
    lon = -104.883 / rad;
    latgd = 39.883 / rad;
    alt = 2.102;
    rho = 0.0186569;
    az = -0.3501725;
    el = -0.5839385;
    drho = 0.6811410;
    daz = -0.4806057;
    del = 0.6284403;

    [rho, az, el, drho, daz, del] = rv2razel(recef, vecef, latgd, lon, alt);
    fprintf(fid,'rv razel  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', recef(1), recef(2), recef(3), vecef(1), vecef(2), vecef(3));
    fprintf(fid,'rhosez  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', rho, az, el, drho, daz, del);

    [recef, vecef] = razel2rv(latgd, lon, alt, rho, az, el, drho, daz, del);
    fprintf(fid,'rv razel  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', recef(1), recef(2), recef(3), vecef(1), vecef(2), vecef(3));
end

function testrv_tradec(fid)
    rijk = [ 4066.716; -2847.545; 3994.302 ];
    vijk = [ -1.56825429; -3.70234891; -6.47948395 ];
    rsijk = [ -1605.79221660; -570.22951108; 193.05319896 ];
    vsijk = [0.0; 0.0; 0.0];
    rho = 0.2634728;
    trtasc = -0.1492353;
    tdecl = 0.0519525;
    drho = 0.3072265;
    dtrtasc = 0.2045751;
    dtdecl = -0.7510033;

    [rho, trtasc, tdecl, drho, dtrtasc, dtdecl] = rv2tradec(rijk, vijk, rsijk, vsijk);
    fprintf(fid,'rv tradec  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', rijk(1), rijk(2), rijk(3), vijk(1), vijk(2), vijk(3));
    fprintf(fid,'tradec %15.11f  %15.11f  %15.11f  %15.11f  %15.11f  %15.11f\n', rho, trtasc, tdecl, drho, dtrtasc, dtdecl);

    [rijk, vijk] = tradec2rv(rho, trtasc, tdecl, drho, dtrtasc, dtdecl, rsijk, vsijk);
    fprintf(fid,'rv tradec  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', rijk(1), rijk(2), rijk(3), vijk(1), vijk(2), vijk(3));
end


function testrvsez_razel(fid)
    rhosez = [ -605.79221660; -5870.22951108; 3493.05319896 ];
    drhosez = [ -1.56825429; -3.70234891; -6.47948395 ];
    rho = 0.0186569;
    az = -0.3501725;
    el = -0.5839385;
    drho = 0.6811410;
    daz = -0.4806057;
    del = 0.6284403;
    % 
    % [rho, az, el, drho, daz, del] = rvsez2razel(rhosez, drhosez);
    % fprintf(fid,'rv rhosez  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', rhosez(1), rhosez(2), rhosez(3), drhosez(1), drhosez(2), drhosez(3));
    % fprintf(fid,'rhosez  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', rho, az, el, drho, daz, del);
    % 
    % [rhosez, drhosez, ] = razel2rvsez( rho, az, el, drho, daz, del);
    % fprintf(fid,'rv rhosez  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', rhosez(1), rhosez(2), rhosez(3), drhosez(1), drhosez(2), drhosez(3));
end


function testrv2rsw(fid)
    r = [ -605.79221660; -5870.22951108; 3493.05319896 ];
    v = [ -1.56825429; -3.70234891; -6.47948395 ];

    [rrsw, vrsw, tm] = rv2rsw( r, v );

    fprintf(fid,'rv2rsw   %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', rrsw(1), rrsw(2), rrsw(3), vrsw(1), vrsw(2), vrsw(3));
end

function testrv2pqw(fid)
    r = [ -605.79221660; -5870.22951108; 3493.05319896 ];
    v = [ -1.56825429; -3.70234891; -6.47948395 ];

    [rpqw, vpqw] = rv2pqw( r, v );

    fprintf(fid,'rv2pqw   %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', rpqw(1), rpqw(2), rpqw(3), vpqw(1), vpqw(2), vpqw(3));
end

function testrv2coe(fid)
    constastro;
    for i=1:21

        if (i == 1)

            r = [ -605.79221660; -5870.22951108; 3493.05319896 ];
            v = [ -1.56825429; -3.70234891; -6.47948395 ];
        end
        if (i == 2)

            fprintf(fid,' coe test ----------------------------\n');
            r = [ 6524.834; 6862.875; 6448.296 ];
            v = [ 4.901327; 5.533756; -1.976341 ];
        end

        % ------- elliptical orbit tests -------------------
        if (i == 3)

            fprintf(fid,' coe test elliptical ----------------------------\n');
            r = [ 1.1372844 * re; -1.0534274 * re; -0.8550194 * re ];
            v = [ 0.6510489 * velkmps; 0.4521008 * velkmps; 0.0381088 * velkmps ];
        end
        if (i == 4)

            fprintf(fid,' coe test elliptical ----------------------------\n');
            r = [ 1.056194 * re; -0.8950922 * re; -0.0823703 * re ];
            v = [ -0.5981066 * velkmps; -0.6293575 * velkmps; 0.1468194 * velkmps ];
        end

        % ------- circular inclined orbit tests -------------------
        if (i == 5)

            fprintf(fid,' coe test near circular inclined ----------------------------\n');
            r = [ -0.422277 * re; 1.0078857 * re; 0.7041832 * re ];
            v = [ -0.5002738 * velkmps; -0.5415267 * velkmps; 0.4750788 * velkmps ];
        end
        if (i == 6)

            fprintf(fid,' coe test near circular inclined ----------------------------\n');
            r = [ -0.7309361 * re; -0.6794646 * re; -0.8331183 * re ];
            v = [ -0.6724131 * velkmps; 0.0341802 * velkmps; 0.5620652 * velkmps ];
        end

        if (i == 7) % -- CI u = 45 deg

            fprintf(fid,' coe test circular inclined ----------------------------\n');
            r = [ -2693.34555010128; 6428.43425355863; 4491.37782050409 ];
            v = [ -3.95484712246016; -4.28096585381370; 3.75567104538731 ];
        end
        if (i == 8) % -- CI u = 315 deg

            fprintf(fid,' coe test circular inclined ----------------------------\n');
            r = [ -7079.68834483379; 3167.87718823353; -2931.53867301568 ];
            v = [ 1.77608080328182; 6.23770933190509; 2.45134017949138 ];
        end

        % ------- elliptical equatorial orbit tests -------------------
        if (i == 9)

            fprintf(fid,' coe test elliptical near equatorial ----------------------------\n');
            r = [ 21648.6109280739; -14058.7723188698; -0.0003598029 ];
            v = [ 2.16378060719980; 3.32694348486311; 0.00000004164788 ];
        end
        if (i == 10)

            fprintf(fid,' coe test elliptical near equatorial ----------------------------\n');
            r = [ 7546.9914487222; 24685.1032834356; -0.0003598029 ];
            v = [ 3.79607016047138; -1.15773520476223; 0.00000004164788 ];
        end

        if (i == 11) % -- EE w = 20 deg

            fprintf(fid,' coe test elliptical equatorial ----------------------------\n');
            r = [ -22739.1086596208; -22739.1086596208; 0.0 ];
            v = [ 2.48514004188565; -2.02004112073465; 0.0 ];
        end
        if (i == 12) % -- EE w = 240 deg

            fprintf(fid,' coe test elliptical equatorial ----------------------------\n');
            r = [ 28242.3662822040; 2470.8868808397; 0.0 ];
            v = [ 0.66575215057746; -3.62533022188304; 0.0 ];
        end

        % ------- circular equatorial orbit tests -------------------
        if (i == 13)

            fprintf(fid,' coe test circular near equatorial ----------------------------\n');
            r = [ -2547.3697454933; 14446.8517254604; 0.000 ];
            v = [ -5.13345156333487; -0.90516601477599; 0.00000090977789 ];
        end
        if (i == 14)

            fprintf(fid,' coe test circular near equatorial ----------------------------\n');
            r = [ 7334.858850000; -12704.3481945462; 0.000 ];
            v = [ -4.51428154312046; -2.60632166411836; 0.00000090977789 ];
        end

        if (i == 15) % -- CE l = 65 deg

            fprintf(fid,' coe test circular equatorial ----------------------------\n');
            r = [ 6199.6905946008; 13295.2793851394; 0.0 ];
            v = [ -4.72425923942564; 2.20295826245369; 0.0 ];
        end
        if (i == 16) % -- CE l = 65 deg i = 180 deg

            fprintf(fid,' coe test circular equatorial ----------------------------\n');
            r = [ 6199.6905946008; -13295.2793851394; 0.0 ];
            v = [ -4.72425923942564; -2.20295826245369; 0.0 ];
        end

        % ------- parabolic orbit tests -------------------
        if (i == 17)

            fprintf(fid,' coe test parabolic ----------------------------\n');
            r = [ 0.5916109 * re; -1.2889359 * re; -0.3738343 * re ];
            v = [ 1.1486347 * velkmps; -0.0808249 * velkmps; -0.1942733 * velkmps ];
        end

        if (i == 18)

            fprintf(fid,' coe test parabolic ----------------------------\n');
            r = [ -1.0343646 * re; -0.4814891 * re; 0.1735524 * re ];
            v = [ 0.1322278 * velkmps; 0.7785322 * velkmps; 1.0532856 * velkmps ];
        end

        if (i == 19)

            fprintf(fid,' coe test hyperbolic ---------------------------\n');
            r = [ 0.9163903 * re; 0.7005747 * re; -1.3909623 * re ];
            v = [ 0.1712704 * velkmps; 1.1036199 * velkmps; -0.3810377 * velkmps ];
        end

        if (i == 20)

            fprintf(fid,' coe test hyperbolic ---------------------------\n');
            r = [ 12.3160223 * re; -7.0604653 * re; -3.7883759 * re ];
            v = [ -0.5902725 * velkmps; 0.2165037 * velkmps; 0.1628339 * velkmps ];
        end

        if (i == 21)

            fprintf(fid,' coe test rectilinear --------------------------\n');
            r = [ -1984.03023322569; 1525.27235370582; 6364.76955283447 ];
            v = [ -1.60595491095; 1.23461759098; 5.15190381139 ];
            v = [-1.60936089585; 1.23723602618; 5.16283021192];  % 196
        end

        fprintf(fid,' start r   %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', r(1), r(2), r(3), v(1), v(2), v(3));

        % --------  coe2rv       - classical elements to posisiotn and velocity
        % --------  rv2coe       - position and velocity vectors to classical elements
        [p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper ] = rv2coe (r, v);
        fprintf(fid,'           p km       a km      ecc      incl deg     raan deg     argp deg      nu deg      m deg      arglat   truelon    lonper\n');
        fprintf(fid,'ansr coes  %11.4f %11.4f %13.9f %13.7f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f\n', p, a, ecc, incl * rad, ...
            raan * rad, argp * rad, nu * rad, m * rad, arglat * rad, truelon * rad, lonper * rad);

        % rectilinear orbits have sign(a) determines orbit type, arglat
        % is nu, but the magnitude is off...?
        if (abs(ecc - 1.0) < 0.0000001)
            p = mag(r) * 1.301;
        end
        [r1, v1] = coe2rv(p, ecc, incl, raan, argp, nu, arglat, truelon, lonper);
        fprintf(fid,' end r   %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', r1(1), r1(2), r1(3), v1(1), v1(2), v1(3));
    end  % through for
end  % testrv2coe


function testfindc2c3(fid)
    % --------  findc2c3     - find c2 c3 parameters for f and g battins method
    znew = -39.47842;
    [c2new, c3new] = findc2c3(znew);
    fprintf(fid,'findc2c3 z %15.11f  %15.11f  %15.11f\n', znew, c2new, c3new);

    znew = 0.0;
    [c2new, c3new] = findc2c3(znew);
    fprintf(fid,'findc2c3 z %15.11f  %15.11f  %15.11f\n', znew, c2new, c3new);

    znew = 0.57483;
    [c2new, c3new] = findc2c3(znew);
    fprintf(fid,'findc2c3 z %15.11f  %15.11f  %15.11f\n', znew, c2new, c3new);

    znew = 39.47842;
    [c2new, c3new] = findc2c3(znew);
    fprintf(fid,'findc2c3 z %15.11f  %15.11f  %15.11f\n', znew, c2new, c3new);
end


function testcoe2rv(fid)
    rad = 180.0 /pi;

    % alt test various combinations of coe/eq and rv
    for j = 1:2
        if j == 1
            fprintf(fid,'coe tests ----------------------------\n' );
        else
            fprintf(fid,'\n\neq tests ----------------------------\n' );
            %pause;
        end
        for i = 1:21
            if i == 1
                r=[ 6524.834;6862.875;6448.296];
                v=[ 4.901327;5.533756;-1.976341];
            end
            if i == 2
                fprintf(fid,'coe test ----------------------------\n' );
                r=[ 6524.834;6862.875;6448.296];
                v=[ 4.901327;5.533756;-1.976341];
            end

            % ------- elliptical orbit tests -------------------
            if i == 3
                fprintf(fid,'coe test elliptical ----------------------------\n' );
                r=[ 1.1372844; -1.0534274; -0.8550194]*6378.137;
                v=[0.6510489;  0.4521008;  0.0381088]*7.905366149846;
            end
            if i == 4
                fprintf(fid,'coe test elliptical ----------------------------\n' );
                r=[  1.0561942;-0.8950922;-0.0823703]*6378.137;
                v=[  -0.5981066;-0.6293575; 0.1468194]*7.905366149846;
            end

            % ------- circular inclined orbit tests -------------------
            if i == 5
                fprintf(fid,'coe test near circular inclined ----------------------------\n' );
                r=[ -0.4222777; 1.0078857; 0.7041832]*6378.137;
                v=[  -0.5002738;-0.5415267; 0.4750788]*7.905366149846;
            end
            if i == 6
                fprintf(fid,'coe test near circular inclined ----------------------------\n' );
                r=[ -0.7309361;-0.6794646;-0.8331183]*6378.137;
                v=[  -0.6724131; 0.0341802; 0.5620652]*7.905366149846;
            end

            if i == 7 % -- CI u = 45 deg
                fprintf(fid,'coe test circular inclined ----------------------------\n' );
                r = [-2693.34555010128  6428.43425355863  4491.37782050409];
                v = [   -3.95484712246016  -4.28096585381370  3.75567104538731];
            end
            if i == 8 % -- CI u = 315 deg
                fprintf(fid,'coe test circular inclined ----------------------------\n' );
                r = [-7079.68834483379;  3167.87718823353; -2931.53867301568];
                v = [    1.77608080328182;  6.23770933190509; 2.45134017949138];
            end

            % ------- elliptical equatorial orbit tests -------------------
            if i == 9
                fprintf(fid,'coe test elliptical near equatorial ----------------------------\n' );
                r=[ 21648.6109280739; -14058.7723188698; -0.0003598029];
                v=[ 2.16378060719980; 3.32694348486311; 0.00000004164788 ];
            end
            if i == 10
                fprintf(fid,'coe test elliptical near equatorial ----------------------------\n' );
                r=[  7546.9914487222;  24685.1032834356; -0.0003598029];
                v=[ 3.79607016047138; -1.15773520476223; 0.00000004164788 ];
            end

            if i == 11 % -- EE w = 20 deg
                fprintf(fid,'coe test elliptical equatorial ----------------------------\n' );
                r = [-22739.1086596208;  -22739.1086596208 ;    0.0];
                v = [    2.48514004188565;  -2.02004112073465 ; 0.0];
            end
            if i == 12 % -- EE w = 240 deg
                fprintf(fid,'coe test elliptical equatorial ----------------------------\n' );
                r = [ 28242.3662822040;    2470.8868808397 ;   0.0];
                v = [    0.66575215057746 ; -3.62533022188304 ; 0.0];
            end

            % ------- circular equatorial orbit tests -------------------
            if i == 13
                fprintf(fid,'coe test circular near equatorial ----------------------------\n' );
                r=[ -2547.3697454933; 14446.8517254604; 0.000 ];
                v=[  -5.13345156333487; -0.90516601477599; 0.00000090977789 ];
            end
            if i == 14
                fprintf(fid,'coe test circular near equatorial ----------------------------\n' );
                r=[  7334.858850000; -12704.3481945462;   0.000 ];
                v=[  -4.51428154312046; -2.60632166411836; 0.00000090977789 ];
            end

            if i == 15 % -- CE l = 65 deg
                fprintf(fid,'coe test circular equatorial ----------------------------\n' );
                r = [ 6199.6905946008; 13295.2793851394;      0.0];
                v = [ -4.72425923942564; 2.20295826245369;    0.0];
            end
            if i == 16 % -- CE l = 65 deg i = 180 deg
                fprintf(fid,'coe test circular equatorial ----------------------------\n' );
                r = [ 6199.6905946008; -13295.2793851394;      0.0];
                v = [ -4.72425923942564; -2.20295826245369;    0.0];
            end

            % ------- parabolic orbit tests -------------------
            if i == 17
                fprintf(fid,'coe test parabolic ----------------------------\n' );
                r=[  0.5916109;-1.2889359;-0.3738343]*6378.137;
                v=[   1.1486347;-0.0808249;-0.1942733]*7.905366149846;
            end

            if i == 18
                fprintf(fid,'coe test parabolic ----------------------------\n' );
                r=[-1.0343646; -0.4814891;  0.1735524]*6378.137;
                v=[ 0.1322278; 0.7785322; 1.0532856  ]*7.905366149846;
            end

            if i == 19
                fprintf(fid,'coe test hyperbolic ---------------------------\n' );
                r=[0.9163903; 0.7005747; -1.3909623  ]*6378.137;
                v=[0.1712704; 1.1036199; -0.3810377  ]*7.905366149846;
            end

            if i == 20
                fprintf(fid,'coe test hyperbolic ---------------------------\n' );
                r=[12.3160223; -7.0604653; -3.7883759]*6378.137;
                v=[-0.5902725; 0.2165037; 0.1628339  ]*7.905366149846;
            end

            if i == 21
                fprintf(fid,'coe test rectilinear --------------------------\n' );
                r = [-1984.03023322569; 1525.27235370582; 6364.76955283447];
                v = [-1.60595491095; 1.23461759098; 5.15190381139];  % 201?
                %v = [-1.60936089585; 1.23723602618; 5.16283021192];  % 196
            end

            fprintf(fid,'start %15.9f %15.9f %15.9f',r );
            fprintf(fid,' v  %15.10f %15.10f %15.10f\n',v );

            if j == 1
                % --------  rv2coe       - position and velocity vectors to classical elements
                [p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper ] = rv2coe (r, v);
                fprintf(fid,'          p km         a km         ecc        incl deg     raan deg     argp deg      nu deg      m deg      arglat   truelon    lonper\n');
                fprintf(fid,'coes %11.4f %11.4f %13.9f %13.7f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f\n',...
                    p,a,ecc,incl*rad,raan*rad,argp*rad,nu*rad,m*rad, ...
                    arglat*rad,truelon*rad,lonper*rad );

                % --------  coe2rv       - classical elements to position and velocity
                % rectilinear orbits have sign(a) determines orbit type, arglat
                % is nu, but the magnitude is off...?
                if abs(ecc-1.0) < 0.0000001
                    p = mag(r)*1.301;
                end
                [rn,vn] = coe2rv(p,ecc,incl,raan,argp,nu,arglat,truelon,lonper);
                fprintf(fid,'rn    %15.9f %15.9f %15.9f',rn );
                fprintf(fid,' vn %15.10f %15.10f %15.10f\n',vn );
                dr(1) = r(1) - rn(1);
                dr(2) = r(2) - rn(2);
                dr(3) = r(3) - rn(3);
                if mag(dr) > 0.01
                    fprintf(fid,'ERROR in this case dr = %15.11f\n', mag(dr));
                end
            else
                % --------  rv2eq       - position and velocity vectors to classical elements
                [ a, n, af, ag, chi, psi, meanlonM, meanlonNu, fr ] = rv2eq (r,v);
                fprintf(fid,'       fr     a km         n rad      af           ag         chi          psi      meanlonnu deg   meanlonm deg\n');
                fprintf(fid,'eqs    %2d %11.4f %11.4f %13.9g %13.7g %11.5g %11.5g %11.5f %11.5f\n',...
                    fr, a, n, af, ag, chi, psi, meanlonNu*rad, meanlonM*rad );

                % --------  eq2rv       - classical elements to position and velocity
                [rn,vn] = eq2rv( a, af, ag, chi, psi, meanlonM, fr);
                fprintf(fid,'rn    %15.9f %15.9f %15.9f',rn );
                fprintf(fid,' vn %15.10f %15.10f %15.10f\n',vn );
                dr(1) = r(1) - rn(1);
                dr(2) = r(2) - rn(2);
                dr(3) = r(3) - rn(3);
                if mag(dr) > 0.01
                    fprintf(fid,'ERROR in this case dr = %15.11f\n', mag(dr));
                end
            end

        % reci = [1525.9870698051157; -5867.209915411114; 3499.601587508083]';
        % veci = [1.4830443958075603; -7.093267951700349; 0.9565730381487033]';
        % rmag = 7000; % km
        % vmag = 7.546;  % km/s
        % latgc = pi / 6;  % 30 degrees
        % lon = pi / 2;  % 90 degrees
        % fpa = -pi / 6;  % -30 degrees
        % az = pi / 4;  % 45 degrees
        %
        % conv = pi / (180.0*3600.0);
        % ttt = 0.042623631888994;
        % jdut1 = 2.45310150e+06;
        % lod = 0.0015563;
        % xp = -0.140682 * conv;
        % yp = 0.333309 * conv;
        % eqeterms = 2;
        % ddpsi = -0.052195 * conv;
        % ddeps = -0.003875 * conv;
        % % ---- flight elements
        % [lon, latgc, rtasc, decl, fpa, az, magr, magv] = rv2flt ( reci,veci,ttt,jdut1,lod,xp,yp,2,ddpsi,ddeps );
        % fprintf(fid,'         rmag km       vmag km/s     latgc deg       lon deg       fpa deg       az deg\n');
        % fprintf(fid,'flt  %14.7f%14.7f%14.7f%15.7f%14.7f%14.7f\n',rmag,vmag,...
        %         latgc*rad,lon*rad,fpa*rad,az*rad );
        % [r,v] = flt2rv ( rmag,vmag,latgc,lon,fpa,az,ttt,jdut1,lod,xp,yp,2,ddpsi,ddeps );
        % fprintf(fid,'r    %15.9f%15.9f%15.9f',r );
        % fprintf(fid,' v %15.10f%15.10f%15.10f\n',v );
        %
        % % ----  adbarv elements
        % [rmag,vmag,rtasc,decl,fpav,az] = rv2adbar ( r,v );
        % fprintf(fid,'          rmag km      vmag km/s     rtasc deg       decl deg      fpav deg      az deg\n');
        % fprintf(fid,'adb  %14.7f%14.7f%14.7f%15.7f%14.7f%14.7f\n',rmag,vmag,...
        %          rtasc*rad,decl*rad,fpav*rad,az*rad );
        % [r,v] = adbar2rv ( rmag,vmag,rtasc,decl,fpav,az );
        % fprintf(fid,'r    %15.9f%15.9f%15.9f',r );
        % fprintf(fid,' v %15.10f%15.10f%15.10f\n',v );
        %
        % % ---- radial, along-track, cross-track
        % [rrac,vrac,transrmat] = rv2rac(r,v);
        % fprintf(fid,'rac  %15.9f%15.9f%15.9f',rrac );
        % fprintf(fid,' v %15.10f%15.10f%15.10f\n',vrac );
        %
        % % ---- in-radial, velocity, cross-track
        % [rivc,vivc,transrmat] = rv2ivc(r,v);
        % fprintf(fid,'ivc  %15.9f%15.9f%15.9f',rivc );
        % fprintf(fid,' v %15.10f%15.10f%15.10f\n',vivc );

        end  % for
    end % for through coe/eq tests


    fprintf(fid,'\n\n\n tests\n');
    r = [4942.74746831; 4942.74746831; 0.];
    v = [-5.34339547; 5.34339547; 0.02137362];
    fprintf(fid,'r    %15.9f %15.9f %15.9f',r );
    fprintf(fid,' v %15.10f %15.10f %15.10f\n',v );
    [ a, n, af, ag, chi, psi, meanlonM, meanlonNu, fr ] = rv2eq (r, v);
    fprintf(fid,'       fr     a km         n rad      af           ag         chi          psi      meanlonnu deg   meanlonm deg\n');
    fprintf(fid,'eqs    %2d %11.4f %11.4f %13.9f %13.7f %11.5f %11.5f %11.5f %11.5f\n',...
        fr, a, n, af, ag, chi, psi, meanlonNu*rad, meanlonM*rad );
    [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe (r,v);
    fprintf(fid,'          p km         a km         ecc        incl deg     raan deg     argp deg      nu deg      m deg      arglat   truelon    lonper\n');
    fprintf(fid,'coes %11.4f %11.4f %13.9f %13.7f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f\n',...
        p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad, ...
        arglat*rad,truelon*rad,lonper*rad );

    fprintf(fid,'\n\ STK ? tests\n');
    r = [4942.72769736; -4942.72769736;  19.77095033];
    v = [-5.34341685; -5.34341685; 0];
    fprintf(fid,'r    %15.9f %15.9f %15.9f',r );
    fprintf(fid,' v %15.10f %15.10f %15.10f\n',v );
    [ a, n, af, ag, chi, psi, meanlonM, meanlonNu, fr ] = rv2eq (r, v);
    fprintf(fid,'       fr     a km         n rad      af           ag         chi          psi      meanlonnu deg   meanlonm deg\n');
    fprintf(fid,'eqs    %2d %11.4f %11.4f %13.9f %13.7f %11.5f %11.5f %11.5f %11.5f\n',...
        fr, a, n, af, ag, chi, psi, meanlonNu*rad, meanlonM*rad );
    [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe (r,v);
    fprintf(fid,'          p km         a km         ecc        incl deg     raan deg     argp deg      nu deg      m deg      arglat   truelon    lonper\n');
    fprintf(fid,'coes %11.4f %11.4f %13.9f %13.7f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f\n',...
        p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad, ...
        arglat*rad,truelon*rad,lonper*rad );


    fprintf(fid,'\n other tests\n');
    [rn,vn] = eq2rv( 7000.0, 0.001, 0.001, 0.001, 0.001, 45.0/rad, fr);
    fprintf(fid,'rn    %15.9f %15.9f %15.9f',rn );
    fprintf(fid,' vn %15.10f %15.10f %15.10f\n',vn );
    [ a, n, af, ag, chi, psi, meanlonM, meanlonNu, fr ] = rv2eq (rn,vn);
    fprintf(fid,'       fr     a km         n rad      af           ag         chi          psi      meanlonnu deg   meanlonm deg\n');
    fprintf(fid,'eqs    %2d %11.4f %11.4f %13.9f %13.7f %11.5f %11.5f %11.5f %11.5f\n',...
        fr, a, n, af, ag, chi, psi, meanlonNu*rad, meanlonM*rad );
end % testcoe2rv


function testlon2nu(fid)
      rad = 180.0 / pi;
  jdut1 = 2449470.5;
    incl = 35.324598 / rad;
    lon = -121.3487 / rad;
    raan = 45.0 / rad;
    argp = 34.456798 / rad;

    lon = lon2nu(jdut1, lon, incl, raan, argp);
end


% faster version?
function testnewtonmx(fid)
    rad = 180.0 / pi;
    ecc = 0.4;
    m = 334.566986 / rad;

    % [eccanom, nu] = newtonmx(ecc, m);
    % 
    % fprintf(fid,' newtonmx ecc %15.11f m  %15.11f  eccanom %15.11f  nu %15.11f\n', ecc, m * rad, eccanom * rad, nu * rad);
end

% --------  newtonm      - find eccentric and true anomaly given ecc and mean anomaly
function testnewtonm(fid)
    rad = 180.0 / pi;
    ecc = 0.4;
    eccanom = 334.566986 / rad;
    [m, nu] = newtone(ecc, eccanom);

    fprintf(fid,' newtone ecc %15.11f eccanom  %15.11f  m %15.11f nu %15.11f\n', ecc, eccanom * rad, m * rad, nu * rad);

    ecc = 0.34;
    m = 235.4 / rad;
    [eccanom, nu] = newtonm(ecc, m);
    fprintf(fid,' newtonm ecc %15.11f m  %15.11f  eccanom %15.11f nu %15.11f\n', ecc, m * rad, eccanom * rad, nu * rad);
end


% --------  newtone      - find true and mean anomaly given ecc and eccentric anomaly
function testnewtone(fid)
    rad = 180.0 / pi;
    ecc = 0.34;
    eccanom = 334.566986 / rad;
    [m, nu] = newtone(ecc, eccanom);

    fprintf(fid,' newtone ecc %15.11f eccanom  %15.11f  m %15.11f nu %15.11f\n', ecc, eccanom * rad, m * rad, nu * rad);
end

% --------  newtonnu     - find eccentric and mean anomaly given ecc and true anomaly
function testnewtonnu(fid)
    rad = 180.0 / pi;
    ecc = 0.34;
    nu = 134.567001 / rad;

    [eccanom, m] = newtonnu(ecc, nu);

    fprintf(fid,' newtonnu ecc %15.11f nu  %15.11f  eccanom %15.11f  m %15.11f\n', ecc, nu * rad, eccanom * rad, m * rad);
end


%
% [r2, v2, c2, c3, x, z] = keplerc2c3(r1, v1, dtsec)
%
function [r2, v2, c2, c3, x, z] = keplerc2c3(r1, v1, dtsec)
    constastro;

    % -------------------------  implementation   -----------------
    show = 'n';
    c2 = 0.0;
    c3 = 0.0;
    x = 0.0;
    z = 0.0;

    for ii = 1:3
        rx(ii) = 0.0;
        vx(ii) = 0.0;
    end

    r2 = rx;  % seems to be the only way to get these variables out
    v2 = vx;
    xnew = 0.0;
    c2new = 0.0;
    c3new = 0.0;

    small = 0.000000001;
    twopi = 2.0 * pi;
    halfpi = pi * 0.5;

    % -------------------------  implementation   -----------------
    % set constants and intermediate printouts
    numiter = 50;

    if (show == 'y')
        %            printf(' r1 %16.8f %16.8f %16.8f ER\n',r1(1)/re,r1(2)/re,r1(3)/ );
        %            printf(' vo %16.8f %16.8f %16.8f ER/TU\n',vo(1)/velkmps, vo(2)/velkmps, vo(3)/velkmps );
    end

    % --------------------  initialize values   -------------------
    ktr = 0;
    xold = 0.0;
    znew = 0.0;
    %errork = '      ok';
    dtseco = dtsec;
    mulrev = 0;

    if (abs(dtseco) > small)

        magro = mag(r1);
        magvo = mag(v1);
        rdotv = dot(r1, v1);

        % -------------  find sme, alpha, and a  ------------------
        sme = ((magvo * magvo) * 0.5) - (mu / magro);
        alpha = -sme * 2.0 / mu;

        if (abs(sme) > small)
            a = -mu / (2.0 * sme);
        else
            a = 999999.9;
            if (abs(alpha) < small)   % parabola
                alpha = 0.0;
            end
        end
        if (show == 'y')

            %           printf(' sme %16.8f  a %16.8f alp  %16.8f ER\n',sme/(mu/), a/re, alpha *  );
            %           printf(' sme %16.8f  a %16.8f alp  %16.8f km\n',sme, a, alpha );
            %           printf(' ktr      xn        psi           r2          xn+1        dtn\n' );
        end

        % ------------   setup initial guess for x  ---------------
        % -----------------  circle and ellipse -------------------
        if (alpha >= small)

            period = twopi * sqrt(abs(a * a * a) / mu);
            % ------- next if needed for 2body multi-rev ----------
            if (abs(dtseco) > abs(period))
                % including the truncation will produce vertical lines that are parallel
                % (plotting chi vs time)
                %                    dtsec = rem( dtseco,period );
                mulrev = int(dtseco / period);
            end
            if (abs(alpha - 1.0) > small)
                xold = sqrt(mu) * dtsec * alpha;
            else
                % - first guess can't be too close. ie a circle, r2=a
                xold = sqrt(mu) * dtsec * alpha * 0.97;
            end
        else

            % --------------------  parabola  ---------------------
            if ( abs( alpha ) < small )
                h = cross( ro,vo );
                magh = mag(h);
                p= magh*magh/mu;
                s= 0.5  * (halfpi - atan( 3.0 *sqrt( mu / (p*p*p) )* dtsec ) );
                w= atan( tan( s )^(1.0 /3.0 ) );
                xold = sqrt(p) * ( 2.0 *cot(2.0 *w) );
                alpha= 0.0;
            else
                % ------------------  hyperbola  ------------------
                temp= -2.0 * mu * dtsec / ...
                    ( a*( rdotv + sign(dtsec)*sqrt(-mu*a)* ...
                    (1.0 -magro*alpha) ) );
                xold= sign(dtsec) * sqrt(-a) *log(temp);
            end
        end % if alpha

        ktr = 1;
        dtnew = -10.0;
         tmp = 1.0 / sqrt(mu);
        while ((abs(dtnew * tmp - dtsec) >= small) && (ktr < numiter))

            xoldsqrd = xold * xold;
            znew = xoldsqrd * alpha;

            % ------------- find c2 and c3 functions --------------
            [c2new, c3new] = findc2c3(znew);

            % ------- use a newton iteration for new values -------
            rval = xoldsqrd * c2new + rdotv * tmp * xold * (1.0 - znew * c3new) + ...
            magro * (1.0 - znew * c2new);
            dtnew = xoldsqrd * xold * c3new + rdotv * tmp * xoldsqrd * c2new + ...
            magro * xold * (1.0 - znew * c3new);

            % ------------- calculate new value for x -------------
            xnew = xold + (dtsec * sqrt(mu) - dtnew) / rval;

            % ----- check if the univ param goes negative. if so, use bissection
            if (xnew < 0.0)
                xnew = xold * 0.5;
            end
            if (show == 'y')

                %  printf('%3i %15.11f %15.11f %15.11f %15.11f %15.11f\n', ktr,xold,znew,rval,xnew,dtnew);
                %  printf('%3i %15.11f %15.11f %15.11f %15.11f %15.11f\n', ktr,xold/sqrt(fid),znew,rval/re,xnew/sqrt(fid),dtnew/sqrt(mu));
            end

            ktr = ktr + 1;
            xold = xnew;
        end  % while

        if (ktr >= numiter)

            %errork = 'knotconv';
            %           printf('not converged in %2i iterations\n',numiter );
            for i=1:3
                v2(i) = 0.0;
                r2(i) = v2(i);
            end
        else

            % --- find position and velocity vectors at new time --
            xnewsqrd = xnew * xnew;
            f = 1.0 - (xnewsqrd * c2new / magro);
            g = dtsec - xnewsqrd * xnew * c3new / sqrt(mu);

            for i=1:3
                r2(i) = f * r1(i) + g * v1(i);
                magr = mag(r2);
                gdot = 1.0 - (xnewsqrd * c2new / magr);
                fdot = (sqrt(mu) * xnew / (magro * magr)) * (znew * c3new - 1.0);
            end
            for i=1:3
                v2(i) = fdot * r1(i) + gdot * v1(i);
                magv = mag(v2);
                temp = f * gdot - fdot * g;
                %if (abs(temp - 1.0) > 0.00001)
                %    errork = 'fandg';
            end
            if (show == 'y')

                %           printf('f %16.8f g %16.8f fdot %16.8f gdot %16.8f\n',f, g, fdot, gdot );
                %           printf('f %16.8f g %16.8f fdot %16.8f gdot %16.8f\n',f, g, fdot, gdot );
                %           printf('r1 %16.8f %16.8f %16.8f ER\n',r2(1)/re,r2(2)/re,r2(3)/ );
                %           printf('v1 %16.8f %16.8f %16.8f ER/TU\n',v(1)/velkmps, v(2)/velkmps, v(3)/velkmps );
            end
        end % if fabs
    else
        % ----------- set vectors to incoming since 0 time --------
        for i=1:3

            r2(i) = r1(i);
            v2(i) = v1(i);
        end

        %       fprintf( fid,'%11.5f  %11.5f %11.5f  %5i %3i ',znew, dtseco/60.0, xold/(rad), ktr, mulrev );
    end
end  % keplerc2c3


function testfindfandg(fid)
    r1 = [ 4938.49830042171; -1922.24810472241; 4384.68293292613 ];
    v1 = [ 0.738204644165659; 7.20989453238397; 2.32877392066299 ];
    r2 = [ -1105.78023519582; 2373.16130661458; 6713.89444816503 ];
    v2 = [ 5.4720951867079; -4.39299050886976; 2.45681739563752 ];
    dtsec = 6372.69272563561; % 1ld
    dtsec = 60.0; % must be small step sizes!!

    %% dan hyperbolic test
    %r1 = [ 4070.5942270000000; 3786.8271570000002; 4697.0576309999997 ];
    %%v1 = [ -32553.559100851671; -37563.543526937596; -37563.543526937596 ];
    %% exact opposite from r velocity
    %v1 = [ -34845.69531184976; -32416.55098811211; -40208.43885307875 ];
    %r2 = [ -1105.78023519582; 2373.16130661458; 6713.89444816503 ];
    %v2 = [ 5.4720951867079; -4.39299050886976; 2.45681739563752 ];
    %dtsec = 0.25; % must be small step sizes!!

    fprintf(fid,'r1  %15.11f  %15.11f  %15.11f v1 %15.11f  %15.11f  %15.11f\n', r1(1), r1(2), r1(3), v1(1), v1(2), v1(3));

    for i = 0:5

        if (i == 0)
            dtsec = 60.0;
        end
        if (i == 1)
            dtsec = 0.1;
        end
        if (i == 2)
            dtsec = 1.0;
        end
        if (i == 3)
            dtsec = 10.0;
        end
        if (i == 4)
            dtsec = 100.0;
        end
        if (i == 5)
            dtsec = 500.0;
        end

        [r2, v2, c2, c3, x, z] = keplerc2c3(r1, v1, dtsec);
        fprintf(fid,'r2  %15.11f  %15.11f  %15.11f v2 %15.11f  %15.11f  %15.11f\n', r2(1), r2(2), r2(3), v2(1), v2(2), v2(3));
        fprintf(fid,'c2 %11.7f c3 %15.11f x %15.11f z %15.11f dtsec %15.11f\n', c2, c3, x, z,dtsec);

        opt = 'pqw';
        [f, g, fdot, gdot] = findfandg(r1, v1, r2, v2, dtsec, x, c2, c3, z, opt);
         ansr = f * gdot - g * fdot;
        fprintf(fid,'f and g pqw  %15.11f %15.11f %15.11f %15.11f %15.11f\n', f, g, fdot, gdot, ansr);

        opt = 'series';  %  pqw, series, c2c3
        [f, g, fdot, gdot] = findfandg(r1, v1, r2, v2, dtsec, x, c2, c3, z, opt);
        ansr = f * gdot - g * fdot;
        fprintf(fid,'f and g series  %15.11f %15.11f %15.11f %15.11f %15.11f\n', f, g, fdot, gdot, ansr);

        opt = 'c2c3';  %  pqw, series, c2c3
        [f, g, fdot, gdot] = findfandg(r1, v1, r2, v2, dtsec, x, c2, c3, z, opt);
        ansr = f * gdot - g * fdot;
        fprintf(fid,'f and g c2c3  %15.11f %15.11f %15.11f %15.11f %15.11f\n', f, g, fdot, gdot, ansr);
    end
end


function testcheckhitearth(fid)
    constastro;

    nrev = 0;
    r1 = [ 2.500000 * re; 0.000000; 0.000000 ];
    r2 = [ 1.9151111 * re; 1.6069690 * re; 0.000000 ];
    % assume circular initial orbit for vel calcs
    v1t = [ 0.0; sqrt(mu / r1(1)); 0.0 ];
    ang = atan(r2(2) / r2(1));
    v2t = [ -sqrt(mu / r2(2)) * cos(ang); sqrt(mu / r2(1)) * sin(ang); 0.0 ];
    altpad = 100.0; % km

    magr1 = mag(r1);
    magr2 = mag(r2);
    cosdeltanu = dot(r1, r2) / (magr1 * magr2);

    [hitearth, hitearthstr] = checkhitearth ( altpad, r1, v1t, r2, v2t, nrev );

    fprintf(fid,'hitearth? %s, %15.11f\n', hitearthstr, (cos(cosdeltanu) * 180.0 / pi));
end

function testcheckhitearthc(fid)
    constastro;

    nrev = 0;
    r1c = [ 2.500000; 0.000000; 0.000000 ];
    r2c = [ 1.9151111; 1.6069690; 0.000000 ];
    % assume circular initial orbit for vel calcs
    v1tc = [ 0.0; sqrt(1.0 / r1c(1)); 0.0 ];
    ang = atan(r2c(2) / r2c(1));
    v2tc = [ -sqrt(1.0 / r2c(2)) * cos(ang); sqrt(1.0 / r2c(1)) * sin(ang); 0.0 ];
    altpadc = 100.0 / re; % er

    magr1c = mag(r1c);
    magr2c = mag(r2c);
    cosdeltanu = dot(r1c, r2c) / (magr1c * magr2c);
    [hitearth, hitearthstr] = checkhitearthc(altpadc, r1c, v1tc, r2c, v2tc, nrev);

    fprintf(fid,'hitearth? %s, %15.11f\n', hitearthstr, (cos(cosdeltanu) * 180.0 / pi));
end


function testgibbs(fid)
    constastro;
    rad = 180.0 / pi;

    r1 = [ 0.0000000; 0.000000; re ];
    r2 = [ 0.0000000; -4464.696; -5102.509 ];
    r3 = [ 0.0000000; 5740.323; 3189.068 ];

    [v2, theta, theta1, copa, errorstr] = gibbs(r1, r2, r3);

    fprintf(fid,'testgibbs %15.11f  %15.11f  %15.11f\n', v2(1), v2(2), v2(3));
    fprintf(fid,'testgibbs %15.11f  %15.11f  %15.11f\n', theta * rad, theta1 * rad, copa * rad);
end


function testhgibbs(fid)
    constastro;
    rad = 180.0 / pi;

    r1 = [ 0.0000000; 0.000000; re ];
    r2 = [ 0.0000000; -4464.696; -5102.509 ];
    r3 = [ 0.0000000; 5740.323; 3189.068 ];
    jd1 = 2451849.5;
    jd2 = jd1 + 1.0 / 1440.0 + 16.48 / 86400.0;
    jd3 = jd1 + 2.0 / 1440.0 + 33.04 / 86400.0;
    [v2, theta, theta1, copa, errorstr] = hgibbs(r1, r2, r3, jd1, jd2, jd3);

    fprintf(fid,'testherrgibbs %15.11f  %15.11f  %15.11f\n', v2(1), v2(2), v2(3));
    fprintf(fid,'testherrgibbs %15.11f  %15.11f  %15.11f\n', (theta * rad), (theta1 * rad), (copa * rad));
end



function testgeo(fid)
    constastro;
    rad = 180.0 / pi;

    % misc test
    dt = 86400.0;  % 1 day in second
    c22 = 1.57461532572292E-06;
    s22 = -9.03872789196567E-07;
    j22 = sqrt(c22 * c22 + s22 * s22);
    % stable longitude point
    lons = 0.5 * atan2(s22, c22);
    omegaearth = 0.000072921158553;   % rad /s
    z = 6.6017;  % rad
    % initial longitude with 0 initial drift
    lonp = 12.0 / rad;
    lonp = lons;

    lona = lonp - 1.0 / rad;
    lona = lons;
    londot = 0.0;
    for jj=0:400

        lona = lona + londot * dt;
        londot = 3.0 * omegaearth / z * sqrt(2.0 * j22) * sqrt(cos(2.0 * (lona - lons)) - cos(2.0 * (lonp - lons)));
        fprintf(fid,'%11.7f %11,7f %11.7f %11.7f\n', jj, lona * rad, (lonp - lons) * rad, londot * rad / 86400.0);
    end % for through all the tracks testing rtasc/decl rates
end


function doangles(jd, jdf, latgd, lon, alt, trtasc, tdecl, initguess, ansrlongstr, fida, fidas)
    constastro;
    eqeterms = 2;

    fileLoc = 'D:\Codes\LIBRARY\DataLib\';
    [iau80arr] = iau80in(fileLoc);

    % read existing data - this does not find x, y, s!
    eopFileName = 'D:\Codes\LIBRARY\DataLib\EOP-All-v1.1_2023-01-01.txt';
    [eoparr] = readeop(eopFileName);
  
    asecef1 = [0.0; 0.0; 0.0];
    asecef2 = [0.0; 0.0; 0.0];
    asecef3 = [0.0; 0.0; 0.0];
    
    %jd1 = jd(1) + jdf(1);
    %jd2 = jd(2) + jdf(2);
    %jd3 = jd(3) + jdf(3);

    [rsecef1, vsecef1] = site(latgd(1), lon(1), alt(1));
    [dut1, dat, lod, xp, yp, ddpsi, ddeps, ddx, ddy] = findeopparam(jd(1), jdf(1), 's', eoparr);
    %convtime(year(1), mon(1), day(1), hr(1), minute(1), second(1), 0, dut1, dat,
    %    out ut1, out tut1, out jdut1, out jdut1frac, out utc, out tai,
    %    out tt, out ttt, out jdtt, out jdttfrac, out tdb, out ttdb, out jdtdb, out jdtdbfrac);
    jdtt = jd(1);
    jdftt = jdf(1) + (dat + 32.184) / 86400.0;
    ttt = (jdtt + jdftt - 2451545.0) / 36525.0;
    % note you have to use tdb for time of interst AND j2000 (when dat = 32)
    %  ttt = (jd + jdFrac + (dat + 32.184) / 86400.0 - 2451545.0 - (32 + 32.184) / 86400.0) / 36525.0;
    jdut1 = jd(1) + jdf(1) + dut1 / 86400.0;
    [rseci1, vseci1, aseci1] = ecef2eci(rsecef1, vsecef1, asecef1, iau80arr, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps );

    [rsecef2, vsecef2] = site(latgd(2), lon(2), alt(2));
    [dut1, dat, lod, xp, yp, ddpsi, ddeps, ddx, ddy] = findeopparam(jd(2), jdf(2), 's', eoparr);
    %convtime(year(2), mon(2), day(2), hr(2), minute(2), second(2), 0, dut1, dat,
    %    out ut1, out tut1, out jdut1, out jdut1frac, out utc, out tai,
    %    out tt, out ttt, out jdtt, out jdttfrac, out tdb, out ttdb, out jdtdb, out jdtdbfrac);
    jdtt = jd(2);
    jdftt = jdf(2) + (dat + 32.184) / 86400.0;
    ttt = (jdtt + jdftt - 2451545.0) / 36525.0;
    jdut1 = jd(2) + jdf(2) + dut1 / 86400.0;
    [rseci2, vseci2, aseci2] = ecef2eci(rsecef2, vsecef2, asecef2, iau80arr, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps );
    [lst, gst] = lstime(lon(2), jdut1);
    fprintf(fida,'\nlst  %11.7f  %11.7f\n', lst, (lst * rad));


    [rsecef3, vsecef3] = site(latgd(3), lon(3), alt(3));
    [dut1, dat, lod, xp, yp, ddpsi, ddeps, ddx, ddy] = findeopparam(jd(3), jdf(3), 's', eoparr);
    jdtt = jd(3);
    jdftt = jdf(3) + (dat + 32.184) / 86400.0;
    ttt = (jdtt + jdftt - 2451545.0) / 36525.0;
    jdut1 = jd(3) + jdf(3) + dut1 / 86400.0;
    [rseci3, vseci3, aseci3] = ecef2eci(rsecef3, vsecef3, asecef3, iau80arr, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps );

    if (abs(latgd(1) - latgd(2)) < 0.001 && abs(latgd(1) - latgd(3)) < 0.001 ...
            && abs(lon(1) - lon(2)) < 0.001 && abs(lon(1) - lon(3)) < 0.001)
        diffsites = 'n';
    else
        diffsites = 'y';
    end

    % write output
    fprintf(fida,'rseci1  %15.11f  %15.11f  %15.11f\n', rseci1);
    fprintf(fida,'rseci2  %15.11f  %15.11f  %15.11f\n', rseci2);
    fprintf(fida,'rseci3  %15.11f  %15.11f  %15.11f\n', rseci3);

    los1(1) = cos(tdecl(1)) * cos(trtasc(1));
    los1(2) = cos(tdecl(1)) * sin(trtasc(1));
    los1(3) = sin(tdecl(1));

    los2(1) = cos(tdecl(2)) * cos(trtasc(2));
    los2(2) = cos(tdecl(2)) * sin(trtasc(2));
    los2(3) = sin(tdecl(2));

    los3(1) = cos(tdecl(3)) * cos(trtasc(3));
    los3(2) = cos(tdecl(3)) * sin(trtasc(3));
    los3(3) = sin(tdecl(3));

    fprintf(fida,'los1  %15.11f  %15.11f  %15.11f  %15.11f\n', los1, mag(los1));
    fprintf(fida,'los2  %15.11f  %15.11f  %15.11f  %15.11f\n', los2, mag(los2));
    fprintf(fida,'los3  %15.11f  %15.11f  %15.11f  %15.11f\n', los3, mag(los3));

    % to get initial guess, take measurements (1/2 and 2/3), assume circular orbit
    % find velocity and compare - just distinguish between LEO, GPS and GEO for now
    dt1 = (jd(2) - jd(1)) * 86400.0 + (jdf(2) - jdf(1)) * 86400.0;
    dt2 = (jd(3) - jd(2)) * 86400.0 + (jdf(3) - jdf(2)) * 86400.0;
    dtrtasc1 = (trtasc(2) - trtasc(1)) / dt1;
    dtrtasc2 = (trtasc(3) - trtasc(2)) / dt2;
    dtdecl1 = (tdecl(2) - tdecl(1)) / dt1;
    dtdecl2 = (tdecl(3) - tdecl(2)) / dt2;

    fprintf(fida,'rtasc  %15.11f  %15.11f  %15.11f\n', trtasc(1)*rad, trtasc(3)*rad, trtasc(3)*rad);
    fprintf(fida,'decl  %15.11f  %15.11f  %15.11f\n', tdecl(1)*rad, tdecl(3)*rad, tdecl(3)*rad);

    fprintf(fida,'');
    fprintf(fida,'Laplace -----------------------------------');
    fprintf(fidas,'Laplace -----------------------------------');

    fprintf(fida,'\n\ninputs:\n');
    fprintf(fida,'Site obs1  %15.11f  %15.11f  %15.11f km lat %15.11f lon %15.11f  %15.11f\n', rseci1(1), rseci1(2), rseci1(3), latgd(1)*rad, lon(1)*rad, alt(1));
    fprintf(fida,'Site obs2  %15.11f  %15.11f  %15.11f km lat %15.11f lon %15.11f  %15.11f\n', rseci2(1), rseci2(2), rseci2(3), latgd(2)*rad, lon(2)*rad, alt(2));
    fprintf(fida,'Site obs3  %15.11f  %15.11f  %15.11f km lat %15.11f lon %15.11f  %15.11f\n', rseci3(1), rseci3(2), rseci3(3), latgd(3)*rad, lon(3)*rad, alt(3));

    [iyear1, imon1, iday1, ihr1, iminute1, isecond1] = invjday(jd(1), jdf(1));
    fprintf(fida,'obs#1  %d %d %d %d:%d:%f  %15.11f  %15.11f\n', iyear1, imon1, iday1, ihr1, iminute1, isecond1, trtasc(1)*rad, tdecl(1)*rad);
    [iyear2, imon2, iday2, ihr2, iminute2, isecond2] = invjday(jd(2), jdf(2));
    fprintf(fida,'obs#2  %d %d %d %d:%d:%f  %15.11f  %15.11f\n', iyear2, imon2, iday2, ihr2, iminute2, isecond2, trtasc(2)*rad, tdecl(2)*rad);
    [iyear3, imon3, iday3, ihr3, iminute3, isecond3] = invjday(jd(3), jdf(3));
    fprintf(fida,'obs#3  %d %d %d %d:%d:%f  %15.11f  %15.11f\n', iyear3, imon3, iday3, ihr3, iminute3, isecond3, trtasc(3)*rad, tdecl(3)*rad);
    %if (caseopt == 2)
    %    diffsites = 'y';
    %else
    %diffsites = 'n';

    [r2, v2] = anglesl(tdecl(1), tdecl(2), tdecl(3), trtasc(1), trtasc(2), trtasc(3), ...
         jd(1), jdf(1), jd(2), jdf(2), jd(3), jdf(3), diffsites, rseci1, rseci2, rseci3, fida);   
    
    %fprintf(fida,errstr);
    fprintf(fida,'r2 %15.11f  %15.11f  %15.11f v2 %15.11f  %15.11f  %15.11f\n', r2(1), r2(2), r2(3), v2(1), v2(2), v2(3));
    [p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper ] = rv2coe (r2, v2);
    fprintf(fida,'\nlaplace coes a= %15.11f e = %15.11f  i = %15.11f  %15.11f  %15.11f  %15.11f  %15.11f  %15.11f\n', ...
          a, ecc, incl*rad, raan*rad, argp*rad, nu*rad, m*rad, arglat*rad); %
    fprintf(fida,'%s\n', ansrlongstr);
    fprintf(fidas,'\nlaplace coes a= %15.11f e = %15.11f  i = %15.11f  %15.11f  %15.11f  %15.11f  %15.11f  %15.11f\n', ...
          a, ecc, incl*rad, raan*rad, argp*rad, nu*rad, m*rad, arglat*rad); %
    %strbuildallsum.AppendLine(ansr);

    fprintf(fida,'Gauss  -----------------------------------');
    fprintf(fidas,'Gauss  -----------------------------------');
    % if (caseopt == 23)
    %     % curtis example -many mistakes!
    %     rseci1 = [ 3489.8; 3430.2; 4078.5 ];
    %     rseci2 = [ 3460.1; 3460.1; 4078.5 ];
    %     rseci3 = [ 3429.9; 3490.1; 4078.5 ];
    % end
    [r2, v2] = anglesg(tdecl(1), tdecl(2), tdecl(3), trtasc(1), trtasc(2), trtasc(3), ...
         jd(1), jdf(1), jd(2), jdf(2), jd(3), jdf(3), diffsites, rseci1, rseci2, rseci3, fida);   
    %fprintf(fida,errstr);
    fprintf(fida,'r2 %15.11f  %15.11f  %15.11f v2 %15.11f  %15.11f  %15.11f\n', r2(1), r2(2), r2(3), v2(1), v2(2), v2(3));
    [p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper ] = rv2coe (r2, v2);
    fprintf(fida,'\nguass coes a= %15.11f e = %15.11f  i = %15.11f  %15.11f  %15.11f  %15.11f  %15.11f  %15.11f\n', ...
          a, ecc, incl*rad, raan*rad, argp*rad, nu*rad, m*rad, arglat*rad); %
    fprintf(fida,'%s\n', ansrlongstr);
    fprintf(fidas,'\nguass coes a= %15.11f e = %15.11f  i = %15.11f  %15.11f  %15.11f  %15.11f  %15.11f  %15.11f\n', ...
          a, ecc, incl*rad, raan*rad, argp*rad, nu*rad, m*rad, arglat*rad); %
    %strbuildallsum.AppendLine(ansr);

     pctchg = 0.05;
    fprintf(fida,'double-r -----------------------------------' );
    fprintf(fidas,'double-r -----------------------------------' );
    % initial guesses needed for -r and Gooding
    % use result from Gauss as it's usually pretty good
    % this seems to really help Gooding!!
    [bigr2x] = getGaussRoot(tdecl(1), tdecl(2), tdecl(3), trtasc(1), trtasc(2), trtasc(3), ...
              jd(1), jdf(1), jd(2), jdf(2), jd(3), jdf(3), rseci1, rseci2, rseci3);
    initguess = bigr2x;

    rng1 = initguess;  % old 12500 needs to be in km!! seems to do better when all the same? if too far off (*2) NAN
    rng2 = initguess * 1.02;  % 1.02 might be better? make the initial guess a bit different
    rng3 = initguess * 1.08;
    [r2, v2] = anglesdr(tdecl(1), tdecl(2), tdecl(3), trtasc(1), trtasc(2), trtasc(3), ...
         jd(1), jdf(1), jd(2), jdf(2), jd(3), jdf(3), diffsites, rseci1, rseci2, rseci3, rng1, rng2, pctchg);   

    %fprintf(fida,errstr);
    fprintf(fida,'r2 %15.11f  %15.11f  %15.11f v2 %15.11f  %15.11f  %15.11f\n', r2(1), r2(2), r2(3), v2(1), v2(2), v2(3));
    [p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper ] = rv2coe (r2, v2);
    fprintf(fida,'\ndouble r coes a= %15.11f e = %15.11f  i = %15.11f  %15.11f  %15.11f  %15.11f  %15.11f  %15.11f\n', ...
          a, ecc, incl*rad, raan*rad, argp*rad, nu*rad, m*rad, arglat*rad); %
    fprintf(fida,'%s\n', ansrlongstr);
    fprintf(fidas,'\ndouble r coes a= %15.11f e = %15.11f  i = %15.11f  %15.11f  %15.11f  %15.11f  %15.11f  %15.11f\n', ...
          a, ecc, incl*rad, raan*rad, argp*rad, nu*rad, m*rad, arglat*rad); %
    %strbuildallsum.AppendLine(ansr);


    fprintf(fida,'Gooding -----------------------------------\n');
    fprintf(fidas,'Gooding -----------------------------------\n');
    numhalfrev = 0;

    % [bigr2x] = getGaussRoot(tdecl(1), tdecl(2), tdecl(3), trtasc(1), trtasc(2), trtasc(3), ...
    %           jd(1), jdf(1), jd(2), jdf(2), jd(3), jdf(3), rseci1, rseci2, rseci3);
    % initguess = bigr2x;
    % 
    % rng1 = initguess;  % old 12500 needs to be in km!! seems to do better when all the same? if too far off (*2) NAN
    % rng2 = initguess * 1.02;  % 1.02 might be better? make the initial guess a bit different
    % rng3 = initguess * 1.08;
    % 
    % [r2, v2] = anglesgood(tdecl(1), tdecl(2), tdecl(3), trtasc(1), trtasc(2), trtasc(3), ...
    %      jd(1), jdf(1), jd(2), jdf(2), jd(3), jdf(3), diffsites, rseci1, rseci2, rseci3, rng1, rng2, pctchg, fid);   

    %fprintf(fida,errstr);
    % fprintf(fida,'r2 %15.11f  %15.11f  %15.11f v2 %15.11f  %15.11f  %15.11f\n', r2(1), r2(2), r2(3), v2(1), v2(2), v2(3));
    % [p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper ] = rv2coe (r2, v2);
    % fprintf(fida,'\ngooding coes a= %15.11f e = %15.11f  i = %15.11f  %15.11f  %15.11f  %15.11f  %15.11f  %15.11f\n', ...
    %       a, ecc, incl*rad, raan*rad, argp*rad, nu*rad, m*rad, arglat*rad); %
    % %fprintf(fida,ansr);
    % fprintf(fid,'\ngooding coes a= %15.11f e = %15.11f  i = %15.11f  %15.11f  %15.11f  %15.11f  %15.11f  %15.11f\n', ...
    %       a, ecc, incl*rad, raan*rad, argp*rad, nu*rad, m*rad, arglat*rad); %
    %strbuildallsum.AppendLine(ansr);

    % directory = 'D:\Codes\LIBRARY\cs\TestAll\';
    % fprintf(fid,'angles only tests case results written to ' + directory + 'testall-Angles.out ');
    % fprintf(fid,'geo data for chap 9 plot written to D:\faabook\current\excel\testgeo.out for ch9 plot ');

    % File.WriteAllText(directory + 'testall-Angles.out', strbuildall);
    % File.WriteAllText(directory + 'testall-Anglessum.out', strbuildallsum);
end  % doangles


% test angles-only routines
% output these results separately to the testall directory
function testangles(fid)
    %conv = pi / (180.0 * 3600.0);
    rad = 180.0 / pi;
    errstr = '';
    diffsites = 'n';
    %StringBuilder strbuildall = new StringBuilder(fid);
    %StringBuilder strbuildallsum = new StringBuilder(fid);

    directory = 'D:\Codes\LIBRARY\matlab\';
    fida = fopen(strcat(directory,'testall-Anglesm.out'), 'wt');
    fidas = fopen(strcat(directory,'testall-Anglessumm.out'), 'wt');

    fileLoc = 'D:\Codes\LIBRARY\DataLib\';
    [iau80arr] = iau80in(fileLoc);

    eopFileName = 'D:\Codes\LIBRARY\DataLib\EOP-All-v1.1_2023-01-01.txt';
    [eoparr] = readeop(eopFileName);

    % gooding tests cases from Gooding paper (1997 CMDA)

    % read input data
    % note the input data has a # line between each case
    % special case needed so ony 3 obs for each case...
    infilename = 'D:\Codes\LIBRARY\DataLib\';
    infile = fopen(append(infilename, 'anglestestmat.dat'), 'r');

    % --- read obs data in
    %caseopt = 0;  % set this for whichever case to run

    % find mins
    % orbits only need to be close
    % rtol = 500.0; % km
    % vtol = 0.1;   % km/s
    % atol = 500.0; % km
    % ptol = 500.0; % km
    % etol = 0.1;  %
    itol = 5.0 / rad;   % rad

    strt = 1;
    stp = 32;
    initguess = 0.0;
    for caseopt = strt: stp
        fprintf(1,'caseopt %d\n', caseopt);
        fprintf(fida,'caseopt %d\n', caseopt);

        % get first heard line
        longstr = fgets(infile);

        while (~feof(infile))

            longstr = fgets(infile);
            ansrlongstr = longstr;

            fprintf(fida,'\n\n ================================ case number %d ================================\n', caseopt);
            fprintf(fid,'\n\n ================================ case number %d ================================\n', caseopt);
            fprintf(fid,'%s \n', longstr);

            % # 0 ansr a 12246.023  e 0.2000  i 40.00  W 330.000  w 0.0  nu 0.0
            % 20,  8,  2012,  11,  40,  28.00,    40.000,    -110.000,     2.0000,     0.939913  ,    18.667717, x
            % 20,  8,  2012,  11,  48,  28.00,    40.000,    -110.000,     2.0000,     45.025748 ,    35.664741, x
            % 20,  8,  2012,  11,  52,  28.00,    40.000,    -110.000,     2.0000,     67.886655 ,    36.996583, x
            if ~contains(longstr, '#')
                for i=1:3
                    obsktr = i;
                    linesplt = split(longstr, ',');  % returns strings
                    day(obsktr) = str2num(cell2mat(linesplt(1)));
                    mon(obsktr) = str2num(cell2mat(linesplt(2)));
                    year(obsktr) = str2num(cell2mat(linesplt(3)));
                    hr(obsktr) = str2num(cell2mat(linesplt(4)));
                    minute(obsktr) = str2num(cell2mat(linesplt(5)));
                    second(obsktr) = str2num(cell2mat(linesplt(6)));
                    [jd(obsktr), jdf(obsktr)] = jday(year(obsktr), mon(obsktr), day(obsktr), hr(obsktr), minute(obsktr), second(obsktr) );

                    latgd(obsktr) = str2num(cell2mat(linesplt(7))) / rad;
                    lon(obsktr) = str2num(cell2mat(linesplt(8))) / rad;
                    alt(obsktr) = str2num(cell2mat(linesplt(9))) / rad;

                    trtasc(obsktr) = str2num(cell2mat(linesplt(10))) / rad;
                    tdecl(obsktr) = str2num(cell2mat(linesplt(11))) / rad;
                    % if (obsktr == 0)
                    %     initguess(tmpcase) = str2num(cell2mat(linesplt(12)));  % initial guess in km
                    % end

                    longstr = fgets(infile);
                    obsktr = obsktr + 1;
                end

                % write summary results to fid
                doangles(jd, jdf, latgd, lon, alt, trtasc, tdecl, initguess, ansrlongstr, fida, fid);

            end  % if

        end  % while not eof

    end  % for caseopt


end   % testangles





function testlambertumins(fid)
   constastro;
    r1 = [ 2.500000 * re; 0.000000; 0.000000 ];
    r2 = [ 1.9151111 * re; 1.6069690 * re; 0.000000 ];
    dm = 'S';
    %char de = 'L';
    nrev = 0;

    for i=0:10  % 1000

        for j=1:5
            [tmin, tminp, tminenergy] = lambertminT(r1, r2, dm, 'L', nrev);
        end

        for j=1:5
            [tmin, tminp, tminenergy] = lambertminT(r1, r2, dm, 'H', nrev);
        end

        [minenergyv, aminenergy, tminenergy, tminabs] = lambertmin ( r1, r2, 'd', 0 );
        fprintf(fid,' minenergyv %16.8f %16.8f %16.8f a %15.11f  dt %15.11f  %15.11f\n', minenergyv, aminenergy, tminenergy, tminabs );

        [minenergyv, aminenergy, tminenergy, tminabs] = lambertmin ( r1, r2, 'r', 0 );
        fprintf(fid,' minenergyv %16.8f %16.8f %16.8f a %15.11f  dt %15.11f  %15.11f\n', minenergyv, aminenergy, tminenergy, tminabs );
    end

    % timing of routines
    for i=0:10  %00

        [ tbi, tbil] = lambgettbiu(r1, r2, 5);
        fprintf(fid,' r1 %16.8f %16.8f %16.8f\n',r1 );
        fprintf(fid,' r2 %16.8f %16.8f %16.8f\n',r2 );
        fprintf(fid,'From universal variables \n%15.11f %15.11f s\n',tbi(1,1),tbi(1,2));
        fprintf(fid,'%15.11f %15.11f s\n',tbi(2,1),tbi(2,2));
        fprintf(fid,'%15.11f %15.11f s\n',tbi(3,1),tbi(3,2));
        fprintf(fid,'%15.11f %15.11f s\n',tbi(4,1),tbi(4,2));
        fprintf(fid,'%15.11f %15.11f s \n\n',tbi(5,1),tbi(5,2));
        fprintf(fid,'%15.11f %15.11f s\n',tbil(1,1),tbil(1,2));
        fprintf(fid,'%15.11f %15.11f s\n',tbil(2,1),tbil(2,2));
        fprintf(fid,'%15.11f %15.11f s\n',tbil(3,1),tbil(3,2));
        fprintf(fid,'%15.11f %15.11f s\n',tbil(4,1),tbil(4,2));
        fprintf(fid,'%15.11f %15.11f s\n',tbil(5,1),tbil(5,2));

        [kbi, tof] = lambertumins(r1, r2, 1, 'S');
        tbiSu(2, 2) = kbi;
        tbiSu(2, 3) = tof;
        [kbi, tof] = lambertumins(r1, r2, 2, 'S');
        tbiSu(3, 2) = kbi;
        tbiSu(3, 3) = tof;
        [kbi, tof] = lambertumins(r1, r2, 3, 'S');
        tbiSu(4, 2) = kbi;
        tbiSu(4, 3) = tof;
        [kbi, tof] = lambertumins(r1, r2, 4, 'S');
        tbiSu(5, 2) = kbi;
        tbiSu(5, 3) = tof;
        [kbi, tof] = lambertumins(r1, r2, 5, 'S');
        tbiSu(6, 2) = kbi;
        tbiSu(6, 3) = tof;

        [kbi, tof] = lambertumins(r1, r2, 1, 'L');
        tbiLu(2, 2) = kbi;
        tbiLu(2, 3) = tof;
        [kbi, tof] = lambertumins(r1, r2, 2, 'L');
        tbiLu(3, 2) = kbi;
        tbiLu(3, 3) = tof;
        [kbi, tof] = lambertumins(r1, r2, 3, 'L');
        tbiLu(4, 2) = kbi;
        tbiLu(4, 3) = tof;
        [kbi, tof] = lambertumins(r1, r2, 4, 'L');
        tbiLu(5, 2) = kbi;
        tbiLu(5, 3) = tof;
        [kbi, tof] = lambertumins(r1, r2, 5, 'L');
        tbiLu(6, 2) = kbi;
        tbiLu(6, 3) = tof;
    end
    
    tusec = 806.8111238242922;
    ootusec = 1.0 / tusec;
    
    [s, tau] = lambertkmins1st(r1, r2);
    
    % for general cases, use 'x' for dm to get the tof/kbi values
    for i=0:10 % 00
    
        [kbi, tof] = lambertkmins(s, tau, 1, 'x', 'L');
        tbidk(2, 2) = kbi;
        tbidk(2, 3) = tof * ootusec;
        [kbi, tof] = lambertkmins(s, tau, 2, 'x', 'L');
        tbidk(3, 2) = kbi;
        tbidk(3, 3) = tof * ootusec;
        [kbi, tof] = lambertkmins(s, tau, 3, 'x', 'L');
        tbidk(4, 2) = kbi;
        tbidk(4, 3) = tof * ootusec;
        [kbi, tof] = lambertkmins(s, tau, 4, 'x', 'L');
        tbidk(5, 2) = kbi;
        tbidk(5, 3) = tof * ootusec;
        [kbi, tof] = lambertkmins(s, tau, 5, 'x', 'L');
        tbidk(6, 2) = kbi;
        tbidk(6, 3) = tof * ootusec;
    
        [kbi, tof] = lambertkmins(s, tau, 1, 'x', 'H');
        tbirk(2, 2) = kbi;
        tbirk(2, 3) = tof * ootusec;
        [kbi, tof] = lambertkmins(s, tau, 2, 'x', 'H');
        tbirk(3, 2) = kbi;
        tbirk(3, 3) = tof * ootusec;
        [kbi, tof] = lambertkmins(s, tau, 3, 'x', 'H');
        tbirk(4, 2) = kbi;
        tbirk(4, 3) = tof * ootusec;
        [kbi, tof] = lambertkmins(s, tau, 4, 'x', 'H');
        tbirk(5, 2) = kbi;
        tbirk(5, 3) = tof * ootusec;
        [kbi, tof] = lambertkmins(s, tau, 5, 'x', 'H');
        tbirk(6, 2) = kbi;
        tbirk(6, 3) = tof * ootusec;
    end
    
    fprintf(fid,'time for Lambert kmin\n' );
    end
    
function testlambertminT(fid)
    constastro;
    r1 = [ 2.500000 * re; 0.000000; 0.000000 ];
    r2 = [ 1.9151111 * re; 1.6069690 * re; 0.000000 ];
    dm = 'S';
    de = 'L';
    nrev = 0;
    
    [tmin, tminp, tminenergy] = lambertminT(r1, r2, dm, de, nrev);
    fprintf(fid,'lambertmtmin  s %15.11f  minp  %15.11f minenr  %15.11f\n', tmin, tminp, tminenergy);
    
    [tmin, tminp, tminenergy] = lambertminT(r1, r2, dm, 'H', nrev);
    fprintf(fid,'lambertmtmin  s %15.11f  minp  %15.11f minenr  %15.11f\n', tmin, tminp, tminenergy);
end
    

function testlambhodograph(fid)
    constastro;
    rad = 180.0 / pi;

    r1 = [ 2.500000 * re; 0.000000; 0.000000 ];
    r2 = [ 1.9151111 * re; 1.6069690 * re; 0.000000 ];
    % assume circular initial orbit for vel calcs
    v1 = [ 0.0; sqrt(mu / r1(1)); 0.0 ];
    p = 12345.235;  % km
    ecc = 0.023487;
    dnu = 34.349128 / rad;
    dtsec = 92854.234;
    
    [v1t, v2t] = lambhodograph(r1, r2, v1, p, ecc, dnu, dtsec);
    
    fprintf(fid,'lamb hod %15.11f  %15.11f  %15.11f \nlamb hod %15.11f  %15.11f  %15.11f\n', ...
         v1t(1), v1t(2), v1t(3), v2t(1), v2t(2), v2t(3));
end
    

function testlambertbattin(fid)
    constastro;

    r1 = [ 2.500000 * re; 0.000000; 0.000000 ];
    r2 = [ 1.9151111 * re; 1.6069690 * re; 0.000000 ];
    v1 = [ 0.0; sqrt(mu / r1(1)); 0.0 ];
    dm = 'S';
    de = 'L';
    nrev = 0;
    dtsec = 76.0 * 60.0;
    altpadc = 100.0 / re;  %er
    dtwait = 0.0;
    
    [v1t, v2t] = lambertb(r1, r2, v1, dm, de, nrev, dtsec);
    
    fprintf(fid,'lambertbattin %15.11f  %15.11f  %15.11f \nllambertbattin %15.11f  %15.11f  %15.11f\n', ...
         v1t(1), v1t(2), v1t(3), v2t(1), v2t(2), v2t(3));
end
    

function testeq2rv(fid)
    a = 7236.346;
    af = 0.23457;
    ag = 0.47285;
    chi = 0.23475;
    psi = 0.28374;
    meanlon = 2.230482378;
    fr = 1;

    [r, v] = eq2rv(a, af, ag, chi, psi, meanlon, fr);

    fprintf(fid,'eq2rv %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', r(1), r(2), r(3), v(1), v(2), v(3));
end


function testrv2eq(fid)
    constastro;
    r = [ 2.500000 * re; 0.000000; 0.000000 ];
    % assume circular initial orbit for vel calcs
    v = [ 0.0; sqrt(mu / r(1)); 0.0 ];

    [ a, n, af, ag, chi, psi, meanlonM, meanlonNu, fr ] = rv2eq (r, v);

    fprintf(fid,'rv2eq   a km         n rad      af           ag         chi          psi      meanlonnu deg   meanlonm deg\n');
    fprintf(fid,'rv2eq   a %11.4f n  %11.4f af %13.9f ag %13.7f chi %11.5f psi %11.5f meanlonNu %11.5f meanlonM %11.5f\n',...
        a, n, af, ag, chi, psi, meanlonNu*rad, meanlonM*rad );
end


function testAllLamb(fid)
    infilename = 'D:\Codes\LIBRARY\DataLib\';
    infile = fopen(append(infilename, 'lamberttest.dat'), 'r');

    methodType = 'lambertu';

    strt = 1;
    stp = 82;

    cc = "x";
    runopt = "all";
    dmin = 'x';
    dein = 'x';
    obsktr = 0;

    directory = 'D:\Codes\LIBRARY\matlab\';
    fida = fopen(strcat(directory,'testall-Lambertm.out'), 'wt');
    %fidas = fopen(strcat(directory,'testall-Lambertsumm.out'), 'wt');


    % run all the tests
    for caseopt = strt: stp
        fprintf(fid,'caseopt %d\n', caseopt);

        while (~feof(infile))

            longstr = fgets(infile);
            %    1
            % 4938.49830042171, -1922.24810472241, 4384.68293292613
            % 0.738204644165659, 7.20989453238397, 2.32877392066299
            % -1105.78023519582, 2373.16130661458, 6713.89444816503
            % 5.4720951867079, -4.39299050886976, 2.45681739563752
            % 6372.69272563561

            if ~contains(longstr, '#')
                nrev = longstr(1);
                longstr = fgets(infile);
                linesplt = split(longstr, ',');  % returns strings
                r1(1) = str2num(cell2mat(linesplt(1)));
                r1(2) = str2num(cell2mat(linesplt(2)));
                r1(3) = str2num(cell2mat(linesplt(3)));
                longstr = fgets(infile);
                v1(1) = str2num(cell2mat(linesplt(1)));
                v1(2) = str2num(cell2mat(linesplt(2)));
                v1(3) = str2num(cell2mat(linesplt(3)));
                longstr = fgets(infile);
                r2(1) = str2num(cell2mat(linesplt(1)));
                r2(2) = str2num(cell2mat(linesplt(2)));
                r2(3) = str2num(cell2mat(linesplt(3)));
                longstr = fgets(infile);
                v2(1) = str2num(cell2mat(linesplt(1)));
                v2(2) = str2num(cell2mat(linesplt(2)));
                v2(3) = str2num(cell2mat(linesplt(3)));
                longstr = fgets(infile);
                dtsec = str2num(cell2mat(linesplt(1)));

                dtwait = 0.0;
                
                obsktr = obsktr + 1;

                fprintf(fid,'\n\n ================================ case number %d ================================\n', caseopt);

                % , ref strbuildall, ref strbuildsum, out caseerr
                methodType = 'lambertu';
                runopt = 'all';
                % write summary results to fid
                dolamberttests(r1, v1, r2, v2, dtwait, dtsec, nrev, dmin, dein, runopt, methodType, fida, fid);

            end  % if

        end  % while not eof

    end  % for caseopt

end   % testAllLamb


% test building the lambert envelope
function dolamberttests(r1, v1, r2, v2, dtwait, dtsec, nrev, dmin, dein, runopt, methodType, fida, fidas);
    constastro;
    detailSum = "";
    detailAll = "";
    caseerr = "";

    show = 'y';  % display results to form

    % implementation
    altpadc = 200.0 / re;  % set 200 km for altitude you set as the over limit.
    tusec = 806.8111238242922;
    numiter = 16;

    v1 = [0.0 0.0 0.0];
    v2 = [0.0 0.0 0.0];

    magr1 = mag(r1);
    magr2 = mag(r2);

    % this value stays constant in all calcs, vara changes with df
    cosdeltanu = dot(r1, r2) / (magr1 * magr2);

    %fprintf(fid,'now do findtbi calcs\n');
    %fprintf(fid,'iter       y         dtnew          psiold      psinew   psinew-psiold   dtdpsi      dtdpsi2    lower    upper\n');

    %     [s, tau] = lambertkmins1st(r1, r2);
    % fprintf(fid,' s " + s.ToString(fmt) + " tau " + tau.ToString(fmt));

    % [kbi, tof] = lambertkmins(s, tau, nrev, 'x', 'L');
    % tbidk[1, 1] = kbi;
    % tbidk[1, 2] = tof / tusec;
    % 
    % [,] tbirk = new [10, 3];
    % [kbi, tof] = lambertkmins(s, tau, nrev, 'x', 'H');
    % tbirk[1, 1] = kbi;
    % tbirk[1, 2] = tof / tusec;
    % 
    % fprintf(fid,'From k variables\n');
    % fprintf(fid,' " + tbidk(2,2), "  " + (tbidk[1, 2] * tusec).ToString("0.00000000000") + " s " + (tbidk[1, 2]).ToString("0.00000000000") + " tu\n');
    % fprintf(fid,'\n');
    % fprintf(fid,' " + tbirk(2,2), + "  " + (tbirk[1, 2] * tusec).ToString("0.00000000000") + " s " + (tbirk[1, 2]).ToString("0.00000000000") + " tu\n');


    % 1 normal prints
    % 3 sum prints
    caseopt = 1; % temp only

    fprintf(fida,'lambertTest  %d %15.11f  %15.11f  %15.11f  \n %15.11f  %15.11f  %15.11f\n %15.11f  %15.11f  %15.11f\n %15.11f  %15.11f  %15.11f  \n %15.11f  %15.11f  %15.11f %15.11f %15.11f\n', ...
        caseopt, r1(1),  r1(2),  r1(3),  v1(1),  v1(2),  v1(3),  ...
        r2(1),  r2(2),  r2(3),  v2(1),  v2(2),  v2(3),  ...
        nrev, dmin, dein, dtwait, dtsec);

    %AstroLibr.lambertminT(r1, r2, 'S', 'L', 1, out tmin, out tminp, out tminenergy);
    %fprintf(fid,'mint S " + tmin.ToString("0.0000") + " minp " + tminp.ToString("0.0000") + " minener " + tminenergy.ToString("0.0000"));
    %AstroLibr.lambertminT(r1, r2, 'S', 'L', 2, out tmin, out tminp, out tminenergy);
    %fprintf(fid,'mint S " + tmin.ToString("0.0000") + " minp " + tminp.ToString("0.0000") + " minener " + tminenergy.ToString("0.0000"));
    %AstroLibr.lambertminT(r1, r2, 'S', 'L', 3, out tmin, out tminp, out tminenergy);
    %fprintf(fid,'mint S " + tmin.ToString("0.0000") + " minp " + tminp.ToString("0.0000") + " minener(" + tminenergy.ToString("0.0000"));

    %AstroLibr.lambertminT(r1, r2, 'L', 'H', 1, out tmin, out tminp, out tminenergy);
    %fprintf(fid,'mint L " + tmin.ToString("0.0000") + " minp " + tminp.ToString("0.0000") + " minener " + tminenergy.ToString("0.0000"));
    %AstroLibr.lambertminT(r1, r2, 'L', 'H', 2, out tmin, out tminp, out tminenergy);
    %fprintf(fid,'mint L " + tmin.ToString("0.0000") + " minp " + tminp.ToString("0.0000") + " minener " + tminenergy.ToString("0.0000"));
    %AstroLibr.lambertminT(r1, r2, 'L', 'H', 3, out tmin, out tminp, out tminenergy);
    %fprintf(fid,'mint L " + tmin.ToString("0.0000") + " minp " + tminp.ToString("0.0000") + " minener " + tminenergy.ToString("0.0000"));

    modecon = 'n';  % 'c' to shortcut bad cases (hitearth) at iter 3 or 'n'
    rad = 180.0 / pi;

    fprintf(fida,' TEST ------------------ s/l  L  0 rev ------------------\n');
    hitearth = ' ';
    dm = 'S';
    de = 'L';
    if (methodType == "lambertk" || (runopt == "all" || (dein == de && nrev == 0)))
        tbi(1,1) = 0.0;
        tbi(1,2) = 0.0;
        [v1tk, v2tk, errorsum, errorout] = lambertk ( r1, r2, v1, dm, de, nrev, 0.0, dtsec, tbi, show, 'n' )
        fprintf(fida,' %s\n',detailAll);
        fprintf(fidas,' %s\n', errorout);
        %fprintf(fida,'k#" + caseopt  detailSum + " diffs " + MathTimeLib::mag(dr).ToString("0.00000000000"));
        fprintf(fida,'lamk v1t %15.11f  %15.11f  %15.11f\n', v1tk(1),  v1tk(2),  v1tk(3));
        fprintf(fidas,'lamk v2t %15.11f  %15.11f  %15.11f\n', v2tk(1),  v2tk(2),  v2tk(3));
        %fprintf(fid,'magv1t.ToString("0.0000000").PadLeft(12)  magv2t.ToString("0.0000000").PadLeft(12));

        [r3h, v3h] = kepler(r1, v1tk, dtsec);
        for j=1:3
            dr(j) = r2(j) - r3h(j);
        end
        % switch for debugs
        fprintf(fida,'r3h %15.11f  %15.11f  %15.11f dr %15.11f\n', r3h(1),  r3h(2),  r3h(3), mag(dr));
        if (mag(dr) > 0.05)
            fprintf(fida,'dm, de nnrev %3d velb does not get to r2 position %15.11f km \n\n', nrev, mag(dr));
            fprintf(fidas,'dm, de nnrev %3d velb does not get to r2 position %15.11f km \n\n', nrev, mag(dr));
        end
        [p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper ] = rv2coe (r1, v1tk);

        tmpstr = sprintf(' %s  %15.11f %15.11f  %15.11f\n', hitearth, mag(dr), a, ecc);
        casearr = strrep(errorout,'\s+', tmpstr)
    end

    if (methodType == "lambertu" || (runopt == "all" || (dein == de && nrev == 0)))
        [kbi, tof] = lambertumins( r1, r2, 0, dm ) ;   
        [v1tu, v2tu, errorl] = lambertu(r1, r2, v1, dm, de, 0, dtsec, kbi, fida );
        %fprintf(fid,'detailSum);
        fprintf(fida,'univ v1t %15.11f  %15.11f  %15.11f\n', v1tu(1),  v1tu(2),  v1tu(3));
        fprintf(fidas,'univ v2t %15.11f  %15.11f  %15.11f\n', v2tu(1),  v2tu(2),  v2tu(3));
        [r3h, v3h] = kepler(r1, v1tu, dtsec);
        for j=1:3
            dr(j) = r2(j) - r3h(j);
        end
        if (mag(dr) > 0.05)
            fprintf(fida,'dm, de nnrev %3d velb does not get to r2 position %15.11f km \n\n', nrev, mag(dr));
            fprintf(fidas,'dm, de nnrev %3d velb does not get to r2 position %15.11f km \n\n', nrev, mag(dr));
        end

        for j=1:3
            dv1(j) = v1tk(j) - v1tu(j);
            dv2(j) = v2tk(j) - v2tu(j);
        end
        if (mag(dv1) > 0.01 || mag(dv2) > 0.01)
            fprintf(fida,'velk does not match velu \n\n');

            [p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper ] = rv2coe (r1, v1tu);
        end
    end

    if (methodType == "lambertb" || (runopt == "all" || (dein == de && nrev == 0)))
        [v1tb, v2tb, hitearth, errorsum, errorout] = lambertb(r1, r2, v1, dm, de, 0, 0.0, dtsec, altpadc, 'y');
        %fprintf(fid,'detailSum);
        fprintf(fida,'batt v1t %15.11f  %15.11f  %15.11f\n', v1tb(1),  v1tb(2),  v1tb(3));
        fprintf(fidas,'batt v2t %15.11f  %15.11f  %15.11f\n', v2tb(1),  v2tb(2),  v2tb(3));
        [r3h, v3h] = kepler(r1, v1tb, dtsec);
        for j=1:3
            dr(j) = r2(j) - r3h(j);
        end
        if (mag(dr) > 0.05)
            fprintf(fida,'dm, de nnrev %3d velb does not get to r2 position %15.11f km \n\n', nrev, mag(dr));
            fprintf(fidas,'dm, de nnrev %3d velb does not get to r2 position %15.11f km \n\n', nrev, mag(dr));
        end
        %fprintf(fid,'diffs " + mag(dr).ToString("0.00000000000"));

        for j=1:3
            dv1(j) = v1tk(j) - v1tb(j);
            dv2(j) = v2tk(j) - v2tb(j);
        end
        if (mag(dv1) > 0.01 || mag(dv2) > 0.01)
            fprintf(fida,'velk does not match velb \n\n');
        end
        %fprintf(fid,'diffs " + MathTimeLib::mag(dr).ToString("0.00000000000"));

        [p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper ] = rv2coe (r1, v1tb);
    end

    fprintf(fida,' TEST ------------------ s/l H 0 rev ------------------\n');
    dm = 'L';
    de = 'H';
    if (methodType == "lambertk" || (runopt == "all" || (dein == de && nrev == 0)))
        % k near 180 is about 53017 while battin is 30324!
        [v1tk, v2tk, errorsum, errorout] = lambertk ( r1, r2, v1, dm, de, 0, 0.0, dtsec, 0.0, show, 'n' )
        fprintf(fida,' %s\n',detailSum);
        fprintf(fidas,' %s\n', errorout);
        %fprintf(fid,'k#" + caseopt  detailSum + " diffs " + mag(dr).ToString("0.00000000000"));
        fprintf(fida,'lamk v1t %15.11f  %15.11f  %15.11f\n', v1tk(1),  v1tk(2),  v1tk(3));
        fprintf(fidas,'lamk v2t %15.11f  %15.11f  %15.11f\n', v2tk(1),  v2tk(2),  v2tk(3));
        %fprintf(fid,'magv1t.ToString("0.0000000").PadLeft(12)  magv2t.ToString("0.0000000").PadLeft(12));

        [r3h, v3h] = kepler(r1, v1tk, dtsec);
        for j=1:3
            dr(j) = r2(j) - r3h(j);
        end
        % switch for debugs
        fprintf(fida,'r3h %15.11f  %15.11f  %15.11f dr %15.11f\n', r3h(1),  r3h(2),  r3h(3), mag(dr));
        if (mag(dr) > 0.05)
            fprintf(fida,'dm, de nnrev %3d velb does not get to r2 position %15.11f km \n\n', nrev, mag(dr));
            fprintf(fidas,'dm, de nnrev %3d velb does not get to r2 position %15.11f km \n\n', nrev, mag(dr));
        end

        [p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper ] = rv2coe (r1, v1tk);
        tmpstr = sprintf(' %s  %15.11f %15.11f  %15.11f\n', hitearth, mag(dr), a, ecc);
        casearr = strrep(errorout,'\s+', tmpstr)
    end

    if (methodType == "lambertu" || (runopt == "all" || (dein == de && nrev == 0)))
        [kbi, tof] = lambertumins( r1, r2, 0, dm ) ;   
        [v1tu, v2tu, errorl] = lambertu(r1, r2, v1, dm, de, 0, dtsec, kbi, fida );
        %fprintf(fid,'detailSum);
        fprintf(fida,'univ v1t %15.11f  %15.11f  %15.11f\n', v1tu(1),  v1tu(2),  v1tu(3));
        fprintf(fidas,'univ v2t %15.11f  %15.11f  %15.11f\n', v2tu(1),  v2tu(2),  v2tu(3));
        [r3h, v3h] = kepler(r1, v1tu, dtsec);
        for j=1:3
            dr(j) = r2(j) - r3h(j);
        end
        if (mag(dr) > 0.05)
            fprintf(fida,'dm, de nnrev %3d velb does not get to r2 position %15.11f km \n\n', nrev, mag(dr));
            fprintf(fidas,'dm, de nnrev %3d velb does not get to r2 position %15.11f km \n\n', nrev, mag(dr));
        end

        for j=1:3
            dv1(j) = v1tk(j) - v1tu(j);
            dv2(j) = v2tk(j) - v2tu(j);
        end
        if (mag(dv1) > 0.01 || mag(dv2) > 0.01)
            fprintf(fida,'velk does not match velu \n\n');
        end

        [p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper ] = rv2coe (r1, v1tu);
    end
    if (methodType == "lambertb" || (runopt == "all" || (dein == de && nrev == 0)))
        [v1tb, v2tb, hitearth, errorsum, errorout] = lambertb(r1, r2, v1, dm, de, 0, 0.0, dtsec, altpadc, 'y');

        %fprintf(fid,'detailSum);
        fprintf(fida,'batt v1t %15.11f  %15.11f  %15.11f\n', v1tb(1),  v1tb(2),  v1tb(3));
        fprintf(fidas,'batt v2t %15.11f  %15.11f  %15.11f\n', v2tb(1),  v2tb(2),  v2tb(3));
        [r3h, v3h] = kepler(r1, v1tb, dtsec);
        for j=1:3
            dr(j) = r2(j) - r3h(j);
        end
        if (mag(dr) > 0.05)
            fprintf(fida,'dm, de nnrev %3d velb does not get to r2 position %15.11f km \n\n', nrev, mag(dr));
            fprintf(fidas,'dm, de nnrev %3d velb does not get to r2 position %15.11f km \n\n', nrev, mag(dr));
        end
        %fprintf(fid,'diffs " + mag(dr).ToString("0.00000000000"));

        for j=1:3
            dv1(j) = v1tk(j) - v1tb(j);
            dv2(j) = v2tk(j) - v2tb(j);
        end
        if (mag(dv1) > 0.01 || mag(dv2) > 0.01)
            fprintf(fida,'velk does not match velb \n\n');
        end
        %fprintf(fid,'diffs " + MathTimeLib::mag(dr).ToString("0.00000000000"));

        [p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper ] = rv2coe (r1, v1tb);
    end



    %watch.Stop(fid);
    %var elapsedMs = watch.ElapsedMilliseconds;
    %Console.WriteLine(watch.ElapsedMilliseconds);

    % use random nrevs, but check if nrev = 0 and set to 1
    % but then you have to check that there is enough time for 1 rev
    int nnrev = nrev;
    if (nnrev == 0)
        nnrev = 1;
    end

    %AstroLibr.lambertminT(r1, r2, 'S', 'L', nnrev, out tmin, out tminp, out tminenergy);
    %fprintf(fid,'mint S " + tmin.ToString("0.0000") + " minp " + tminp.ToString("0.0000") + " minener " + tminenergy.ToString("0.0000"));
    %AstroLibr.lambertminT(r1, r2, 'L', 'L', nnrev, out tmin, out tminp, out tminenergy);
    %fprintf(fid,'mint L " + tmin.ToString("0.0000") + " minp " + tminp.ToString("0.0000") + " minener " + tminenergy.ToString("0.0000"));

    fprintf(fida,' TEST ------------------ S  L " + nnrev, " rev ------------------\n');
    dm = 'S';
    de = 'L';
    if (methodType == "lambertk" || (runopt == "all" || (dmin == dm && dein == de && nrev == nnrev)))
        [kbi, tof] = lambertkmins(s, tau, nnrev, dm, de);
        [v1tk, v2tk, errorsum, errorout] = lambertk ( r1, r2, v1, dm, de, nnrev, 0.0, dtsec, tof, show, 'n' )

        fprintf(fida,' %s\n',detailSum);
        fprintf(fidas,' %s\n', errorout);
        if (~detailAll.Contains("not enough time"))
            %fprintf(fid,'k#" + caseopt  detailSum + " diffs " + mag(dr).ToString("0.00000000000"));
            fprintf(fid,'lamk v1t %15.11f  %15.11f  %15.11f\n', v1tk(1),  v1tk(2),  v1tk(3));
            fprintf(fid,'lamk v2t %15.11f  %15.11f  %15.11f\n', v2tk(1),  v2tk(2),  v2tk(3));
            %fprintf(fid,'magv1t.ToString("0.0000000").PadLeft(12)  magv2t.ToString("0.0000000").PadLeft(12));

            [r3h, v3h] = kepler(r1, v1tk, dtsec);
            for j=1:3
                dr(j) = r2(j) - r3h(j);
            end
            if (mag(dr) > 0.05)
                fprintf(fid,'dm, de nnrev %3d velb does not get to r2 position %15.11f km \n\n', nnrev, mag(dr));
                fprintf(2,'dm, de nnrev %3d velb does not get to r2 position %15.11f km \n\n', nnrev, mag(dr));
            end

            [p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper ] = rv2coe (r1, v1tk);
            tmpstr = sprintf(' %s  %15.11f %15.11f  %15.11f\n', hitearth, mag(dr), a, ecc);
            casearr = strrep(errorout,'\s+', tmpstr)
        end
    else
        fprintf(fida,'dm  %c de %c nrev %d %s\n', dm, de, nnrev, detailAll);
    end
    if (methodType == "lambertu" || (runopt == "all" || (dmin == dm && dein == de && nrev == nnrev)))
        [kbi, tof] = lambertumins( r1, r2, nnrev, dm ) ;   
        [v1tu, v2tu, errorl] = lambertu(r1, r2, v1, dm, de, nnrev, dtsec, kbi, fida );
        %fprintf(fid,'detailSum);
        fprintf(fida,'univ v1t %15.11f  %15.11f  %15.11f\n', v1tu(1),  v1tu(2),  v1tu(3));
        fprintf(fidas,'univ v2t %15.11f  %15.11f  %15.11f\n', v2tu(1),  v2tu(2),  v2tu(3));
        [r3h, v3h] = kepler(r1, v1tu, dtsec);
        for j=1:3
            dr(j) = r2(j) - r3h(j);
        end
        if (mag(dr) > 0.05)
            fprintf(fida,'dm, de nnrev %3d velb does not get to r2 position %15.11f km \n\n', nnrev, mag(dr));
            fprintf(fidas,'dm, de nnrev %3d velb does not get to r2 position %15.11f km \n\n', nnrev, mag(dr));
        end

        for j=1:3
            dv1(j) = v1tk(j) - v1tu(j);
            dv2(j) = v2tk(j) - v2tu(j);
        end
        if (mag(dv1) > 0.01 || mag(dv2) > 0.01)
            fprintf(fida,'velk does not match velu \n\n');
        end

        [p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper ] = rv2coe (r1, v1tu);
    end
    if (methodType == "lambertb" || (runopt == "all" || (dmin == dm && dein == de && nrev == nnrev)))
        [v1tb, v2tb, hitearth, errorsum, errorout] = lambertb(r1, r2, v1, dm, de, nnrev, 0.0, dtsec, altpadc, 'y');
        %fprintf(fid,'detailSum);
        fprintf(fida,'batt v1t %15.11f  %15.11f  %15.11f\n', v1tb(1),  v1tb(2),  v1tb(3));
        fprintf(fidas,'batt v2t %15.11f  %15.11f  %15.11f\n', v2tb(1),  v2tb(2),  v2tb(3));
        [r3h, v3h] = kepler(r1, v1tb, dtsec);
        for j=1:3
            dr(j) = r2(j) - r3h(j);
        end
        if (mag(dr) > 0.05)
            fprintf(fida,'dm, de nnrev %3d velb does not get to r2 position %15.11f km \n\n', nnrev, mag(dr));
            fprintf(fidas,'dm, de nnrev %3d velb does not get to r2 position %15.11f km \n\n', nnrev, mag(dr));
        end
        %fprintf(fid,'diffs " + mag(dr).ToString("0.00000000000"));

        for j=1:3
            dv1(j) = v1tk(j) - v1tb(j);
            dv2(j) = v2tk(j) - v2tb(j);
        end
        if (mag(dv1) > 0.01 || mag(dv2) > 0.01)
            fprintf(fida,'velk does not match velb \n\n');
        end
        %fprintf(fid,'diffs " + MathTimeLib::mag(dr).ToString("0.00000000000"));

        [p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper ] = rv2coe (r1, v1tb);
    end

    fprintf(fida,' TEST ------------------ L  L " + nnrev, " rev ------------------\n');
    dm = 'L';
    de = 'L';
    % switch tdi!!  tdidk to tdirk 'H'
    if (methodType == "lambertk" || (runopt == "all" || (dmin == dm && dein == de && nrev == nnrev)))
        [kbi, tof] = lambertkmins(s, tau, nnrev, dm, de);  % 'H'

        % tofk1, kbik2, tofk2, kbik1;
        %string outstr;
        %getmins(1, 'k', nrev, r1, r2, s, tau, dm, de, out tofk1, out kbik1, out tofk2, out kbik2, out outstr);
        [v1tk, v2tk, errorsum, errorout] = lambertk ( r1, r2, v1, dm, de, nnrev, 0.0, dtsec, tof, show, 'n' )
        fprintf(fida,' %s\n',detailSum);
        fprintf(fidas,' %s\n', errorout);
        if (~detailAll.Contains("not enough time"))
            %fprintf(fid,'k#" + caseopt  detailSum + " diffs " + mag(dr).ToString("0.00000000000"));
            fprintf(fida,'lamk v1t %15.11f  %15.11f  %15.11f\n', v1tk(1),  v1tk(2),  v1tk(3));
            fprintf(fidas,'lamk v2t %15.11f  %15.11f  %15.11f\n', v2tk(1),  v2tk(2),  v2tk(3));
            %fprintf(fid,'magv1t.ToString("0.0000000").PadLeft(12)  magv2t.ToString("0.0000000").PadLeft(12));

            [r3h, v3h] = kepler(r1, v1tk, dtsec);
            for j=1:3
                dr(j) = r2(j) - r3h(j);
            end
            if (mag(dr) > 0.05)
                fprintf(fida,'dm, de nnrev %3d velb does not get to r2 position %15.11f km \n\n', nnrev, mag(dr));
                fprintf(fidas,'dm, de nnrev %3d velb does not get to r2 position %15.11f km \n\n', nnrev, mag(dr));
            end

            [p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper ] = rv2coe (r1, v1tk);
            tmpstr = sprintf(' %s  %15.11f %15.11f  %15.11f\n', hitearth, mag(dr), a, ecc);
            casearr = strrep(errorout,'\s+', tmpstr)
        end
    else
        fprintf(fida,'dm  %c de %c nrev %d %s\n', dm, de, nnrev, detailAll);
        if (methodType == "lambertu" || (runopt == "all" || (dmin == dm && dein == de && nrev == nnrev)))
            [kbi, tof] = lambertumins( r1, r2, nnrev, dm ) ;   
        [v1tu, v2tu, errorl] = lambertu(r1, r2, v1, dm, de, nnrev, dtsec, kbi, fida );
            %fprintf(fid,'detailSum);
            fprintf(fida,'univ v1t %15.11f  %15.11f  %15.11f\n', v1tu(1),  v1tu(2),  v1tu(3));
            fprintf(fidas,'univ v2t %15.11f  %15.11f  %15.11f\n', v2tu(1),  v2tu(2),  v2tu(3));
            [r3h, v3h] = kepler(r1, v1tu, dtsec);
            for j=1:3
                dr(j) = r2(j) - r3h(j);
            end
            if (mag(dr) > 0.05)
                fprintf(fida,'dm, de nnrev %3d velb does not get to r2 position %15.11f km \n\n', nnrev, mag(dr));
                fprintf(fidas,'dm, de nnrev %3d velb does not get to r2 position %15.11f km \n\n', nnrev, mag(dr));
            end

            for j=1:3
                dv1(j) = v1tk(j) - v1tu(j);
                dv2(j) = v2tk(j) - v2tu(j);
            end
            if (mag(dv1) > 0.01 || mag(dv2) > 0.01)
                fprintf(fida,'velk does not match velu \n\n');
            end

            [p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper ] = rv2coe (r1, v1tu);
        end
        if (methodType == "lambertb" || (runopt == "all" || (dmin == dm && dein == de && nrev == nnrev)))
            [v1tb, v2tb, hitearth, errorsum, errorout] = lambertb(r1, r2, v1, dm, de, nnrev, 0.0, dtsec, altpadc, 'y');
            %fprintf(fid,'detailSum);
            fprintf(fida,'batt v1t %15.11f  %15.11f  %15.11f\n', v1tb(1),  v1tb(2),  v1tb(3));
            fprintf(fidas,'batt v2t %15.11f  %15.11f  %15.11f\n', v2tb(1),  v2tb(2),  v2tb(3));
            [r3h, v3h] = kepler(r1, v1tb, dtsec);
            for j=1:3
                dr(j) = r2(j) - r3h(j);
            end
            if (mag(dr) > 0.05)
                fprintf(fida,'dm, de nnrev %3d velb does not get to r2 position %15.11f km \n\n', nnrev, mag(dr));
                fprintf(fidas,'dm, de nnrev %3d velb does not get to r2 position %15.11f km \n\n', nnrev, mag(dr));
            end
            %fprintf(fid,'diffs " + mag(dr).ToString("0.00000000000"));

            for j=1:3
                dv1(j) = v1tk(j) - v1tb(j);
                dv2(j) = v2tk(j) - v2tb(j);
            end
            if (mag(dv1) > 0.01 || mag(dv2) > 0.01)
                fprintf(fida,'velk does not match velb \n\n');
            end
            %fprintf(fid,'diffs " + MathTimeLib::mag(dr).ToString("0.00000000000"));

            [p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper ] = rv2coe (r1, v1tb);
        end
    end

    fprintf(fida,' TEST ------------------ S  H " + nnrev, " rev ------------------\n');
    dm = 'S';
    de = 'H';
    % switch tdi!!  tdirk to tdidk  'L'
    if (methodType == "lambertk" || (runopt == "all" || (dmin == dm && dein == de && nrev == nnrev)))
        [kbi, tof] = lambertkmins(s, tau, nnrev, dm, de);  % 'L'
        [v1tk, v2tk, errorsum, errorout] = lambertk ( r1, r2, v1, dm, de, nnrev, 0.0, dtsec, tof, show, 'n' )
        fprintf(fid,' %s\n',detailSum);
        fprintf(fida,' %s\n', errorout);
        if (~detailAll.Contains("not enough time"))
            %fprintf(fid,'k#" + caseopt  detailSum + " diffs " + mag(dr).ToString("0.00000000000"));
            fprintf(fida,'lamk v1t %15.11f  %15.11f  %15.11f\n', v1tk(1),  v1tk(2),  v1tk(3));
            fprintf(fidas,'lamk v2t %15.11f  %15.11f  %15.11f\n', v2tk(1),  v2tk(2),  v2tk(3));
            %fprintf(fid,'magv1t.ToString("0.0000000").PadLeft(12)  magv2t.ToString("0.0000000").PadLeft(12));

            [r3h, v3h] = kepler(r1, v1tk, dtsec);
            for j=1:3
                dr(j) = r2(j) - r3h(j);
            end
            if (mag(dr) > 0.05)
                fprintf(fida,'dm, de nnrev %3d velb does not get to r2 position %15.11f km \n\n', nnrev, mag(dr));
                fprintf(fidas,'dm, de nnrev %3d velb does not get to r2 position %15.11f km \n\n', nnrev, mag(dr));
            end

            [p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper ] = rv2coe (r1, v1tk);
            tmpstr = sprintf(' %s  %15.11f %15.11f  %15.11f\n', hitearth, mag(dr), a, ecc);
            casearr = strrep(errorout,'\s+', tmpstr)
        end
    else
        fprintf(fida,'dm  %c de %c nrev %d %s\n', dm, de, nnrev, detailAll);
        if (methodType == "lambertu" || (runopt == "all" || (dmin == dm && dein == de && nrev == nnrev)))
            [kbi, tof] = lambertumins( r1, r2, nnrev, dm ) ;   
        [v1tu, v2tu, errorl] = lambertu(r1, r2, v1, dm, de, nnrev, dtsec, kbi, fida );
            %fprintf(fid,'detailSum);
            fprintf(fida,'univ v1t %15.11f  %15.11f  %15.11f\n', v1tu(1),  v1tu(2),  v1tu(3));
            fprintf(fidas,'univ v2t %15.11f  %15.11f  %15.11f\n', v2tu(1),  v2tu(2),  v2tu(3));
            [r3h, v3h] = kepler(r1, v1tu, dtsec);
            for j=1:3
                dr(j) = r2(j) - r3h(j);
            end
            if (mag(dr) > 0.05)
                fprintf(fida,'dm, de nnrev %3d velb does not get to r2 position %15.11f km \n\n', nnrev, mag(dr));
                fprintf(fidas,'dm, de nnrev %3d velb does not get to r2 position %15.11f km \n\n', nnrev, mag(dr));
            end

            for j=1:3
                dv1(j) = v1tk(j) - v1tu(j);
                dv2(j) = v2tk(j) - v2tu(j);
            end
            if (mag(dv1) > 0.01 || mag(dv2) > 0.01)
                fprintf(fida,'velk does not match velu \n\n');
            end

            [p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper ] = rv2coe (r1, v1tu);
        end
        if (methodType == "lambertb" || (runopt == "all" || (dmin == dm && dein == de && nrev == nnrev)))
            [v1tb, v2tb, hitearth, errorsum, errorout] = lambertb(r1, r2, v1, dm, de, nnrev, 0.0, dtsec, altpadc, 'y');
            %fprintf(fid,'detailSum);
            fprintf(fida,'batt v1t %15.11f  %15.11f  %15.11f\n', v1tb(1),  v1tb(2),  v1tb(3));
            fprintf(fidas,'batt v2t %15.11f  %15.11f  %15.11f\n', v2tb(1),  v2tb(2),  v2tb(3));
            [r3h, v3h] = kepler(r1, v1tb, dtsec);
            for j=1:3
                dr(j) = r2(j) - r3h(j);
            end
            if (mag(dr) > 0.05)
                fprintf(fida,'dm, de nnrev %3d velb does not get to r2 position %15.11f km \n\n', nnrev, mag(dr));
                fprintf(fidas,'dm, de nnrev %3d velb does not get to r2 position %15.11f km \n\n', nnrev, mag(dr));
            end
            %fprintf(fid,'diffs " + mag(dr).ToString("0.00000000000"));

            for j=1:3
                dv1(j) = v1tk(j) - v1tb(j);
                dv2(j) = v2tk(j) - v2tb(j);
            end
            if (mag(dv1) > 0.01 || mag(dv2) > 0.01)
                fprintf(fida,'velk does not match velb \n\n');
            end
            %fprintf(fid,'diffs " + MathTimeLib::mag(dr).ToString("0.00000000000"));

            [p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper ] = rv2coe (r1, v1tb);
        end
    end

    fprintf(fida,' TEST ------------------ L  H " + nnrev, " rev ------------------\n');
    dm = 'L';
    de = 'H';
    if (methodType == "lambertk" || (runopt == "all" || (dmin == dm && dein == de && nrev == nnrev)))
        [kbi, tof] = lambertkmins(s, tau, nnrev, dm, de);
        [v1tk, v2tk, errorsum, errorout] = lambertk ( r1, r2, v1, dm, de, nnrev, 0.0, dtsec, tof, show, 'n' )
        fprintf(fida,' %s\n',detailSum);
        fprintf(fidsa,' %s\n', errorout);
        if (~detailAll.Contains("not enough time"))
            %fprintf(fid,'k#" + caseopt  detailSum + " diffs " + mag(dr).ToString("0.00000000000"));
            fprintf(fida,'lamk v1t %15.11f  %15.11f  %15.11f\n', v1tk(1),  v1tk(2),  v1tk(3));
            fprintf(fidas,'lamk v2t %15.11f  %15.11f  %15.11f\n', v2tk(1),  v2tk(2),  v2tk(3));
            %fprintf(fid,'magv1t.ToString("0.0000000").PadLeft(12)  magv2t.ToString("0.0000000").PadLeft(12));

            [r3h, v3h] = kepler(r1, v1tk, dtsec);
            for j=1:3
                dr(j) = r2(j) - r3h(j);
            end
            if (mag(dr) > 0.05)
                fprintf(fida,'dm, de nnrev %3d velb does not get to r2 position %15.11f km \n\n', nnrev, mag(dr));
                fprintf(2,'dm, de nnrev %3d velb does not get to r2 position %15.11f km \n\n', nnrev, mag(dr));
            end

            [p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper ] = rv2coe (r1, v1tk);
            tmpstr = sprintf(' %s  %15.11f %15.11f  %15.11f\n', hitearth, mag(dr), a, ecc);
            casearr = strrep(errorout,'\s+', tmpstr)
        end
    else
        fprintf(fida,'dm  %c de %c nrev %d %s\n', dm, de, nnrev, detailAll);
    end
    if (methodType == "lambertu" || (runopt == "all" || (dmin == dm && dein == de && nrev == nnrev)))
        [kbi, tof] = lambertumins( r1, r2, nnrev, dm ) ;
        [v1tu, v2tu, errorl] = lambertu(r1, r2, v1, dm, de, nnrev, dtsec, kbi, fida );
        %fprintf(fid,'detailSum);
        fprintf(fida,'univ v1t %15.11f  %15.11f  %15.11f\n', v1tu(1),  v1tu(2),  v1tu(3));
        fprintf(fidas,'univ v2t %15.11f  %15.11f  %15.11f\n', v2tu(1),  v2tu(2),  v2tu(3));
        [r3h, v3h] = kepler(r1, v1tu, dtsec);
        for j=1:3
            dr(j) = r2(j) - r3h(j);
            if (mag(dr) > 0.05)
                fprintf(fida,'dm, de nnrev %3d velb does not get to r2 position %15.11f km \n\n', nnrev, mag(dr));
                fprintf(2,'dm, de nnrev %3d velb does not get to r2 position %15.11f km \n\n', nnrev, mag(dr));
            end

            for j=1:3
                dv1(j) = v1tk(j) - v1tu(j);
                dv2(j) = v2tk(j) - v2tu(j);
            end
            if (mag(dv1) > 0.01 || mag(dv2) > 0.01)
                fprintf(fida,'velk does not match velu \n\n');
            end

            [p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper ] = rv2coe (r1, v1tu);
        end
        if (methodType == "lambertb" || (runopt == "all" || (dmin == dm && dein == de && nrev == nnrev)))
            [v1tb, v2tb, hitearth, errorsum, errorout] = lambertb(r1, r2, v1, dm, de, nnrev, 0.0, dtsec, altpadc, 'y');

            %fprintf(fid,' %s\n',detailSum);
            fprintf(fida,'batt v1t %15.11f  %15.11f  %15.11f\n', v1tb(1),  v1tb(2),  v1tb(3));
            fprintf(fidas,'batt v2t %15.11f  %15.11f  %15.11f\n', v2tb(1),  v2tb(2),  v2tb(3));
            [r3h, v3h] = kepler(r1, v1tb, dtsec);
            for j=1:3
                dr(j) = r2(j) - r3h(j);
                if (mag(dr) > 0.05)
                    fprintf(fida,'dm, de nnrev %3d velb does not get to r2 position %15.11f km \n\n', nnrev, mag(dr));
                    fprintf(fidas,'dm, de nnrev %3d velb does not get to r2 position %15.11f km \n\n', nnrev, mag(dr));
                end
                %fprintf(fid,'diffs " + mag(dr).ToString("0.00000000000"));

                for j=1:3
                    dv1(j) = v1tk(j) - v1tb(j);
                    dv2(j) = v2tk(j) - v2tb(j);
                end
                if (mag(dv1) > 0.01 || mag(dv2) > 0.01)
                    fprintf(fida,'velk does not match velb \n\n');
                end
                %fprintf(fid,'diffs " + MathTimeLib::mag(dr).ToString("0.00000000000"));

                [p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper ] = rv2coe (r1, v1tb);

            end  % dolamberttests
        end
    end

end  % dolamberttests



% ------------------------------------------------------------------------------
%                                    getmins
%
% find tbi mins for k, lambert, etc. does either universal and k approach.
%   universal seems to find them better though... could store in an array, but
%   faster to do 1 at a time. also note that the k lambert needs SH - LL, and LL - SH
%   reversal and it is done here.
%
%  author        : david vallado             davallado@gmail.com  10 oct 2019
%
%  inputs        description                                   range / units
%    loopktr     - ktr for whether or not to write output           0 (writes)
%    app         - which approach to use                           'u', 'k'
%    nrev        - number of revolultions
%    r1          - first position vector                            km
%    r2          - second position vector                           km
%    s           - parameter for k only, not needed for univ
%    tau         - parameter for k only, not needed for univ
%    dm          - parameter for k only, not needed for univ        'S', 'L'
%    de          - parameter for k only, not needed for univ        'L', 'H'
%
%  outputs       :
%    tof         - minute tof for the specified nrev                    s
%    kbi         - minute psi/k/etc for the given nrev
%    outstr      - output string if case 0
%
%  locals :
%
%  coupling      :
% [tof1, kbi1, tof2, kbi2, outstr] = getmins (loopktr, app, nrev, r1, r2, s, tau, dm, de)
% ------------------------------------------------------------------------------*/

% function [tof1, kbi1, tof2, kbi2, outstr] = getmins (loopktr, app, nrev, r1, r2, s, tau, dm, de)
% 
%     tof1 = 0.0;
%     kbi1 = 0.0;
%     tof2 = 0.0;
%     kbi2 = 0.0;
% 
%     %tusec = 806.8111238242922;
%     if (nrev > 0)
% 
% 
%         if (app == 'u')
% 
%             % universal variable approach
%             [kbi1, tof1] = lambertumins(r1, r2, nrev, 'S');
%             [kbi2, tof2] = lambertumins(r1, r2, nrev, 'L');
%         end
%     else
% 
%         % -----------do these calcs one time to save time
%         % call this outside getmins
%         %     [s, tau] = lambertkmins1st(r1, r2);
% 
%         % k value approaches
%         [kbi1, tof1] = lambertkmins(s, tau, nrev, 'x', 'L');
%         %tof1 = tof1 / tusec;  % in tu, tof is in secs
%         [kbi2, tof2] = lambertkmins(s, tau, nrev, 'x', 'H');
%         %tof2 = tof2 / tusec;
%         % switch these here so it's not needed elsewhere
%         % no, switch in kmins
%         %if ((dm == 'S' && de == 'H') || ((dm == 'L' && de == 'L')))
%         %
%         %    temp = kbi1;
%         %    kbi1 = kbi2;
%         %    kbi2 = temp;
%         %    temp = tof1;
%         %    tof1 = tof2;
%         %    tof2 = temp;
%         %end
%     end
% 
%     % -----------------------------put minute values etc points on plot ---------------------------------
%     % writeout just one time
%     if (loopktr == 0)
% 
%         fprintf(fid,'Lambertumin');
%         detailSum = 'S   L   1  0.000 ' + tof1 + ' psimin ' + kbi1 + ' -0\n';
%         fprintf(fid,detailSum);
% 
%         fprintf(fid,'Lambertkmin');
%         detailSum = 'S   L   1  0.000 ' + (tof1) + ' kmin ' + kbi1 + ' -0\n';
%         fprintf(fid,detailSum);
% 
%         % -----------------------------put minute values etc points on plot ---------------------------------
%         fprintf(fid,'Lamberttmin');
%         [tmin, tminp, tminenergy] = lambertminT(r1, r2, 'S', 'L', 1);
%         detailSum = 'S   L   tminp  0.000 ' + tminp + '5.0000000' ' -0';
%         fprintf(fid,detailSum);
%         detailSum = 'S   L   tminp  0.000 ' + tminp + '-5.0000000' ' -0\n';
%         fprintf(fid,detailSum);
% 
%         detailSum = 'S   L   tminenergy  0.000 ' + tminenergy + '5.0000000' ' -0';
%         fprintf(fid,detailSum);
%         detailSum = 'S   L   tminenergy  0.000 ' + tminenergy + '-5.0000000' ' -0\n';
%         fprintf(fid,detailSum);
% 
%         [kbish, tofsh] = lambertumins(r1, r2, 1, 'S');
%         detailSum = 'S   L   1  0.000 ' + tmin + (tofsh - 5) + ' -0';
%         fprintf(fid,detailSum);
%         detailSum = 'S   L   1  0.000 ' + tmin + (tofsh + 5) + ' -0\n';
%         fprintf(fid,detailSum);
%         [kbish, tofsh] = lambertumins(r1, r2, 2, 'S');
%         [tmin, tminp, tminenergy] = lambertminT(r1, r2, 'S', 'L', 2);
%         detailSum = 'S   L   2  0.000 ' + tmin + (tofsh - 5) + ' -0';
%         fprintf(fid,detailSum);
%         detailSum = 'S   L   2  0.000 ' + tmin + (tofsh + 5) + ' -0\n';
%         fprintf(fid,detailSum);
%         [kbish, tofsh] = lambertumins(r1, r2, 3, 'S');
%         [tmin, tminp, tminenergy] = lambertminT(r1, r2, 'S', 'L', 3);
%         detailSum = 'S   L   3  0.000 ' + tmin + (tofsh - 5) + ' -0';
%         fprintf(fid,detailSum);
%         detailSum = 'S   L   3  0.000 ' + tmin + (tofsh + 5) + ' -0\n';
%         fprintf(fid,detailSum);
% 
%         [tmin, tminp, tminenergy] = lambertminT(r1, r2, 'L', 'H', 1);
%         detailSum = 'L   H   tminp  0.000 ' + tminp + '5.0000000' ' -0';
%         fprintf(fid,detailSum);
%         detailSum = 'L   H   tminp  0.000 ' + tminp + '-5.0000000' ' -0\n';
%         fprintf(fid,detailSum);
% 
%         detailSum = 'L   H   tminenergy  0.000 ' + tminenergy + '5.0000000' ' -0';
%         fprintf(fid,detailSum);
%         detailSum = 'L   H   tminenergy  0.000 ' + tminenergy + '-5.0000000' ' -0\n';
%         fprintf(fid,detailSum);
% 
%         [kbilg, toflg] = lambertumins(r1, r2, 1, 'L');
%         detailSum = 'L   H   1  0.000 ' + tmin + (toflg - 5) + ' -0';
%         fprintf(fid,detailSum);
%         detailSum = 'L   H   1  0.000 ' + tmin + (toflg + 5) + ' -0\n';
%         fprintf(fid,detailSum);
%         [kbilg, toflg] = lambertumins(r1, r2, 2, 'L');
%         [tmin, tminp, tminenergy] = lambertminT(r1, r2, 'S', 'H', 2);
%         detailSum = 'L   H   2  0.000 ' + tmin + (toflg - 5) + ' -0';
%         fprintf(fid,detailSum);
%         detailSum = 'L   H   2  0.000 ' + tmin + (toflg + 5) + ' -0\n';
%         fprintf(fid,detailSum);
%         [kbilg, toflg] = lambertumins(r1, r2, 3, 'L');
%         [tmin, tminp, tminenergy] = lambertminT(r1, r2, 'S', 'H', 3);
%         detailSum = 'L   H   3  0.000 ' + tmin + (toflg - 5) + ' -0';
%         fprintf(fid,detailSum);
%         detailSum = 'L   H   3  0.000 ' + tmin + (toflg + 5) + ' -0\n';
%         fprintf(fid,detailSum);
% 
%         % find max rp values for each nrev
%          tmaxrp;
%         [tmaxrp, v1t] = lambertTmaxrp(r1, r2, 'S', 0);
%         fprintf(fid,'x   x   0  0.000 ' + tmaxrp,
%         v1t(1), v1t(2), v1t(3).ToString('0.0000000').PadLeft(15));
%         fprintf(fid,'x   x   0  0.000 ' + tmaxrp,
%         v1t(1), v1t(2), v1t(3), '\n');
%         [tmaxrp, v1t] = lambertTmaxrp(r1, r2, 'S', 1);
%         fprintf(fid,'x   x   1  0.000 ' + tmaxrp,
%         v1t(1), v1t(2), v1t(3).ToString('0.0000000').PadLeft(15));
%         fprintf(fid,'x   x   1  0.000 ' + tmaxrp,
%         v1t(1), v1t(2), v1t(3), '\n');
%         [tmaxrp, v1t] = lambertTmaxrp(r1, r2, 'S', 2);
%         fprintf(fid,'x   x   2  0.000 ' + tmaxrp,
%         v1t(1), v1t(2), v1t(3).ToString('0.0000000').PadLeft(15));
%         fprintf(fid,'x   x   2  0.000 ' + tmaxrp,
%         v1t(1), v1t(2), v1t(3), '\n');
%     end  % if nrev > 0
% 
%     outstr = strbuild;
% end   % getmins
% 
% 
% % ------------------------------------------------------------------------------
% %                                      makesurf
% %
% %  make a surface from a fixdat result with numbers takes a text file of number of points,
% %  then all the points, and cross-hatches it to get a surface. you run fixdat first
% %  in most cases.
% %
% %  author        : david vallado             davallado@gmail.com  10 oct 2019
% %
% %  inputs        description                                   range / units
% %    infilename  - in filename
% %    outfilename - out filename
% %
% %  outputs       :
% %
% %  locals :
% %
% %  coupling      :
% %
% % ------------------------------------------------------------------------------*/
% 
% function makesurf( infilename, outfilename )
%     Restoflgine = '';
% 
%     fileData = File.ReadAllLines(infilename);
% 
%     % process all the x lines
%     numPerLine = 0;
%     ktr = 0;     % reset the file
%     NumLines = 0;
%     while (ktr < fileData.Count(fid) - 1)  % not eof
% 
%         line = fileData[ktr];
%         linesplt = line.Split(' ');
%         numPerLine = str2num(linesplt(1));
%         % matlab uses Inf or Nan to start a new line
%         % needs to be in each col as well
%         fprintf(fid,' Nan Nan NaN NaN NaN NaN');  % numPerLine,
%         ktr = ktr + 1;
%         NumLines = NumLines + 1;
% 
%         for i = 1 : numPerLine
% 
%             line = fileData[ktr];
%             line1 = Regex.Replace(line, @'\s+', ' ');
%             linesplt = line1.Split(' ');
%             %posrest = line1.IndexOf(linesplt(3), 15); % start at position 3
%             %Restoflgine = line1.Substring(posrest -1, line1.Length -posrest);
%             Restoflgine = linesplt(3), linesplt(4),
%             linesplt(5), linesplt(6);
%             fprintf(fid,linesplt(1), linesplt(2), Restoflgine);
%             ktr = ktr + 1;
%         end
%     end
% 
%     % ------process the y lines-------
%     % 'process y lines ');
%     % the number of lines needs to be constant!!
%     %numPerLine = 0;
%     numinrow = 0;  % position of each y line
%     % go through each poin initial line
%     while (numinrow < numPerLine)
% 
%         % ---get the nth pofrom the first row---
%         ktr = 0;  % reset the file
%         %line = fileData[ktr];
%         %linesplt = line.Split(' ');
%         %k = str2num(linesplt(1));
%         fprintf(fid,' Nan Nan NaN NaN NaN NaN');  % k,
%         ktr = ktr + 1 + numinrow;  % get to first poof data
%         line = fileData[ktr];
%         line1 = Regex.Replace(line, '\s+', ' ');
%         linesplt = line1.Split(' ');
%         Restoflgine = linesplt(3), linesplt(4),
%         linesplt(5), linesplt(6);
%         fprintf(fid,linesplt(1), linesplt(2), Restoflgine);
% 
%         % ---get nth number from each other segment---
%         % since they are all evenly spaced, simply add the delta until the end of file
%         Int32 ktr0 = ktr + 1;
%         for (j = 1; j < NumLines; j++)
% 
%             ktr = ktr0 + j * numPerLine;
%             line = fileData[ktr];
%             line1 = Regex.Replace(line, '\s+', ' ');
%             linesplt = line1.Split(' ');
%             Restoflgine = linesplt(3), linesplt(4),
%             linesplt(5), linesplt(6);
%             fprintf(fid,linesplt(1), linesplt(2), Restoflgine);
%             ktr0 = ktr0 + 1;
%         end
% 
%         numinrow = numinrow + 1;
%     end  % while
% 
%     directory = 'd:\codes\library\matlab\';
%     File.WriteAllText(directory + 'surf.out', strbuild);
% end  % makesurf
% 
% 
% 
% % ------------------------------------------------------------------------------
% %                                      fixdat
% %
% %  fix the blank lines in a datafile. let 4 values be taken depending on the
% %  intindex values inserts the number of points for each segment, then it can be
% %  used in makesurf.
% %
% %  author        : david vallado             davallado@gmail.com  10 oct 2019
% %
% %  inputs        description                                   range / units
% %    infilename  - in filename
% %    outfilename - out filename
% %    intindxx    - which indices to use, 1.2.3.4.5, 2.5.7.8 etc
% %
% %  outputs       :
% %
% %  locals :
% %
% %  ------------------------------------------------------------------------------*/
% 
% function fixdat( infilename,outfilename,intindx1, intindx2, intindx3, intindx4, intindx5, intindx6 )
%     fileData = File.ReadAllLines(infilename);
% 
%     i = 0;
%     ktr = 0;
%     while (ktr < fileData.Count(fid) - 1)  % not eof
% 
%         LongString = fileData[ktr];
% 
%         if ((LongString.Contains('xx')) || LongString.Length < 10 || (i == 2000))
% 
%             % ----Put a mandatory break at 2000----
%             if ((i == 2000) && (!LongString.Contains('xx')))
%                 fprintf(fid,(i + 1), ' xx broken');
%             else
% 
%                 if (i > 0)
%                     fprintf(fid,i, ' xx ');
%                 end
%             end
% 
%             % ----Write out all the data to this point----
%             for j = 1:i
%                 fprintf(fid,DatArray(j));
%             end
% 
%             % ----Write out crossover pofor the mandatory break ---
%             if ((i == 2000) && ((!LongString.Contains('x')) || LongString.Length > 10))
%                 fprintf(fid,LongString);
%                 % write out last one in loop !
%                 i = 1;
%                 DatArray(i) = LongString;
%             else
%                 i = 0;
%                 ktr = ktr + 1;   % file counter
%             end
%         end
%     else
%         % ----Add new data ----
% 
%         i = i + 1;
%         string line1 = Regex.Replace(LongString, @'\s+', ' ');
%         string[] linesplt = line1.Split(' ');
% 
%         DatArray(i) = linesplt[intindx1], linesplt[intindx2], linesplt[intindx3], linesplt[intindx4]
%         , linesplt[intindx5], linesplt[intindx6];
%         ktr = ktr + 1;   % file counter
% 
%     end  % while through file
% 
%     % ----Write out the last section----
%     if (i > 1)
% 
%         if ((i == 2000) && (!LongString.Contains('x')) && LongString.Length > 10)
%             fprintf(fid,(i + 1), ' xx broken');
%         else
%             if (i > 0)
%                 fprintf(fid,i, ' xx ');
%             end
%         end
%         for (j = 1; j <= i; j++)
%             fprintf(fid,DatArray(j));
%         end
%     end
% 
%     directory = 'd:\codes\library\matlab\';
%     File.WriteAllText(directory + 't.out', strbuild);
% end  % fixdat


function testlambertuniv(fid)
    constastro;

    errorsum = '';
    errorout = '';
    show = 'n';     % for test180, show = n, show180 = y
    % show180 = 'n';  % for testlamb known show = y, show180 = n, n/n for envelope
    nrev = 0;
    r1 = [ 2.500000 * re; 0.000000; 0.000000 ];
    r2 = [ 1.9151111 * re; 1.6069690 * re; 0.000000 ];
    % assume circular initial orbit for vel calcs
    v1 = [ 0.0; sqrt(mu / r1(1)); 0.0 ];
    ang = atan(r2(2) / r2(1));
    v2 = [ -sqrt(mu / r2(2)) * cos(ang); sqrt(mu / r2(1)) * sin(ang); 0.0 ];
    dtsec = 76.0 * 60.0;
    altpad = 100.0;  % 100 km

    [kbi, tof] = lambertumins( r1, r2, 1, 'S' ) ;
    tbiSu(2, 2) = kbi;
    tbiSu(2, 3) = tof;
    [kbi, tof] = lambertumins( r1, r2, 2, 'S' ) ;
    tbiSu(3, 2) = kbi;
    tbiSu(3, 3) = tof;
    [kbi, tof] = lambertumins( r1, r2, 3, 'S' ) ;
    tbiSu(4, 2) = kbi;
    tbiSu(4, 3) = tof;
    [kbi, tof] = lambertumins( r1, r2, 4, 'S' ) ;
    tbiSu(5, 2) = kbi;
    tbiSu(5, 3) = tof;
    [kbi, tof] = lambertumins( r1, r2, 5, 'S' ) ;
    tbiSu(6, 2) = kbi;
    tbiSu(6, 3) = tof;

    [kbi, tof] = lambertumins( r1, r2, 1, 'L' ) ;
    tbiLu(2, 2) = kbi;
    tbiLu(2, 3) = tof;
    [kbi, tof] = lambertumins( r1, r2, 2, 'L' ) ;
    tbiLu(3, 2) = kbi;
    tbiLu(3, 3) = tof;
    [kbi, tof] = lambertumins( r1, r2, 3, 'L' ) ;
    tbiLu(4, 2) = kbi;
    tbiLu(4, 3) = tof;
    [kbi, tof] = lambertumins( r1, r2, 4, 'L' ) ;
    tbiLu(5, 2) = kbi;
    tbiLu(5, 3) = tof;
    [kbi, tof] = lambertumins( r1, r2, 5, 'L' ) ;
    tbiLu(6, 2) = kbi;
    tbiLu(6, 3) = tof;

    if (show == 'y')
        dtsec = 21000.0;
        fprintf(fid,' TEST ------------------ S  L  0 rev\n');
        [v1t, v2t, errorl] = lambertu ( r1,  r2, v1, 'S', 'L', 0, dtsec, 0.0, fid );
        fprintf(fid,' v1t %16.8f%16.8f%16.8f\n',v1t );
        fprintf(fid,' v2t %16.8f%16.8f%16.8f\n',v2t );

        fprintf(fid,' TEST ------------------ L  H  0 rev\n');
        [v1t, v2t, errorl] = lambertu ( r1,  r2, v1, 'L', 'H', 0, dtsec, 0.0, fid );
        fprintf(fid,' v1t %16.8f%16.8f%16.8f\n',v1t );
        fprintf(fid,' v2t %16.8f%16.8f%16.8f\n',v2t );

        fprintf(fid,' TEST ------------------ S  L  1 rev\n');
        [kbi, tof] = lambertumins( r1, r2, 1, 'S' ) ;
        [v1t, v2t, errorl] = lambertu ( r1,  r2, v1, 'S', 'L', 1, dtsec, kbi, fid );
        fprintf(fid,' v1t %16.8f%16.8f%16.8f\n',v1t );
        fprintf(fid,' v2t %16.8f%16.8f%16.8f\n',v2t );

        fprintf(fid,' TEST ------------------ S  H  1 rev\n');
        [kbi, tof] = lambertumins( r1, r2, 1, 'S' ) ;
        [v1t, v2t, errorl] = lambertu ( r1,  r2, v1, 'S', 'H', 1, dtsec, kbi, fid );
        fprintf(fid,' v1t %16.8f%16.8f%16.8f\n',v1t );
        fprintf(fid,' v2t %16.8f%16.8f%16.8f\n',v2t );

        fprintf(fid,' TEST ------------------ L  L  1 rev\n');
        [kbi, tof] = lambertumins( r1, r2, 1, 'L' ) ;
        [v1t, v2t, errorl] = lambertu ( r1,  r2, v1, 'L', 'L', 1, dtsec, kbi, fid );
        fprintf(fid,' v1t %16.8f%16.8f%16.8f\n',v1t );
        fprintf(fid,' v2t %16.8f%16.8f%16.8f\n',v2t );

        fprintf(fid,' TEST ------------------ L  H  1 rev\n');
        [kbi, tof] = lambertumins( r1, r2, 1, 'L' ) ;
        [v1t, v2t, errorl] = lambertu ( r1,  r2, v1, 'L', 'H', 1, dtsec, kbi, fid );
        fprintf(fid,' v1t %16.8f%16.8f%16.8f\n',v1t );
        fprintf(fid,' v2t %16.8f%16.8f%16.8f\n',v2t );

        fprintf(fid,'\n-------- lambertb test\n' );
        [v1t, v2t, errorl] = lambertb ( r1,  r2, v1, 'L', 'L', 1, dtsec );
        fprintf(fid,' v1dv %16.8f%16.8f%16.8f\n',v1t );
        fprintf(fid,' v2dv %16.8f%16.8f%16.8f\n',v2t );
    end

    % test chap 7 fig 7-18 runs
    if (show == 'y')
        fprintf(fid,' TEST ------------------ s/l  d  0 rev ------------------');
        hitearth = ' ';
        dm = 'S';
        de = 'L';
        nrev = 0;
        dtwait = 0.0;

        % fig 7-18 fixed tgt and int
        r1 = [ -6518.1083; -2403.8479; -22.1722 ];
        v1 = [ 2.604057; -7.105717; -0.263218 ];
        r2 = [ 6697.4756; 1794.5832; 0.0 ];
        v2 = [ -1.962373; 7.323674; 0.0 ];
        fprintf(fid,'dtwait  dtsec       dv1       dv2 ');
        %           for (i = 0; i < 250; i++)

        i = 0;
        dtsec = i * 60.0;
        [kbi, tof] = lambertumins( r1, r2, 1, 'S' ) ;
        [v1t1, v2t1, errorl] = lambertu ( r1,  r2, v1, 'S', 'L', 1, dtsec, kbi, fid );

        [kbi, tof] = lambertumins( r1, r2, 1, 'S' ) ;
        [v1t2, v2t2, errorl] = lambertu ( r1,  r2, v1, 'S', 'H', 1, dtsec, kbi, fid );

        [kbi, tof] = lambertumins( r1, r2, 1, 'L' ) ;
        [v1t3, v2t3, errorl] = lambertu ( r1,  r2, v1, 'L', 'L', 1, dtsec, kbi, fid );

        [kbi, tof] = lambertumins( r1, r2, 1, 'L' ) ;
        [v1t4, v2t4, errorl] = lambertu ( r1,  r2, v1, 'L', 'H', 1, dtsec, kbi, fid );

        if (errorout.Contains('ok'))
            dv11 = v1t1 - v1;
            dv21 = v2t1 - v2;
            dv12 = v1t2 - v1;
            dv22 = v2t2 - v2;
            dv13 = v1t3 - v1;
            dv23 = v2t3 - v2;
            dv14 = v1t4 - v1;
            dv24 = v2t4 - v2;

            fprintf(fid,'%15.11f %15.11f  %15.11f %15.11f %15.11f %15.11f %15.11f %15.11f %15.11f %15.11f\n',dtwait, dtsec,...
                mag(dv11), mag(dv21),  mag(dv12),  mag(dv22),   mag(dv13), mag(dv23),  mag(dv14), mag(dv24));
        end
    else
        %fprintf(fid,errorsum, errorout);
    end


    % fig 7-19 moving tgt
    r1 = [ 5328.7862; 4436.1273; 101.4720 ];
    v1 = [ -4.864779; 5.816486; 0.240163 ];
    r2 = [ 6697.4756; 1794.5831; 0.0 ];
    v2 = [ -1.962372; 7.323674; 0.0 ];
    fprintf(fid,'dtwait  dtsec       dv1       dv2 ');
    % diff vectors, needs new umins

    %           for (i = 0; i < 250; i++)

    i = 0;
    dtsec = i * 60.0;
    [r3, v3] = kepler(r2, v2, dtsec);
    %                for (j = 0; j < 25; j++)

    j = 0;
    dtwait = j * 10.0;
    dtwait = 0.0;  % set to 0 for now

    [kbi, tof] = lambertumins( r1, r3, 1, 'S' ) ;
    [v1t1, v2t1, errorl] = lambertu ( r1,  r3, v1, 'S', 'L', 1, dtsec, kbi, fid );

    [kbi, tof] = lambertumins( r1, r3, 1, 'S' ) ;
    [v1t2, v2t2, errorl] = lambertu ( r1,  r3, v1, 'S', 'H', 1, dtsec, kbi, fid );

    [kbi, tof] = lambertumins( r1, r3, 1, 'L' ) ;
    [v1t3, v2t3, errorl] = lambertu ( r1,  r3, v1, 'L', 'L', 1, dtsec, kbi, fid );

    [kbi, tof] = lambertumins( r1, r3, 1, 'L' ) ;
    [v1t4, v2t4, errorl] = lambertu ( r1,  r3, v1, 'L', 'H', 1, dtsec, kbi, fid );

    %if (errorout.Contains('ok'))
        dv11 = v1t1 - v1;
        dv21 = v2t1 - v2;
        dv12 = v1t2 - v1;
        dv22 = v2t2 - v2;
        dv13 = v1t3 - v1;
        dv23 = v2t3 - v2;
        dv14 = v1t4 - v1;
        dv24 = v2t4 - v2;

        fprintf(fid,'%15.11f %15.11f  %15.11f %15.11f %15.11f %15.11f %15.11f %15.11f %15.11f %15.11f\n',dtwait, dtsec,...
            mag(dv11), mag(dv21),  mag(dv12),  mag(dv22),   mag(dv13), mag(dv23),  mag(dv14), mag(dv24));
    %end

    % fig 7-21
    r1 = [ -6175.1034; 2757.0706; 1626.6556 ];
    v1 = [ 2.376641; 1.139677; 7.078097 ];
    r2 = [ -6078.007289; 2796.641859; 1890.7135 ];
    v2 = [ 2.654700; 1.018600; 7.015400 ];

    fprintf(1,'dtwait  dtsec       dv1       dv2\n');
    totaldts = 15000;
    totaldtw = 30000;
    step1 = 60;   % 60 orig
    step2 = 600;  % 600 orig
    stop1 = floor(totaldts / step1);
    stop2 = floor(totaldtw / step2);

    for i=0 : 2 %stop1  % orig 250, 15000 s total

        dtsec = i * step1;    % orig 60
        [r4, v4] = kepler(r1, v1, dtsec);
        for j = 0: 2 %stop2   % orig 25 600*25 = 15000 s total

            dtwait = j * step2;   % orig 600
            [r3, v3] =  kepler(r2, v2, dtsec + dtwait);
            [kbi, tof] = lambertumins( r3, r4, 1, 'S' ) ;
            [v1t1, v2t1, errorl] = lambertu (r3,  r4, v1, 'S', 'L', 1, dtsec, kbi, fid );

            [kbi, tof] = lambertumins( r3, r4, 1, 'S' ) ;
            [v1t2, v2t2, errorl] = lambertu ( r3,  r4, v1, 'S', 'H', 1, dtsec, kbi, fid );

            [kbi, tof] = lambertumins( r3, r4, 1, 'L' ) ;
            [v1t3, v2t3, errorl] = lambertu ( r3,  r4, v1, 'L', 'L', 1, dtsec, kbi, fid );

            [kbi, tof] = lambertumins( r3, r4, 1, 'L' ) ;
            [v1t4, v2t4, errorl] = lambertu ( r3,  r4, v1, 'L', 'H', 1, dtsec, kbi, fid );

            %if (errorout.Contains('ok'))
                dv11 = v1t1 - v1;
                dv21 = v2t1 - v2;
                dv12 = v1t2 - v1;
                dv22 = v2t2 - v2;
                dv13 = v1t3 - v1;
                dv23 = v2t3 - v2;
                dv14 = v1t4 - v1;
                dv24 = v2t4 - v2;

                fprintf(fid,'%15.11f %15.11f  %15.11f %15.11f %15.11f %15.11f %15.11f %15.11f %15.11f %15.11f\n',dtwait, dtsec,...
                    mag(dv11), mag(dv21),  mag(dv12),  mag(dv22),   mag(dv13), mag(dv23),  mag(dv14), mag(dv24));
           % else
           %     strbuildFig.AppendLine(errorsum, errorout);
           % end
        end

        % write data out
       % string directory = @'D:\Codes\LIBRARY\Matlab\';
       % File.WriteAllText(directory + 'surfMovingSalltest.out', strbuildFig);
    end
end



function testradecgeo2azel(fid)
    rad = 180.0 / pi;
    rtasc = 294.98914583 / rad;
    decl = -20.8234944 / rad;
    xp = 0.0;
    yp = 0.0;
    lod = 0.0;
    jdut1 = 2453101.82740678310;
    ttt = 0.042623631889;
    ddpsi = -0.052195;
    ddeps = -0.003875;
    rr = 12373.3546098;  % km
    latgd = 39.007 / rad;
    lon = -104.883 / rad;
    alt = 0.3253;

    [az, el] = radecgeo2azel(rtasc, decl, rr, latgd, lon, alt, ttt, jdut1, lod, xp, yp, ddpsi, ddeps);
end

function testijk2ll(fid)
    rad = 180.0 / pi;

    r = [ 1.023 * re; 1.076 * re; 1.011 * re ];

    [latgc, latgd, lon, hellp] = ecef2ll(r);

    fprintf(fid,'ecef2ll %15.11f  %15.11f  %15.11f\n', latgd*rad, lon * rad, hellp);

    [latgc, latgd, lon, hellp] = ecef2llb(r);

    fprintf(fid,'ecef2llb %15.11f  %15.11f  %15.11f\n', latgd*rad, lon * rad, hellp);
end

function testgd2gc(fid)
    rad = 180.0 / pi;
    latgd = 34.173429 / rad;

    [ansr] = gd2gc(latgd);

    fprintf(fid,'gd2gc %15.11f\n', ansr);
end

function testsite(fid)
    rad = 180.0 / pi;
    latgd = 39.007 / rad;
    lon = -104.883 / rad;
    alt = 0.3253;

    [rsecef, vsecef] = site(latgd, lon, alt);
    fprintf(fid,'site gc %16.9f %16.9f %16.9f  %16.9f\n', rsecef, latgd*rad );

    fprintf(fid,'site %16.9f %16.9f %16.9f %16.9f %16.9f %16.9f\n', rsecef(1), rsecef(2), rsecef(3), ...
        vsecef(1), vsecef(2), vsecef(3));
end


% --------  sun          - analytical sun ephemeris
function testsun(fid)
    jd = 2449444.5;

    [rsun, rtasc, decl] = sun(jd);

    fprintf(fid,'sun  %16.9f %16.9f %16.9f\n', rsun(1), rsun(2), rsun(3));


end

% --------  moon         - analytical moon ephemeris
function testmoon(fid)
    jd = 2449470.5;
    [rmoon, rtasc, decl] = moon(jd);

    fprintf(fid,'moon %15.11f  %15.11f  %15.11f\n', rmoon(1), rmoon(2), rmoon(3));
end


function testkepler(fid)
    dtsec = 42397.344;  % s

    r1 = [ 2.500000 * re; 0.000000; 0.000000 ];
    % assume circular initial orbit for vel calcs
    v1 = [ 0.0; sqrt(mu / r1(1)); 0.0 ];
    fprintf(fid,'kepler %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f  %16.8f\n', r1(1), r1(2), r1(3), v1(1), v1(2), v1(3), dtsec);

    [r2, v2] = kepler(r1, v1, dtsec);

    fprintf(fid,'kepler %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f\n', r2(1), r2(2), r2(3), v2(1), v2(2), v2(3));

    % test multi-rev case
    [p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper ] = rv2coe (r1, v1);
    period = 2.0 * pi * sqrt(mag(r1)^3 / mu);

    [r2, v3] = kepler(r1, v1, dtsec+7.0*period);

    fprintf(fid,'kepler %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f\n', r2(1), r2(2), r2(3), ...
           v2(1), v2(2), v2(3), dtsec+7.0*period);

    % alt tests
    % initial coes with more than one period = 6281.815597 second
    rad = 180.0/pi;
    [ro, vo] = coe2rv (7358.39, 0.0, 28.5/rad, 0.0/rad, 30.0/rad, 0.0/rad, 0.0, 0.0, 0.0);
    fprintf(fid,'input:\n' );
    fprintf(fid,'ro %16.8f %16.8f %16.8f km\n',ro );
    fprintf(fid,'vo %16.8f %16.8f %16.8f km/s\n',vo );

    % convert 40 minutes to seconds
    dtsec = 4000.0*60.0;
    dtsec = 1.291007302335531e+03;
    dtsec = 6281.815597;
    fprintf(fid,'dt %16.8f second\n',dtsec );
    fprintf(fid,'intermediate values:\n' );

    [r1,v1] =  kepler ( ro,vo, dtsec );

    % ansrwer in km and km/s
    fprintf(fid,'output:\n' );
    fprintf(fid,'r1 %16.8f %16.8f %16.8f er\n',r1/re );
    fprintf(fid,'r1 %16.8f %16.8f %16.8f km\n',r1 );
    fprintf(fid,'v1 %16.8f %16.8f %16.8f er/tu\n',v1/velkmps );
    fprintf(fid,'v1 %16.8f %16.8f %16.8f km/s\n',v1 );

    ro=[-3244.01178958993; 5561.5015207476; 3181.63137126354];
    vo=[-0.311911476329513; 3.55766787343696; -6.53796978233233];
    dtsec = 240.0;
    fprintf(fid,'dt %16.8f second\n',dtsec );
    fprintf(fid,'intermediate values:\n' );

    [r1,v1] =  kepler ( ro,vo, dtsec );

    % ansrwer in km and km/s
    fprintf(fid,'output:\n' );
    fprintf(fid,'r1 %16.8f %16.8f %16.8f er\n',r1/re );
    fprintf(fid,'r1 %16.8f %16.8f %16.8f km\n',r1 );
    fprintf(fid,'v1 %16.8f %16.8f %16.8f er/tu\n',v1/velkmps );
    fprintf(fid,'v1 %16.8f %16.8f %16.8f km/s\n',v1 );
end

function testsunmoonjpl(fid)
    [jd, jdF] = jday(2017, 5, 11, 3, 51, 42.7657);

    fprintf(fid,' =============================   test sun and moon ephemerides =============================\n');

    % read in jpl sun moon files
    % ansrwers
    fprintf(fid,'2017  5 11  0  96576094.2145 106598001.2476 46210616.7776     151081093.9419  0.9804616 -252296.5509 -302841.7334 -93212.7720\n');
    fprintf(fid,'2017  5 11 12  95604355.9737 107353047.2919 46537942.1006     151098145.9151  0.9802403 -218443.5158 -325897.7785 -102799.8515\n');
    fprintf(fid,'2017  5 12  0  94625783.6875 108100430.4112 46861940.2387     151115133.0492  0.9800199 -182165.5046 -345316.4032 -111246.7742\n');

    % for 1 day centers, need to adjust the initjpl function
    %initjplde(ref jpldearr, 'D:/Codes/LIBRARY/DataLib/', 'sunmooneph_430t.txt', out jdjpldestart, out jdjpldestartFrac);
    infilename = append('D:\Codes\LIBRARY\DataLib\', 'sunmooneph_430t12.txt');
    [jpldearr] = readjplde(infilename);

    [rsun, rsmag, rmoon, rmmag] = findjpldeparam(jd, 0.0, 'l', jpldearr);
    fprintf(fid,'findjpldeephem 0000 hrs l \n %16.8f %16.8f  %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f\n', jd, 0.0, rsun(1), ...
        rsun(2), rsun(3), rmoon(1), rmoon(2), rmoon(3));

    [rsun, rsmag, rmoon, rmmag] = findjpldeparam(jd, 0.0, 's', jpldearr);
    fprintf(fid,'findjpldeephem 0000 hrs s \n %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f\n', jd, 0.0,rsun(1), ...
        rsun(2), rsun(3), rmoon(1), rmoon(2), rmoon(3));

    [rsun, rtascs, decls, rmoon, rtascm, declm] = sunmoonjpl(jd, jdF, 's', jpldearr);
    fprintf(fid,'sunmoon 0000 hrs s\n %16.8f %16.8f  %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f\n',...
        jd, jdF, rsun(1), rsun(2), rsun(3),rmoon(1), rmoon(2), rmoon(3));

    [rsun, rsmag, rmoon, rmmag] = findjpldeparam(jd, jdF, 'l', jpldearr);
    fprintf(fid,'findjpldeephem hrs l\n %16.8f %16.8f  %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f\n',...
        jd, jdF, rsun(1), rsun(2), rsun(3),rmoon(1), rmoon(2), rmoon(3));

    [rsun, rtascs, decls, rmoon, rtascm, declm] = sunmoonjpl(jd, jdF, 'l', jpldearr);
    fprintf(fid,'sunmoon hrs l\n %16.8f %16.8f  %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f\n',...
        jd, jdF, rsun(1), rsun(2), rsun(3),rmoon(1), rmoon(2), rmoon(3));

    [rsun, rsmag, rmoon, rmmag] = findjpldeparam(jd, 1.0, 'l', jpldearr);
    fprintf(fid,'findjpldeephem 2400 hrs s\n %16.8f %16.8f  %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f\n',...
        jd, jdF, rsun(1), rsun(2), rsun(3),rmoon(1), rmoon(2), rmoon(3));


    % ex 8.5 test
    [jd, jdF] = jday(2020, 2, 18, 15, 8, 47.23847);
    [rsun, rsmag, rmoon, rmmag] = findjpldeparam(jd, jdF, 's', jpldearr);
    fprintf(fid,'findjpldeephem 0000 hrs s \n%16.8f  %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f\n', jd, rsun(1), ...
        rsun(2), rsun(3), rmoon(1), rmoon(2), rmoon(3));

    % test interpolation of vectors
    % shows spline is MUCH better - 3 km sun variation in mid day linear, 60m diff with spline.
    [jd, jdF] = jday(2017, 5, 11, 3, 51, 42.7657);
    [jd, jdF] = jday(2000, 1, 1, 0, 0, 0.0);
    fprintf(fid,'findjplde  mfme     rsun x             y                 z             rmoon x             y                z      (km)\n');

    % the code that you want to measure comes here
    % read in jpl sun moon files - seems to be the slowest part (800 msec)
    infilename = append('D:\Codes\LIBRARY\DataLib\', 'sunmooneph_430t12.txt');
    [jpldearr] = readjplde(infilename);

    % for ii = 0: 36500
    %     % seems pretty fast (45 msec)
    %     for jj = 0:24
    % 
    % [rsun, rsmag, rmoon, rmmag] = findjpldeparam(jd + ii, (jj * 1.0) / 24.0, 's', jpldearr);
    %         % the write takes some time (160 msec)
    % fprintf(fid,'sunmoon 0000 hrs s\n %16.8f %16.8f  %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f\n',...
    %     jd,  (ii * 60.0), rsun(1), rsun(2), rsun(3),rmoon(1), rmoon(2), rmoon(3));
    %     end
    % end

end


function testkp2ap(fid)

    for i = 1: 27
        kp = 10.0 * i / 3.0;
        [ap] = kp2ap(kp);
        % get spacing correct, leading 0, front spaces
        fprintf(fid,'%15.11f %15.11f  %15.11f\n', i, 0.1 * kp, ap);
    end
end


function testazel2radec(fid)
    rad = 180.0 / pi;

    fileLoc = 'D:\Codes\LIBRARY\DataLib\';
    [iau80arr] = iau80in(fileLoc);

    xp = 0.0;
    yp = 0.0;
    lod = 0.0;
    ddpsi = -0.052195;
    ddeps = -0.003875;
    dut1 = -0.37816;

    year = 2015;
    mon = 12;
    day = 15;
    hr = 16;
    dat = 36;
    minute = 58;
    second = 50.208;
    eqeterms = 2;
    [jd, jdFrac] = jday(year, mon, day, hr, minute, second);

    % note you have to use tdb for time of ineterst AND j2000 (when dat = 32)
    jdtt = jd;
    jdftt = jdFrac + (dat + 32.184) / 86400.0;
    ttt = (jdtt + jdftt - 2451545.0) / 36525.0;
    jdut1 = jd + jdFrac + dut1 / 86400.0;

    recef = [ -605.79221660; -5870.22951108; 3493.05319896 ];
    recef = [ -100605.79221660; -1005870.22951108; 1003493.05319896 ];
    vecef = [ -1.56825429; -3.70234891; -6.47948395 ];
    aecef = [0.0; 0.0; 0.0];
    asecef = [0.0; 0.0; 0.0];

    [reci, veci, aeci] = ecef2eci(recef, vecef, aecef, iau80arr, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps );
    lon = -104.883 / rad;
    latgd = 39.007 / rad;
    alt = 2.102;
    [rsecef, vsecef] = site(latgd, lon, alt);

    [rseci, vseci, aseci] = ecef2eci(rsecef', vsecef', asecef', iau80arr, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps );

    [lst, gst] = lstime(lon, jdut1);

    % print out initial conditions
    fprintf(fid,'recef  %15.11f %15.11f %15.11f %15.11f %15.11f %15.11f\n',  recef(1), recef(2), recef(3), vecef(1), vecef(2), vecef(3));
    fprintf(fid,'rs ecef  %15.11f %15.11f %15.11f %15.11f %15.11f %15.11f\n',  rsecef(1), rsecef(2), rsecef(3), vsecef(1), vsecef(2), vsecef(3));
    fprintf(fid,'reci  %15.11f %15.11f %15.11f %15.11f %15.11f %15.11f\n',  reci(1), reci(2), reci(3), veci(1), veci(2), veci(3));
    fprintf(fid,'rs eci  %15.11f %15.11f %15.11f\n',  rseci(1), rseci(2), rseci(3));

    [rho, az, el, drho, daz, del] = rv_razel(recef, vecef, latgd, lon, alt );

    [rr, rtasc, decl, drr, drtasc, ddecl ] = rv_radec(reci, veci );

    [rho, trtasc, tdecl, drho, dtrtasc, dtdecl] = rv_tradec(reci, veci, rseci);

    % prout results
    fprintf(fid,'razel  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', rho, az, el, drho, daz, del);
    fprintf(fid,'radec  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', rho, rtasc, decl, drho, drtasc, ddecl);
    fprintf(fid,'tradec  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', rho, trtasc, tdecl, drr, dtrtasc, dtdecl);

    [rtasc, decl, rtasc1] = azel_radec(az, el, lst, latgd);
    fprintf(fid,'radec ' + rtasc + ' rtasc1 ' + rtasc1, decl);
    fprintf(fid,'radec ' + (pi * 2 - rtasc) + ' rtasc1 ' + (pi * 2 - rtasc1), decl);
end


% ------------------------------------------------------------------------------
%
%                           function LegPolyEx
%
%   this function finds the exact (from equations) Legendre polynomials for the gravity field.
%   note that the arrays are indexed from 0 to coincide with the usual nomenclature (eq 8-21
%   in my text). fortran implementations will have indicies of 1 greater as they often
%   start at 1. these are exact expressions derived from mathematica.
%
%  author        : david vallado             davallado@gmail.com  16 dec 2019
%
%  inputs        description                                   range / units
%    latgc       - Geocentric lat of satellite                   pi to pi rad
%    order       - size of gravity field                         1..2160..
%
%  outputs       :
%    LegArr      - [,] array of Legendre polynomials
%
%  locals :
%    L,m         - degree and order indices
%    conv        - conversion to un-unitalize
%
%  coupling      :
%   none
%
%  references :
%    vallado       2013, 597, Eq 8-57
% [LegArrEx] = LegPolyEx(latgc, order)
% ------------------------------------------------------------------------------*/

% function [LegArrEx] = LegPolyEx(latgc, order)
% 
%     % -------------------- exact epxressions ----------------------
%     LegArrEx(3, 1) = 0.5 * (3 * s * s - 1.0);
%     LegArrEx(3, 2) = 3.0 * s * c;
%     LegArrEx(3, 3) = 3.0 * c * c;
% 
%     % include (-1)^m for all the terms
%     LegArrEx(4, 1) = -0.5 * s * (3 - 5 * s * s);
%     LegArrEx(4, 2) = (3.0 / 2) * c * (-1 + 5 * s * s); % 15s^2 - 3
%     LegArrEx(4, 3) = 15 * s * c * c;
%     LegArrEx(4, 4) = 15 * (c * c * c);
% 
%     LegArrEx(5, 1) = 1.0 / 8.0 * (35.0 * s^4) - 30.0 * s^2) + 3.0);
%     LegArrEx(5, 2) = 2.5 * c * (-3 * s + 7 * s^3));
%     LegArrEx(5, 3) = 7.5 * c^2) * (-1 + 7 * s^2));
%     LegArrEx(5, 4) = 105.0 * c^3) * s;
%     LegArrEx(5, 5) = 105.0 * c^4);
% 
%     LegArrEx(6, 1) = (1.0 / 8) * s * (15 - 70 * s^2) + 63 * s^4));
%     LegArrEx(6, 2) = (15.0 / 8) * c * (1 - 14 * s^2) + 21 * s^4));
%     LegArrEx(6, 3) = (105.0 / 2) * c * c * (-s + 3 * s^3));
%     LegArrEx(6, 4) = (105.0 / 2) * c^3) * (-1 + 9 * s^2));
%     LegArrEx(6, 5) = 945.0 * s * c^4);
%     LegArrEx(6, 6) = 945.0 * c^5);
% 
%     LegArrEx[6, 0] = 1.0 / 16 * (-5 + 105 * s^2) - 315 * s^4) + 231 * s^6));
%     LegArrEx[6, 1] = (21.0 / 8) * c * (5 * s - 30 * s^3) + 33 * s^5));
%     LegArrEx[6, 2] = (105.0 / 8) * c * c * (1 - 18 * s^2) + 33 * s^4));
%     LegArrEx[6, 3] = (315.0 / 2) * c^3) * (-3 * s + 11 * s^3));
%     LegArrEx[6, 4] = 945.0 / 2 * Math.Pow(c * c, 2) * (-1 + 11 * s^2));
%     LegArrEx[6, 5] = 10395.0 * s * c^5);
%     LegArrEx[6, 6] = 10395.0 * Math.Pow(c * c, 3);
% 
%     LegArrEx[7, 0] = 1.0 / 16 * (-35 * s + 315 * s^3) - 693 * s^5) + 429 * s^7));
%     LegArrEx[7, 1] = (7.0 / 16) * c * (-5 + 135 * s^2) - 495 * s^4) + 429 * s^6));
%     LegArrEx[7, 2] = (63.0 / 8) * c * c * (15 * s - 110 * s^3) + 143 * s^5));
%     LegArrEx[7, 3] = (315.0 / 8) * c^3) * (3 - 66 * s^2) + 143 * s^4));
%     LegArrEx[7, 4] = 3465.0 / 2 * Math.Pow(c * c, 2) * (-3 * s + 13 * s^3));
%     LegArrEx[7, 5] = (10395.0 / 2) * c^5) * (-1 + 13 * s^2));
%     LegArrEx[7, 6] = 135135.0 * s * Math.Pow(c * c, 3);
%     LegArrEx[7, 7] = 135135.0 * c^7);
% 
%     LegArrEx[8, 0] = 1.0 / 128 * (35 - 1260 * s^2) + 6930 * s^4) - 12012 * s^6) + 6435 * s^8));
%     LegArrEx[8, 1] = (9.0 / 16) * c * (-35 * s + 385 * s^3) - 1001 * s^5) + 715 * s^7));
%     LegArrEx[8, 2] = (315.0 / 16) * c * c * (-1 + 33 * s^2) - 143 * s^4) + 143 * s^6));
%     LegArrEx[8, 3] = (3465.0 / 8) * c^3) * (3 * s - 26 * s^3) + 39 * s^5));
%     LegArrEx[8, 4] = 10395.0 / 8 * Math.Pow(c * c, 2) * (1 - 26 * s^2) + 65 * s^4));
%     LegArrEx[8, 5] = (135135.0 / 2) * c^5) * (-s + 5 * s^3));
%     LegArrEx[8, 6] = (135135.0 / 2) * Math.Pow(c * c, 3) * (-1 + 15 * s^2));
%     LegArrEx[8, 7] = 2027025.0 * s * c^7);
%     LegArrEx[8, 8] = 2027025.0 * Math.Pow(c * c, 4);
% 
%     LegArrEx[9, 0] = 1.0 / 128 * (315 * s - 4620 * s^3) + 18018 * s^5) - 25740 * s^7) + 12155 * s^9));
%     LegArrEx[9, 1] = (45.0 / 128) * c * (7 - 308 * s^2) + 2002 * s^4) - 4004 * s^6) + 2431 * s^8));
%     LegArrEx[9, 2] = (495.0 / 16) * c * c * (-7 * s + 91 * s^3) - 273 * s^5) + 221 * s^7));
%     LegArrEx[9, 3] = (3465.0 / 16) * c^3) * (-1 + 39 * s^2) - 195 * s^4) + 221 * s^6));
%     LegArrEx[9, 4] = 135135.0 / 8 * Math.Pow(c * c, 2) * (s - 10 * s^3) + 17 * s^5));
%     LegArrEx[9, 5] = (135135.0 / 8) * c^5) * (1 - 30 * s^2) + 85 * s^4));
%     LegArrEx[9, 6] = (675675.0 / 2) * Math.Pow(c * c, 3) * (-3 * s + 17 * s^3));
%     LegArrEx[9, 7] = (2027025.0 / 2) * c^7) * (-1 + 17 * s^2));
%     LegArrEx[9, 8] = 34459425.0 * s * Math.Pow(c * c, 4);
%     LegArrEx[9, 9] = 34459425.0 * c^9);
% 
%     LegArrEx[10, 0] = 1.0 / 256 * (-63 + 3465 * s^2) - 30030 * s^4) + 90090 * s^6) - 109395 * s^8) + 46189 * s^10));
%     LegArrEx[10, 1] = (55.0 / 128) * c * (63 * s - 1092 * s^3) + 4914 * s^5) - 7956 * s^7) + 4199 * s^9));
%     LegArrEx[10, 2] = (495.0 / 128) * c * c * (7 - 364 * s^2) + 2730 * s^4) - 6188 * s^6) + 4199 * s^8));
%     LegArrEx[10, 3] = (6435.0 / 16) * c^3) * (-7 * s + 105 * s^3) - 357 * s^5) + 323 * s^7));
%     LegArrEx[10, 4] = 45045.0 / 16 * Math.Pow(c * c, 2) * (-1 + 45 * s^2) - 255 * s^4) + 323 * s^6));
%     LegArrEx[10, 5] = (135135.0 / 8) * c^5) * (15 * s - 170 * s^3) + 323 * s^5));
%     LegArrEx[10, 6] = (675675.0 / 8) * Math.Pow(c * c, 3) * (3 - 102 * s^2) + 323 * s^4));
%     LegArrEx[10, 7] = (11486475.0 / 2) * c^7) * (-3 * s + 19 * s^3));
%     LegArrEx[10, 8] = 34459425.0 / 2 * Math.Pow(c * c, 4) * (-1 + 19 * s^2));
%     LegArrEx[10, 9] = 654729075.0 * s * c^9);
%     LegArrEx[10, 10] = 654729075.0 * Math.Pow(c * c, 5);
% 
%     LegArrEx[11, 0] = 1.0 / 256 * (-693 * s + 15015 * s^3) - 90090 * s^5) + 218790 * s^7) - 230945 * s^9) + 88179 * s^11));
%     LegArrEx[11, 1] = (33.0 / 256) * c * (-21 + 1365 * s^2) - 13650 * s^4) + 46410 * s^6) - 62985 * s^8) + 29393 * s^10));
%     LegArrEx[11, 2] = (2145.0 / 128) * c * c * (21 * s - 420 * s^3) + 2142 * s^5) - 3876 * s^7) + 2261 * s^9));
%     LegArrEx[11, 3] = (45045.0 / 128) * c^3) * (1 - 60 * s^2) + 510 * s^4) - 1292 * s^6) + 969 * s^8));
%     LegArrEx[11, 4] = 135135.0 / 16 * Math.Pow(c * c, 2) * (-5 * s + 85 * s^3) - 323 * s^5) + 323 * s^7));
%     LegArrEx[11, 5] = (135135.0 / 16) * c^5) * (-5 + 255 * s^2) - 1615 * s^4) + 2261 * s^6));
%     LegArrEx[11, 6] = (2297295.0 / 8) * Math.Pow(c * c, 3) * (15 * s - 190 * s^3) + 399 * s^5));
%     LegArrEx[11, 7] = (34459425.0 / 8) * c^7) * (1 - 38 * s^2) + 133 * s^4));
%     LegArrEx[11, 8] = 654729075.0 / 2 * Math.Pow(c * c, 4) * (-s + 7 * s^3));
%     LegArrEx[11, 9] = (654729075.0 / 2) * c^9) * (-1 + 21 * s^2));
%     LegArrEx[11, 10] = 13749310575.0 * s * Math.Pow(c * c, 5);
%     LegArrEx[11, 11] = 13749310575.0 * c^11);
% 
%     LegArrEx[12, 0] = (231.0 - 18018 * s^2) + 225225 * s^4) - 1021020 * s^6) + 2078505 * s^8) - 1939938 * s^10) + 676039 * s^12)) / 1024;
%     LegArrEx[12, 1] = (39.0 / 256) * c * (-231 * s + 5775 * s^3) - 39270 * s^5) + 106590 * s^7) - 124355 * s^9) + 52003 * s^11));
%     LegArrEx[12, 2] = (3003.0 / 256) * c * c * (-3 + 225 * s^2) - 2550 * s^4) + 9690 * s^6) - 14535 * s^8) + 7429 * s^10));
%     LegArrEx[12, 3] = (15015.0 / 128) * c^3) * (45 * s - 1020 * s^3) + 5814 * s^5) - 11628 * s^7) + 7429 * s^9));
%     LegArrEx[12, 4] = 135135.0 / 128 * Math.Pow(c * c, 2) * (5 - 340 * s^2) + 3230 * s^4) - 9044 * s^6) + 7429 * s^8));
%     LegArrEx[12, 5] = (2297295.0 / 16) * c^5) * (-5 * s + 95 * s^3) - 399 * s^5) + 437 * s^7));
%     LegArrEx[12, 6] = (2297295.0 / 16) * Math.Pow(c * c, 3) * (-5 + 285 * s^2) - 1995 * s^4) + 3059 * s^6));
%     LegArrEx[12, 7] = (130945815.0 / 8) * c^7) * (5 * s - 70 * s^3) + 161 * s^5));
%     LegArrEx[12, 8] = 654729075.0 / 8 * Math.Pow(c * c, 4) * (1 - 42 * s^2) + 161 * s^4));
%     LegArrEx[12, 9] = (4583103525.0 / 2) * c^9) * (-3 * s + 23 * s^3));
%     LegArrEx[12, 10] = (13749310575.0 / 2) * Math.Pow(c * c, 5) * (-1 + 23 * s^2));
%     LegArrEx[12, 11] = 316234143225.0 * s * c^11);
%     LegArrEx[12, 12] = 316234143225.0 * Math.Pow(c * c, 6);
% 
%     LegArrEx[13, 0] = (3003.0 * s - 90090 * s^3) + 765765 * s^5) - 2771340 * s^7) + 4849845 * s^9) - 4056234 * s^11) + 1300075 * s^13)) / 1024;
%     LegArrEx[13, 1] = ((91.0 * c * (33 - 2970 * s^2) + 42075 * s^4) - 213180 * s^6) + 479655 * s^8) - 490314 * s^10) + 185725 * s^12)) / 1024));
%     LegArrEx[13, 2] = (1365.0 / 256) * c * c * (-99 * s + 2805 * s^3) - 21318 * s^5) + 63954 * s^7) - 81719 * s^9) + 37145 * s^11));
%     LegArrEx[13, 3] = (15015.0 / 256) * c^3) * (-9 + 765 * s^2) - 9690 * s^4) + 40698 * s^6) - 66861 * s^8) + 37145 * s^10));
%     LegArrEx[13, 4] = 255255.0 / 128 * Math.Pow(c * c, 2) * (45 * s - 1140 * s^3) + 7182 * s^5) - 15732 * s^7) + 10925 * s^9));
%     LegArrEx[13, 5] = (2297295.0 / 128) * c^5) * (5 - 380 * s^2) + 3990 * s^4) - 12236 * s^6) + 10925 * s^8));
%     LegArrEx[13, 6] = (43648605.0 / 16) * Math.Pow(c * c, 3) * (-5 * s + 105 * s^3) - 483 * s^5) + 575 * s^7));
%     LegArrEx[13, 7] = (218243025.0 / 16) * c^7) * (-1 + 63 * s^2) - 483 * s^4) + 805 * s^6));
%     LegArrEx[13, 8] = 4583103525.0 / 8 * Math.Pow(c * c, 4) * (3 * s - 46 * s^3) + 115 * s^5));
%     LegArrEx[13, 9] = (4583103525.0 / 8) * c^9) * (3 - 138 * s^2) + 575 * s^4));
%     LegArrEx[13, 10] = (105411381075.0 / 2) * Math.Pow(c * c, 5) * (-3 * s + 25 * s^3));
%     LegArrEx[13, 11] = (316234143225.0 / 2) * c^11) * (-1 + 25 * s^2));
%     LegArrEx[13, 12] = 7905853580625.0 * s * Math.Pow(c * c, 6);
%     LegArrEx[13, 13] = 7905853580625.0 * c^13);
% 
%     LegArrEx[14, 0] = (-429.0 + 45045 * s^2) - 765765 * s^4) + 4849845 * s^6) - 14549535 * s^8) + 22309287 * s^10) - 16900975 * s^12) + 5014575 * s^14)) / 2048;
%     LegArrEx[14, 1] = ((105.0 * c * (429 * s - 14586 * s^3) + 138567 * s^5) - 554268 * s^7) + 1062347 * s^9) - 965770 * s^11) + 334305 * s^13)) / 1024));
%     LegArrEx[14, 2] = ((1365.0 * c * c * (33 - 3366 * s^2) + 53295 * s^4) - 298452 * s^6) + 735471 * s^8) - 817190 * s^10) + 334305 * s^12)) / 1024));
%     LegArrEx[14, 3] = (23205.0 / 256) * c^3) * (-99 * s + 3135 * s^3) - 26334 * s^5) + 86526 * s^7) - 120175 * s^9) + 58995 * s^11));
%     LegArrEx[14, 4] = 2297295.0 / 256 * Math.Pow(c * c, 2) * (-1 + 95 * s^2) - 1330 * s^4) + 6118 * s^6) - 10925 * s^8) + 6555 * s^10));
%     LegArrEx[14, 5] = (43648605.0 / 128) * c^5) * (5 * s - 140 * s^3) + 966 * s^5) - 2300 * s^7) + 1725 * s^9));
%     LegArrEx[14, 6] = (218243025.0 / 128) * Math.Pow(c * c, 3) * (1 - 84 * s^2) + 966 * s^4) - 3220 * s^6) + 3105 * s^8));
%     LegArrEx[14, 7] = (654729075.0 / 16) * c^7) * (-7 * s + 161 * s^3) - 805 * s^5) + 1035 * s^7));
%     LegArrEx[14, 8] = 4583103525.0 / 16 * Math.Pow(c * c, 4) * (-1 + 69 * s^2) - 575 * s^4) + 1035 * s^6));
%     LegArrEx[14, 9] = (105411381075.0 / 8) * c^9) * (3 * s - 50 * s^3) + 135 * s^5));
%     LegArrEx[14, 10] = (316234143225.0 / 8) * Math.Pow(c * c, 5) * (1 - 50 * s^2) + 225 * s^4));
%     LegArrEx[14, 11] = (7905853580625.0 / 2) * c^11) * (-s + 9 * s^3));
%     LegArrEx[14, 12] = 7905853580625.0 / 2 * Math.Pow(c * c, 6) * (-1 + 27 * s^2));
%     LegArrEx[14, 13] = 213458046676875.0 * s * c^13);
%     LegArrEx[14, 14] = 213458046676875.0 * Math.Pow(c * c, 7);
% 
%     LegArrEx[15, 0] = (-6435.0 * s + 255255 * s^3) - 2909907 * s^5) + 14549535 * s^7) - 37182145 * s^9) + 50702925 * s^11) - 35102025 * s^13) + 9694845 * s^15)) / 2048;
%     LegArrEx[15, 1] = ((15.0 * c * (-429 + 51051 * s^2) - 969969 * s^4) + 6789783 * s^6) - 22309287 * s^8) + 37182145 * s^10) - 30421755 * s^12) + 9694845 * s^14)) / 2048));
%     LegArrEx[15, 2] = ((1785.0 * c * c * (429 * s - 16302 * s^3) + 171171 * s^5) - 749892 * s^7) + 1562275 * s^9) - 1533870 * s^11) + 570285 * s^13)) / 1024));
%     LegArrEx[15, 3] = ((69615.0 * c^3) * (11 - 1254 * s^2) + 21945 * s^4) - 134596 * s^6) + 360525 * s^8) - 432630 * s^10) + 190095 * s^12)) / 1024));
%     LegArrEx[15, 4] = 3968055.0 / 256 * Math.Pow(c * c, 2) * (-11 * s + 385 * s^3) - 3542 * s^5) + 12650 * s^7) - 18975 * s^9) + 10005 * s^11));
%     LegArrEx[15, 5] = (43648605.0 / 256) * c^5) * (-1 + 105 * s^2) - 1610 * s^4) + 8050 * s^6) - 15525 * s^8) + 10005 * s^10));
%     LegArrEx[15, 6] = (218243025.0 / 128) * Math.Pow(c * c, 3) * (21 * s - 644 * s^3) + 4830 * s^5) - 12420 * s^7) + 10005 * s^9));
%     LegArrEx[15, 7] = (654729075.0 / 128) * c^7) * (7 - 644 * s^2) + 8050 * s^4) - 28980 * s^6) + 30015 * s^8));
%     LegArrEx[15, 8] = 15058768725.0 / 16 * Math.Pow(c * c, 4) * (-7 * s + 175 * s^3) - 945 * s^5) + 1305 * s^7));
%     LegArrEx[15, 9] = (105411381075.0 / 16) * c^9) * (-1 + 75 * s^2) - 675 * s^4) + 1305 * s^6));
%     LegArrEx[15, 10] = (1581170716125.0 / 8) * Math.Pow(c * c, 5) * (5 * s - 90 * s^3) + 261 * s^5));
%     LegArrEx[15, 11] = (7905853580625.0 / 8) * c^11) * (1 - 54 * s^2) + 261 * s^4));
%     LegArrEx[15, 12] = 71152682225625.0 / 2 * Math.Pow(c * c, 6) * (-3 * s + 29 * s^3));
%     LegArrEx[15, 13] = (213458046676875.0 / 2) * c^13) * (-1 + 29 * s^2));
%     LegArrEx[15, 14] = 6190283353629375.0 * s * Math.Pow(c * c, 7);
%     LegArrEx[15, 15] = 6190283353629375.0 * c^15);
% 
%     LegArrEx[16, 0] = (6435.0 - 875160 * s^2) + 19399380 * s^4) - 162954792 * s^6) + 669278610 * s^8) - 1487285800 * s^10) + 1825305300 * s^12) - 1163381400 * s^14) + 300540195 * s^16)) / 32768;
%     LegArrEx[16, 1] = ((17.0 * c * (-6435 * s + 285285 * s^3) - 3594591 * s^5) + 19684665 * s^7) - 54679625 * s^9) + 80528175 * s^11) - 59879925 * s^13) + 17678835 * s^15)) / 2048));
%     LegArrEx[16, 2] = ((765.0 * c * c * (-143 + 19019 * s^2) - 399399 * s^4) + 3062059 * s^6) - 10935925 * s^8) + 19684665 * s^10) - 17298645 * s^12) + 5892945 * s^14)) / 2048));
%     LegArrEx[16, 3] = ((101745.0 * c^3) * (143 * s - 6006 * s^3) + 69069 * s^5) - 328900 * s^7) + 740025 * s^9) - 780390 * s^11) + 310155 * s^13)) / 1024));
%     LegArrEx[16, 4] = (1322685.0 * Math.Pow(c * c, 2) * (11 - 1386 * s^2) + 26565 * s^4) - 177100 * s^6) + 512325 * s^8) - 660330 * s^10) + 310155 * s^12)) / 1024);
%     LegArrEx[16, 5] = (3968055.0 / 256) * c^5) * (-231 * s + 8855 * s^3) - 88550 * s^5) + 341550 * s^7) - 550275 * s^9) + 310155 * s^11));
%     LegArrEx[16, 6] = (43648605.0 / 256) * Math.Pow(c * c, 3) * (-21 + 2415 * s^2) - 40250 * s^4) + 217350 * s^6) - 450225 * s^8) + 310155 * s^10));
%     LegArrEx[16, 7] = (5019589575.0 / 128) * c^7) * (21 * s - 700 * s^3) + 5670 * s^5) - 15660 * s^7) + 13485 * s^9));
%     LegArrEx[16, 8] = 15058768725.0 / 128 * Math.Pow(c * c, 4) * (7 - 700 * s^2) + 9450 * s^4) - 36540 * s^6) + 40455 * s^8));
%     LegArrEx[16, 9] = (75293843625.0 / 16) * c^9) * (-35 * s + 945 * s^3) - 5481 * s^5) + 8091 * s^7));
%     LegArrEx[16, 10] = (527056905375.0 / 16) * Math.Pow(c * c, 5) * (-5 + 405 * s^2) - 3915 * s^4) + 8091 * s^6));
%     LegArrEx[16, 11] = (14230536445125.0 / 8) * c^11) * (15 * s - 290 * s^3) + 899 * s^5));
%     LegArrEx[16, 12] = 71152682225625.0 / 8 * Math.Pow(c * c, 6) * (3 - 174 * s^2) + 899 * s^4));
%     LegArrEx[16, 13] = (2063427784543125.0 / 2) * c^13) * (-3 * s + 31 * s^3));
%     LegArrEx[16, 14] = (6190283353629375.0 / 2) * Math.Pow(c * c, 7) * (-1 + 31 * s^2));
%     LegArrEx[16, 15] = 191898783962510625.0 * s * c^15);
%     LegArrEx[16, 16] = 191898783962510625.0 * Math.Pow(c * c, 8);
% 
%     LegArrEx[17, 0] = (109395.0 * s - 5542680 * s^3) + 81477396 * s^5) - 535422888 * s^7) + 1859107250 * s^9) - 3650610600 * s^11) + 4071834900 * s^13) - 2404321560 * s^15) + 583401555 * s^17)) / 32768;
%     LegArrEx[17, 1] = ((153.0 * c * (715 - 108680 * s^2) + 2662660 * s^4) - 24496472 * s^6) + 109359250 * s^8) - 262462200 * s^10) + 345972900 * s^12) - 235717800 * s^14) + 64822395 * s^16)) / 32768));
%     LegArrEx[17, 2] = ((2907.0 * c * c * (-715 * s + 35035 * s^3) - 483483 * s^5) + 2877875 * s^7) - 8633625 * s^9) + 13656825 * s^11) - 10855425 * s^13) + 3411705 * s^15)) / 2048));
%     LegArrEx[17, 3] = ((14535.0 * c^3) * (-143 + 21021 * s^2) - 483483 * s^4) + 4029025 * s^6) - 15540525 * s^8) + 30045015 * s^10) - 28224105 * s^12) + 10235115 * s^14)) / 2048));
%     LegArrEx[17, 4] = (305235.0 * Math.Pow(c * c, 2) * (1001 * s - 46046 * s^3) + 575575 * s^5) - 2960100 * s^7) + 7153575 * s^9) - 8064030 * s^11) + 3411705 * s^13)) / 1024);
%     LegArrEx[17, 5] = ((43648605.0 * c^5) * (7 - 966 * s^2) + 20125 * s^4) - 144900 * s^6) + 450225 * s^8) - 620310 * s^10) + 310155 * s^12)) / 1024));
%     LegArrEx[17, 6] = (1003917915.0 / 256) * Math.Pow(c * c, 3) * (-21 * s + 875 * s^3) - 9450 * s^5) + 39150 * s^7) - 67425 * s^9) + 40455 * s^11));
%     LegArrEx[17, 7] = (3011753745.0 / 256) * c^7) * (-7 + 875 * s^2) - 15750 * s^4) + 91350 * s^6) - 202275 * s^8) + 148335 * s^10));
%     LegArrEx[17, 8] = 75293843625.0 / 128 * Math.Pow(c * c, 4) * (35 * s - 1260 * s^3) + 10962 * s^5) - 32364 * s^7) + 29667 * s^9));
%     LegArrEx[17, 9] = (75293843625.0 / 128) * c^9) * (35 - 3780 * s^2) + 54810 * s^4) - 226548 * s^6) + 267003 * s^8));
%     LegArrEx[17, 10] = (2032933777875.0 / 16) * Math.Pow(c * c, 5) * (-35 * s + 1015 * s^3) - 6293 * s^5) + 9889 * s^7));
%     LegArrEx[17, 11] = (14230536445125.0 / 16) * c^11) * (-5 + 435 * s^2) - 4495 * s^4) + 9889 * s^6));
%     LegArrEx[17, 12] = 412685556908625.0 / 8 * Math.Pow(c * c, 6) * (15 * s - 310 * s^3) + 1023 * s^5));
%     LegArrEx[17, 13] = (6190283353629375.0 / 8) * c^13) * (1 - 62 * s^2) + 341 * s^4));
%     LegArrEx[17, 14] = (191898783962510625.0 / 2) * Math.Pow(c * c, 7) * (-s + 11 * s^3));
%     LegArrEx[17, 15] = (191898783962510625.0 / 2) * c^15) * (-1 + 33 * s^2));
%     LegArrEx[17, 16] = 6332659870762850625.0 * s * Math.Pow(c * c, 8);
%     LegArrEx[17, 17] = 6332659870762850625.0 * c^17);
% 
%     LegArrEx[18, 0] = (-12155.0 + 2078505 * s^2) - 58198140 * s^4) + 624660036 * s^6) - 3346393050 * s^8) + 10039179150 * s^10) - 17644617900 * s^12) + 18032411700 * s^14) - 9917826435 * s^16) + 2268783825 * s^18)) / 65536;
%     LegArrEx[18, 1] = ((171.0 * c * (12155 * s - 680680 * s^3) + 10958948 * s^5) - 78278200 * s^7) + 293543250 * s^9) - 619109400 * s^11) + 738168900 * s^13) - 463991880 * s^15) + 119409675 * s^17)) / 32768));
%     LegArrEx[18, 2] = ((14535.0 * c * c * (143 - 24024 * s^2) + 644644 * s^4) - 6446440 * s^6) + 31081050 * s^8) - 80120040 * s^10) + 112896420 * s^12) - 81880920 * s^14) + 23881935 * s^16)) / 32768));
%     LegArrEx[18, 3] = ((101745.0 * c^3) * (-429 * s + 23023 * s^3) - 345345 * s^5) + 2220075 * s^7) - 7153575 * s^9) + 12096045 * s^11) - 10235115 * s^13) + 3411705 * s^15)) / 2048));
%     LegArrEx[18, 4] = (3357585.0 * Math.Pow(c * c, 2) * (-13 + 2093 * s^2) - 52325 * s^4) + 470925 * s^6) - 1950975 * s^8) + 4032015 * s^10) - 4032015 * s^12) + 1550775 * s^14)) / 2048);
%     LegArrEx[18, 5] = ((77224455.0 * c^5) * (91 * s - 4550 * s^3) + 61425 * s^5) - 339300 * s^7) + 876525 * s^9) - 1051830 * s^11) + 471975 * s^13)) / 1024));
%     LegArrEx[18, 6] = ((1003917915.0 * Math.Pow(c * c, 3) * (7 - 1050 * s^2) + 23625 * s^4) - 182700 * s^6) + 606825 * s^8) - 890010 * s^10) + 471975 * s^12)) / 1024));
%     LegArrEx[18, 7] = (75293843625.0 / 256) * c^7) * (-7 * s + 315 * s^3) - 3654 * s^5) + 16182 * s^7) - 29667 * s^9) + 18879 * s^11));
%     LegArrEx[18, 8] = 75293843625.0 / 256 * Math.Pow(c * c, 4) * (-7 + 945 * s^2) - 18270 * s^4) + 113274 * s^6) - 267003 * s^8) + 207669 * s^10));
%     LegArrEx[18, 9] = (225881530875.0 / 128) * c^9) * (315 * s - 12180 * s^3) + 113274 * s^5) - 356004 * s^7) + 346115 * s^9));
%     LegArrEx[18, 10] = (14230536445125.0 / 128) * Math.Pow(c * c, 5) * (5 - 580 * s^2) + 8990 * s^4) - 39556 * s^6) + 49445 * s^8));
%     LegArrEx[18, 11] = (412685556908625.0 / 16) * c^11) * (-5 * s + 155 * s^3) - 1023 * s^5) + 1705 * s^7));
%     LegArrEx[18, 12] = 2063427784543125.0 / 16 * Math.Pow(c * c, 6) * (-1 + 93 * s^2) - 1023 * s^4) + 2387 * s^6));
%     LegArrEx[18, 13] = (191898783962510625.0 / 8) * c^13) * (s - 22 * s^3) + 77 * s^5));
%     LegArrEx[18, 14] = (191898783962510625.0 / 8) * Math.Pow(c * c, 7) * (1 - 66 * s^2) + 385 * s^4));
%     LegArrEx[18, 15] = (2110886623587616875.0 / 2) * c^15) * (-3 * s + 35 * s^3));
%     LegArrEx[18, 16] = 6332659870762850625.0 / 2 * Math.Pow(c * c, 8) * (-1 + 35 * s^2));
%     LegArrEx[18, 17] = 221643095476699771875.0 * s * c^17);
%     LegArrEx[18, 18] = 221643095476699771875.0 * Math.Pow(c * c, 9);
% 
%     LegArrEx[19, 0] = (-230945.0 * s + 14549535 * s^3) - 267711444 * s^5) + 2230928700 * s^7) - 10039179150 * s^9) + 26466926850 * s^11) - 42075627300 * s^13) + 39671305740 * s^15) - 20419054425 * s^17) + 4418157975 * s^19)) / 65536;
%     LegArrEx[19, 1] = ((95.0 * c * (-2431 + 459459 * s^2) - 14090076 * s^4) + 164384220 * s^6) - 951080130 * s^8) + 3064591530 * s^10) - 5757717420 * s^12) + 6263890380 * s^14) - 3653936055 * s^16) + 883631595 * s^18)) / 65536));
%     LegArrEx[19, 2] = ((5985.0 * c * c * (7293 * s - 447304 * s^3) + 7827820 * s^5) - 60386040 * s^7) + 243221550 * s^9) - 548354040 * s^11) + 695987820 * s^13) - 463991880 * s^15) + 126233085 * s^17)) / 32768));
%     LegArrEx[19, 3] = ((1119195.0 * c^3) * (39 - 7176 * s^2) + 209300 * s^4) - 2260440 * s^6) + 11705850 * s^8) - 32256120 * s^10) + 48384180 * s^12) - 37218600 * s^14) + 11475735 * s^16)) / 32768));
%     LegArrEx[19, 4] = (25741485.0 * Math.Pow(c * c, 2) * (-39 * s + 2275 * s^3) - 36855 * s^5) + 254475 * s^7) - 876525 * s^9) + 1577745 * s^11) - 1415925 * s^13) + 498945 * s^15)) / 2048);
%     LegArrEx[19, 5] = ((77224455.0 * c^5) * (-13 + 2275 * s^2) - 61425 * s^4) + 593775 * s^6) - 2629575 * s^8) + 5785065 * s^10) - 6135675 * s^12) + 2494725 * s^14)) / 2048));
%     LegArrEx[19, 6] = ((1930611375.0 * Math.Pow(c * c, 3) * (91 * s - 4914 * s^3) + 71253 * s^5) - 420732 * s^7) + 1157013 * s^9) - 1472562 * s^11) + 698523 * s^13)) / 1024));
%     LegArrEx[19, 7] = ((25097947875.0 * c^7) * (7 - 1134 * s^2) + 27405 * s^4) - 226548 * s^6) + 801009 * s^8) - 1246014 * s^10) + 698523 * s^12)) / 1024));
%     LegArrEx[19, 8] = 225881530875.0 / 256 * Math.Pow(c * c, 4) * (-63 * s + 3045 * s^3) - 37758 * s^5) + 178002 * s^7) - 346115 * s^9) + 232841 * s^11));
%     LegArrEx[19, 9] = (1581170716125.0 / 256) * c^9) * (-9 + 1305 * s^2) - 26970 * s^4) + 178002 * s^6) - 445005 * s^8) + 365893 * s^10));
%     LegArrEx[19, 10] = (45853950767625.0 / 128) * Math.Pow(c * c, 5) * (45 * s - 1860 * s^3) + 18414 * s^5) - 61380 * s^7) + 63085 * s^9));
%     LegArrEx[19, 11] = (2063427784543125.0 / 128) * c^11) * (1 - 124 * s^2) + 2046 * s^4) - 9548 * s^6) + 12617 * s^8));
%     LegArrEx[19, 12] = 63966261320836875.0 / 16 * Math.Pow(c * c, 6) * (-s + 33 * s^3) - 231 * s^5) + 407 * s^7));
%     LegArrEx[19, 13] = (63966261320836875.0 / 16) * c^13) * (-1 + 99 * s^2) - 1155 * s^4) + 2849 * s^6));
%     LegArrEx[19, 14] = (2110886623587616875.0 / 8) * Math.Pow(c * c, 7) * (3 * s - 70 * s^3) + 259 * s^5));
%     LegArrEx[19, 15] = (2110886623587616875.0 / 8) * c^15) * (3 - 210 * s^2) + 1295 * s^4));
%     LegArrEx[19, 16] = 73881031825566590625.0 / 2 * Math.Pow(c * c, 8) * (-3 * s + 37 * s^3));
%     LegArrEx[19, 17] = (221643095476699771875.0 / 2) * c^17) * (-1 + 37 * s^2));
%     LegArrEx[19, 18] = 8200794532637891559375.0 * s * Math.Pow(c * c, 9);
%     LegArrEx[19, 19] = 8200794532637891559375.0 * c^19);
% 
%     LegArrEx[20, 0] = (1.0 / 262144) * (46189 - 9699690 * s^2) + 334639305 * s^4) - 4461857400 * s^6) + 30117537450 * s^8) - 116454478140 * s^10) + 273491577450 * s^12) - 396713057400 * s^14) + 347123925225 * s^16) - 167890003050 * s^18) + 34461632205 * s^20));
%     LegArrEx[20, 1] = ((105.0 * c * (-46189 * s + 3187041 * s^3) - 63740820 * s^5) + 573667380 * s^7) - 2772725670 * s^9) + 7814045070 * s^11) - 13223768580 * s^13) + 13223768580 * s^15) - 7195285845 * s^17) + 1641030105 * s^19)) / 65536));
%     LegArrEx[20, 2] = ((21945.0 * c * c * (-221 + 45747 * s^2) - 1524900 * s^4) + 19213740 * s^6) - 119399670 * s^8) + 411265530 * s^10) - 822531060 * s^12) + 949074300 * s^14) - 585262485 * s^16) + 149184555 * s^18)) / 65536));
%     LegArrEx[20, 3] = ((1514205.0 * c^3) * (663 * s - 44200 * s^3) + 835380 * s^5) - 6921720 * s^7) + 29801850 * s^9) - 71524440 * s^11) + 96282900 * s^13) - 67856520 * s^15) + 19458855 * s^17)) / 32768));
%     LegArrEx[20, 4] = (77224455.0 * Math.Pow(c * c, 2) * (13 - 2600 * s^2) + 81900 * s^4) - 950040 * s^6) + 5259150 * s^8) - 15426840 * s^10) + 24542700 * s^12) - 19957800 * s^14) + 6486285 * s^16)) / 32768);
%     LegArrEx[20, 5] = ((386122275.0 * c^5) * (-65 * s + 4095 * s^3) - 71253 * s^5) + 525915 * s^7) - 1928355 * s^9) + 3681405 * s^11) - 3492615 * s^13) + 1297257 * s^15)) / 2048));
%     LegArrEx[20, 6] = ((25097947875.0 * Math.Pow(c * c, 3) * (-1 + 189 * s^2) - 5481 * s^4) + 56637 * s^6) - 267003 * s^8) + 623007 * s^10) - 698523 * s^12) + 299367 * s^14)) / 2048));
%     LegArrEx[20, 7] = ((225881530875.0 * c^7) * (21 * s - 1218 * s^3) + 18879 * s^5) - 118668 * s^7) + 346115 * s^9) - 465682 * s^11) + 232841 * s^13)) / 1024));
%     LegArrEx[20, 8] = (1581170716125.0 * Math.Pow(c * c, 4) * (3 - 522 * s^2) + 13485 * s^4) - 118668 * s^6) + 445005 * s^8) - 731786 * s^10) + 432419 * s^12)) / 1024);
%     LegArrEx[20, 9] = (45853950767625.0 / 256) * c^9) * (-9 * s + 465 * s^3) - 6138 * s^5) + 30690 * s^7) - 63085 * s^9) + 44733 * s^11));
%     LegArrEx[20, 10] = (137561852302875.0 / 256) * Math.Pow(c * c, 5) * (-3 + 465 * s^2) - 10230 * s^4) + 71610 * s^6) - 189255 * s^8) + 164021 * s^10));
%     LegArrEx[20, 11] = (21322087106945625.0 / 128) * c^11) * (3 * s - 132 * s^3) + 1386 * s^5) - 4884 * s^7) + 5291 * s^9));
%     LegArrEx[20, 12] = 63966261320836875.0 / 128 * Math.Pow(c * c, 6) * (1 - 132 * s^2) + 2310 * s^4) - 11396 * s^6) + 15873 * s^8));
%     LegArrEx[20, 13] = (2110886623587616875.0 / 16) * c^13) * (-s + 35 * s^3) - 259 * s^5) + 481 * s^7));
%     LegArrEx[20, 14] = (2110886623587616875.0 / 16) * Math.Pow(c * c, 7) * (-1 + 105 * s^2) - 1295 * s^4) + 3367 * s^6));
%     LegArrEx[20, 15] = (14776206365113318125.0 / 8) * c^15) * (15 * s - 370 * s^3) + 1443 * s^5));
%     LegArrEx[20, 16] = 221643095476699771875.0 / 8 * Math.Pow(c * c, 8) * (1 - 74 * s^2) + 481 * s^4));
%     LegArrEx[20, 17] = (8200794532637891559375.0 / 2) * c^17) * (-s + 13 * s^3));
%     LegArrEx[20, 18] = (8200794532637891559375.0 / 2) * Math.Pow(c * c, 9) * (-1 + 39 * s^2));
%     LegArrEx[20, 19] = 319830986772877770815625.0 * s * c^19);
%     LegArrEx[20, 20] = 319830986772877770815625.0 * Math.Pow(c * c, 10);
% 
%     LegArrEx[21, 0] = (1.0 / 262144) * (969969 * s - 74364290 * s^3) + 1673196525 * s^5) - 17210021400 * s^7) + 97045398450 * s^9) - 328189892940 * s^11) + 694247850450 * s^13) - 925663800600 * s^15) + 755505013725 * s^17) - 344616322050 * s^19) + 67282234305 * s^21));
%     LegArrEx[21, 1] = (1.0 / 262144) * 231 * c * (4199 - 965770 * s^2) + 36216375 * s^4) - 521515800 * s^6) + 3780989550 * s^8) - 15628090140 * s^10) + 39070225350 * s^12) - 60108039000 * s^14) + 55599936075 * s^16) - 28345065450 * s^18) + 6116566755 * s^20));
%     LegArrEx[21, 2] = ((26565.0 * c * c * (-4199 * s + 314925 * s^3) - 6802380 * s^5) + 65756340 * s^7) - 339741090 * s^9) + 1019223270 * s^11) - 1829375100 * s^13) + 1933910820 * s^15) - 1109154735 * s^17) + 265937685 * s^19)) / 65536));
%     LegArrEx[21, 3] = ((504735.0 * c^3) * (-221 + 49725 * s^2) - 1790100 * s^4) + 24226020 * s^6) - 160929990 * s^8) + 590076630 * s^10) - 1251677700 * s^12) + 1526771700 * s^14) - 992401605 * s^16) + 265937685 * s^18)) / 65536));
%     LegArrEx[21, 4] = (22713075.0 * Math.Pow(c * c, 2) * (1105 * s - 79560 * s^3) + 1615068 * s^5) - 14304888 * s^7) + 65564070 * s^9) - 166890360 * s^11) + 237497820 * s^13) - 176426952 * s^15) + 53187537 * s^17)) / 32768);
%     LegArrEx[21, 5] = ((5019589575.0 * c^5) * (5 - 1080 * s^2) + 36540 * s^4) - 453096 * s^6) + 2670030 * s^8) - 8306760 * s^10) + 13970460 * s^12) - 11974680 * s^14) + 4091349 * s^16)) / 32768));
%     LegArrEx[21, 6] = ((15058768725.0 * Math.Pow(c * c, 3) * (-45 * s + 3045 * s^3) - 56637 * s^5) + 445005 * s^7) - 1730575 * s^9) + 3492615 * s^11) - 3492615 * s^13) + 1363783 * s^15)) / 2048));
%     LegArrEx[21, 7] = ((225881530875.0 * c^7) * (-3 + 609 * s^2) - 18879 * s^4) + 207669 * s^6) - 1038345 * s^8) + 2561251 * s^10) - 3026933 * s^12) + 1363783 * s^14)) / 2048));
%     LegArrEx[21, 8] = (45853950767625.0 * Math.Pow(c * c, 4) * (3 * s - 186 * s^3) + 3069 * s^5) - 20460 * s^7) + 63085 * s^9) - 89466 * s^11) + 47027 * s^13)) / 1024);
%     LegArrEx[21, 9] = ((45853950767625.0 * c^9) * (3 - 558 * s^2) + 15345 * s^4) - 143220 * s^6) + 567765 * s^8) - 984126 * s^10) + 611351 * s^12)) / 1024));
%     LegArrEx[21, 10] = (4264417421389125.0 / 256) * Math.Pow(c * c, 5) * (-3 * s + 165 * s^3) - 2310 * s^5) + 12210 * s^7) - 26455 * s^9) + 19721 * s^11));
%     LegArrEx[21, 11] = (4264417421389125.0 / 256) * c^11) * (-3 + 495 * s^2) - 11550 * s^4) + 85470 * s^6) - 238095 * s^8) + 216931 * s^10));
%     LegArrEx[21, 12] = 234542958176401875.0 / 128 * Math.Pow(c * c, 6) * (9 * s - 420 * s^3) + 4662 * s^5) - 17316 * s^7) + 19721 * s^9));
%     LegArrEx[21, 13] = (2110886623587616875.0 / 128) * c^13) * (1 - 140 * s^2) + 2590 * s^4) - 13468 * s^6) + 19721 * s^8));
%     LegArrEx[21, 14] = (2110886623587616875.0 / 16) * Math.Pow(c * c, 7) * (-35 * s + 1295 * s^3) - 10101 * s^5) + 19721 * s^7));
%     LegArrEx[21, 15] = (14776206365113318125.0 / 16) * c^15) * (-5 + 555 * s^2) - 7215 * s^4) + 19721 * s^6));
%     LegArrEx[21, 16] = 1640158906527578311875.0 / 8 * Math.Pow(c * c, 8) * (5 * s - 130 * s^3) + 533 * s^5));
%     LegArrEx[21, 17] = (8200794532637891559375.0 / 8) * c^17) * (1 - 78 * s^2) + 533 * s^4));
%     LegArrEx[21, 18] = (106610328924292590271875.0 / 2) * Math.Pow(c * c, 9) * (-3 * s + 41 * s^3));
%     LegArrEx[21, 19] = (319830986772877770815625.0 / 2) * c^19) * (-1 + 41 * s^2));
%     LegArrEx[21, 20] = 13113070457687988603440625.0 * s * Math.Pow(c * c, 10);
%     LegArrEx[21, 21] = 13113070457687988603440625.0 * c^21);
% 
% end % LegPolyEx


% ----------------------------------------------------------------------------
% fukushima method (JG 2018)
%   for very large spherical harmonic expansrions and calcs of unitalized associated
%   Legendre polynomials
%
%   Plm are converted to X-numbers
%   Clm, Slm treated as F-numbers
% -----------------------------------------------------------------------------

% initialize legendre function values
% function [p] = pinit(n, m)
%     if (n == 0)
%         p(1) = 1.0;
%     else if (n == 1)
%             p(1) = 1.7320508075688773;
%     else if (n == 2)
%             if (m == 0)
%                 p(1) = 0.5590169943749474;
%                 p(2) = 1.6770509831248423;
%             end
%     else if (m == 1)
%             p(1) = 0.0;
%             p(2) = 1.9364916731037084;
%     end
%     else if (m == 2)
%             p(1) = 0.9682458365518542;
%             p(2) = -0.9682458365518542;
%     end
%     else if (n == 3)
%             if (m == 0)
%                 p(1) = 0.9921567416492215;
%                 p(2) = 1.6535945694153691;
%             end
%     else if (m == 1)
%             p(1) = 0.4050462936504913;
%             p(2) = 2.0252314682524563;
%     end
%     else if (m == 2)
%             p(1) = 1.2808688457449498;
%             p(2) = -1.2808688457449498;
%     end
%     else if (m == 3)
%             p(1) = 1.5687375497513917;
%             p(2) = -0.5229125165837972;
%     end
%     else if (n == 4)
%             if (m == 0)
%                 p(1) = 0.421875;
%                 p(2) = 0.9375;
%                 p(3) = 1.640625;
%             end
%     else if (m == 1)
%             p(1) = 0.0;
%             p(2) = 0.5929270612815711;
%             p(3) = 2.0752447144854989;
%     end
%     else if (m == 2)
%             p(1) = 0.6288941186718159;
%             p(2) = 0.8385254915624211;
%             p(3) = -1.4674196102342370;
%     end
%     else if (m == 3)
%             p(1) = 0.0;
%             p(2) = 1.5687375497513917;
%             p(3) = -0.7843687748756958;
%     end
%     else if (m == 4)
% 
%             p(1) = 0.8319487194983835;
%             p(2) = -1.1092649593311780;
%             p(3) = 0.2773162398327945;
%     end
%     end
%     end
%     end  % pinit


    % /* ----------------------------------------------------------------------------
    % *      x2f
    % *
    % *      convert x to f number
    % *
    % ------------------------------------------------------------------------------*/
    % function x2f(x, ix)
    % 
    %     IND = 960;
    %     BIG = Math.Pow(2.0, IND);
    %     BIGI = Math.Pow(2.0, -IND);
    % 
    %     if (ix == 0)
    %         x2fv = x;
    %     else if (ix == -1)
    %             x2fv = x * BIGI;
    %     else if (ix == 1)
    %             x2fv = x * BIG;
    %     else if (ix < 0)
    %             x2fv = 0.0;
    %     else if (x < 0)
    %             x2fv = -BIG;
    %     else
    %         x2fv = BIG;
    % 
    %         return x2fv;
    %     end
    %     end
    %     end
    %     end
    %     end
    % end
    % 


    % /* ----------------------------------------------------------------------------
    % *                                  xunit
    % *
    % * uses the 'x' factor approach - value and exponent
    % *
    % ------------------------------------------------------------------------------*/
    % 
    % function xunit ( x, ix )
    %     IND = 960;
    % 
    %      BIG = Math.Pow(2, IND);
    %      BIGI = Math.Pow(2, -IND);
    %      BIGS = Math.Pow(2.0, 480);  % IND / 2
    %      BIGSI = Math.Pow(2.0, -480);  % IND / 2
    % 
    %     w = abs(x);
    %     if (w >= BIGS)
    % 
    %         x = x * BIGI;
    %         ix = ix + 1;
    %     end
    % else if (w < BIGSI)
    % 
    %         x = x * BIG;
    %         ix = ix - 1;
    % end
    % end  % xunit


    % /* ----------------------------------------------------------------------------
    % *                                       xl2sum
    % *
    % * routine to compute the two-term linear sum of X-numbers
    % * with F-number coefficients
    % *
    % ---------------------------------------------------------------------------- */
    % [z, iz] = function xlsum2(f, g, x, y, ix, iy)
    %     IND = 960;
    %     BIGI = Math.Pow(2, -IND);
    % 
    %     id = ix - iy;
    %     if (id == 0)
    % 
    %         z = f * x + g * y;
    %         iz = ix;
    %     end
    % else if (id == 1)
    % 
    %         z = f * x + g * (y * BIGI);
    %         iz = ix;
    % end
    % else if (id == -1)
    % 
    %         z = g * y + f * (x * BIGI);
    %         iz = iy;
    % end
    % else if (id > 1)
    % 
    %         z = f * x;
    %         iz = ix;
    % end
    % else
    % 
    %     z = g * y;
    %     iz = iy;
    % end
    % 
    % xunit(ref z, ref iz);
    % end  % xlsum2


    % /* ----------------------------------------------------------------------------
    % *                                  dpeven
    % *
    % * find Pnnj and Pn,n-1,j when degree n is even and n >= 6. the returned values are
    % * (xp(j), ip(j)) and (xp1(j), ip1(j)),  X-number vectors representing
    % * Pnnj and Pn,n-1,j, respectively. required initial values for Pn-2,n-2,j as
    % * (xpold(j), ipold(j)) are needed.
    % *
    % *  inputs          description                                  range / units
    % *    n           - degree
    % *
    % *
    % *
    % *  outputs       :
    % *    xp1         - value
    % *    ip1         - exponent
    % *
    % *  locals        :
    % *                -
    % *
    % *  coupling      :
    % *    xunit
    % *    xlsum2
    % *
    % *  references    :
    % * Fukushima (2012a)
    % ---------------------------------------------------------------------------- */
    % 
    %     [xp, xp1, ip, ip1] = function dpeven(n, xpold, ipold)
    % 
    %         jx = n / 2;
    %         jxm2 = jx - 2;
    %         jxm1 = jx - 1;
    %         n2 = n * 2;
    %         gamma = sqrt(str2num(n2 + 1) * str2num(n2 - 1) / (str2num(n) * str2num(n - 1))) * 0.125;
    %         gamma2 = gamma * 2.0;
    %         xlsum2(gamma2, xpold(1), -gamma, xpold(2), out xp(1), ipold(1), ipold(2), out ip(1));
    %         xlsum2(-gamma2, xpold(1), gamma2, xpold(2), out xtemp, ipold(1), ipold(2), out itemp);
    %         xlsum2(1.0, xtemp, -gamma, xpold(3), out xp(2), itemp, ipold(3), out ip(2));
    %         j = 2;
    %         while (j <= jxm2)
    % 
    %             jm1 = j - 1;
    %             jp1 = j + 1;
    %             xlsum2(-gamma, xpold[jm1], gamma2, xpold(j), out xtemp, ipold[jm1], ipold(j), out itemp);
    %             xlsum2(1.0, xtemp, -gamma, xpold[jp1], out xp(j), itemp, ipold[jp1], out ip(j));
    %             j = j + 1;
    %         end
    %         xlsum2(-gamma, xpold[jxm2], gamma2, xpold[jxm1], out xp[jxm1], ipold[jxm2], ipold[jxm1], out ip[jxm1]);
    %         xp[jx] = -gamma * xpold[jxm1];
    %         ip[jx] = ipold[jxm1];
    %         xunit(ref xp[jx], ref ip[jx]);
    %         alpha2 = sqrt(2.0 / str2num(n)) * 2.0;
    %         xp1(1) = 0.0;
    %         ip1(1) = 0;
    %         j = 1;
    %         while (j <= jx)
    % 
    %             xp1(j) = -str2num(j) * alpha2 * xp(j);
    %             ip1(j) = ip(j);
    %             xunit(ref xp1(j), ref ip1(j));
    %             j = j + 1;
    %         end
    %     end   % dpeven
    % 
    % 
    %     % ----------------------------------------------------------------------------
    %     %                                 dpodd
    %    %
    %        %  routine to return Pnnj and Pn,n-1,j when n is odd and n >= 5. Same as Table 4
    %        %  but when n is odd and n >= 5.
    %        %---------------------------------------------------------------------------- */
    %        function dpodd
    %            (
    %            Int32 n,
    %            xpold,
    %            out xp,
    %            out xp1,
    %            Int32[] ipold,
    %            out Int32[] ip,
    %            out Int32[] ip1
    %            )
    % 
    %            xpold = new [360];
    %            xp = new [360];
    %            xp1 = new [360];
    %            ipold = new Int32[360];
    %            ip = new Int32[360];
    %            ip1 = new Int32[360];
    % 
    %            Int32 jx, jxm2, jxm1, n2, j, itemp, jm1, jp1;
    %             gamma, gamma2, xtemp, alpha;
    % 
    %            jx = (n - 1) / 2;
    %            jxm2 = jx - 2;
    %            jxm1 = jx - 1;
    %            n2 = n * 2;
    %            gamma = sqrt(str2num(n2 + 1) * str2num(n2 - 1) / (str2num(n) * str2num(n - 1))) * 0.125;
    %            gamma2 = gamma * 2.0;
    %            xlsum2(gamma * 3.0, xpold(1), -gamma, xpold(2), out xp(1), ipold(1), ipold(2), out ip(1));
    %            j = 1;
    %            while (j <= jxm2)
    % 
    %                jm1 = j - 1;
    %                jp1 = j + 1;
    %                xlsum2(-gamma, xpold[jm1], gamma2, xpold(j), out xtemp, ipold[jm1], ipold(j), out itemp);
    %                xlsum2(1.0, xtemp, -gamma, xpold[jp1], out xp(j), itemp, ipold[jp1], out ip(j));
    %                j = j + 1;
    %            end
    %            xlsum2(-gamma, xpold[jxm2], gamma2, xpold[jxm1], out xp[jxm1], ipold[jxm2], ipold[jxm1], out ip[jxm1]);
    %            xp[jx] = -gamma * xpold[jxm1];
    %            ip[jx] = ipold[jxm1];
    %            xunit(ref xp[jx], ref ip[jx]);
    %            alpha = sqrt(2.0 / str2num(n));
    %            j = 0;
    %            while (j <= jx)
    % 
    %                xp1(j) = str2num(2 * j + 1) * alpha * xp(j);
    %                ip1(j) = ip(j);
    %                xunit(ref xp1(j), ref ip1(j));
    %                j = j + 1;
    %            end
    %        end   % dpodd
    % 
    % 
    %        % ----------------------------------------------------------------------------
    %        %                             gpeven
    %        %
    %        % routine to return Pnmj when n is even. The returned values are (xp0(j), ip0(j)),
    %        % a  X-number vector representing Pnmj. We assume that Pn,m+1,j and Pn,m+2,j are
    %        % externally provided as (xp1(j), ip1(j)) and (xp2(j), ip2(j)), respectively.
    %        % The routine internally calls xunit and xlsum2 provided in Tables 7 and 8 of
    %        % Fukushima (2012a).
    %        % ---------------------------------------------------------------------------- */
    % 
    %        function gpeven
    %            (
    %            Int32 jmax,
    %            Int32 n,
    %            Int32 m,
    %            xp2,
    %            xp1,
    %            out xp0,
    %            Int32[] ip2,
    %            Int32[] ip1,
    %            out Int32[] ip0
    %            )
    % 
    %            xp2 = new [360];
    %            xp1 = new [360];
    %            xp0 = new [360];
    %            ip2 = new Int32[360];
    %            ip1 = new Int32[360];
    %            ip0 = new Int32[360];
    % 
    %            Int32 j, m1, m2, modd;
    %             u, alpha2, beta;
    % 
    %            m1 = m + 1;
    %            m2 = m + 2;
    %            modd = m - str2num(m * 0.5) * 2;
    %            if (m == 0)
    %                u = sqrt(0.5 / (str2num(n) * str2num(n + 1)));
    %            else
    %                u = sqrt(1.0 / (str2num(n - m) * str2num(n + m1)));
    % 
    %                alpha2 = 4.0 * u;
    %                beta = sqrt(str2num(n - m1) * str2num(n + m2)) * u;
    %                xp0(1) = beta * xp2(1);
    %                ip0(1) = ip2(1);
    %                xunit(ref xp0(1), ref ip0(1));
    %                if (modd == 0)
    % 
    %                    j = 1;
    %                    while (j <= jmax)
    % 
    %                        xlsum2(str2num(j) * alpha2, xp1(j), beta, xp2(j), out xp0(j), ip1(j), ip2(j), out ip0(j));
    %                        j = j + 1;
    %                    end
    %                end
    %            else
    % 
    %                j = 1;
    %                while (j <= jmax)
    % 
    %                    xlsum2(-str2num(j) * alpha2, xp1(j), beta, xp2(j), out xp0(j), ip1(j), ip2(j), out ip0(j));
    %                    j = j + 1;
    %                end
    %            end
    %        end  % gpeven
    % 
    % 
    %        % ----------------------------------------------------------------------------
    %        %                                 gpodd
    %        %
    %        % Table 7: Fortran routine to return Pnmj when n is odd. Same as
    %        % Table 6 but when n is odd.
    %        % --------------------------------------------------------------------------- */
    %        function gpodd
    %            (
    %            Int32 jmax,
    %            Int32 n,
    %            Int32 m,
    %            xp2,
    %            xp1,
    %            out xp0,
    %            Int32[] ip2,
    %            Int32[] ip1,
    %            out Int32[] ip0
    %            )
    % 
    %            xp2 = new [360];
    %            xp1 = new [360];
    %            xp0 = new [360];
    %            ip2 = new Int32[360];
    %            ip1 = new Int32[360];
    %            ip0 = new Int32[360];
    % 
    %            Int32 j, m1, m2, modd;
    %             u, alpha, beta;
    % 
    %            m1 = m + 1;
    %            m2 = m + 2;
    %            modd = m - str2num(m * 0.5) * 2;
    %            if (m == 0)
    %                u = sqrt(0.50 / (str2num(n) * str2num(n + 1)));
    %            else
    %                u = sqrt(1.0 / (str2num(n - m) * str2num(n + m1)));
    % 
    %                alpha = 2.0 * u;
    %                beta = sqrt(str2num(n - m1) * str2num(n + m2)) * u;
    %                if (modd == 0)
    % 
    %                    j = 0;
    %                    while (j <= jmax)
    % 
    %                        xlsum2(str2num(2 * j + 1) * alpha, xp1(j), beta, xp2(j), out xp0(j), ip1(j), ip2(j), out ip0(j));
    %                        j = j + 1;
    %                    end
    %                end
    %            else
    % 
    %                j = 0;
    %                while (j <= jmax)
    % 
    %                    xlsum2(-str2num(2 * j + 1) * alpha, xp1(j), beta, xp2(j), out xp0(j), ip1(j), ip2(j), out ip0(j));
    %                    j = j + 1;
    %                end
    %            end
    %        end  % gpodd



% Fukushima combined approach to find matrix of unitalized Legendre polynomials
%
%

% function [ALFArr] = LegPolyFF(recef, latgc, order, unitalized, unitArr)
%     constastro;
% 
%     magr = mag(recef);
% 
%     % initial values
%     ALFArr(1, 1) = 1.0;
%     ALFArr(1, 2) = 0.0;
%     ALFArr(2, 1) = sin(latgc);
%     ALFArr(2, 2) = cos(latgc);
%     m = 2;
%     L = m + 1;
%     x = ALFArr(2, 1);
%     y = ALFArr(2, 2);
% 
%     % find zonal terms
%     for (L = 2; L <= order + 1; L++)
% 
%         % find tesseral and sectoral terms
%         %               for (m = 0; m <= order + 1; m++)
% 
%         f = unitArr[L, m, 0] * sin(latgc);
%         g = -unitArr[L, m, 1];
%         z = f * x + g * y;
%         ALFArr[L, m] = z;
%         y = x;
%         x = z;
%     end
% 
% end  % LegPolyFF



           % /* ----------------------------------------------------------------------------
           % *                                   xfsh2f
           % *
           % *  transrfrom the Cnm, Snm, the 4 fully unitalized spherical harmonic coefficients
           % *  of a given function depend on the spherical surface, to (Akm, Bkm), the
           % *  corresponding Fourier series coefficients of the function. in the program,
           % *  (i) x2f and xunit are the Fortran function/routine to handle X-numbers
           % *      (Fukushima, 2012a, tables 6 and 7), and
           % *  (ii) pinit, dpeven, dpodd, gpeven, and gpodd are the Fortran routines listed
           % *  in Tables 3–7, respectively.
           % *
           % *
           % *
           % *    Fukushima 2018 table xx
           % ---------------------------------------------------------------------------- */
           % 
           % function xfsh2f
           %     (
           %     Int32 nmax,
           %     AstroLib.gravityConst gravData,
           %     out [,] a,
           %     out [,] b
           %     )
           % 
           %     Int32 NX = 100;  % 2200;
           %     Int32 j, m, k, L, jmax, n1;
           %     a = new [NX, NX];
           %     b = new [NX, NX];
           %     %     [,] pja = new [NX, NX];  % test to see the values
           %     jmax = 0;
           % 
           % 
           %     % initialize all to 0
           %     for (m = 0; m <= nmax; m++)
           % 
           %         for (k = 0; k <= nmax; k++)
           % 
           %             a[k, m] = 0.0;
           %             b[k, m] = 0.0;
           %         end
           %     end
           % 
           %     % initialize the first 4x4 values
           %     for (L = 0; L <= 4; L += 2)
           % 
           %         jmax = str2num(L * 0.5);
           % 
           %         for (m = 0; m <= L; m++)
           % 
           %             pinit(L, m, ref xp);
           % 
           %             for (j = 0; j <= jmax; j++)
           % 
           %                 k = 2 * j;
           %                 a[k, m] = a[k, m] + xp(j) * gravData.c[L, m];
           %                 b[k, m] = b[k, m] + xp(j) * gravData.s[L, m];
           %             end
           %         end
           %     end
           % 
           %     % even terms
           %     for (j = 0; j <= jmax; j++)
           % 
           %         ip(j) = 0;
           %     end
           % 
           %     for (L = 6; L <= nmax; L += 2)
           % 
           % 
           %         for (j = 0; j <= jmax; j++)
           % 
           %             xpold(j) = xp(j);
           %             ipold(j) = ip(j);
           %         end
           %         jmax = str2num(L * 0.5);
           %         n1 = L - 1;
           %         dpeven(L, xpold, out xp, out xp1, ipold, out ip, out ip1);
           % 
           %         for (j = 0; j <= jmax; j++)
           % 
           %             k = 2 * j;
           %             pj = x2f(xp(j), ip(j));
           %             %       pja[k, n1] = pj;
           %             a[k, L] = a[k, L] + pj * gravData.c[L, L];
           %             b[k, L] = b[k, L] + pj * gravData.s[L, L];
           %             pj = x2f(xp1(j), ip1(j));
           %             a[k, n1] = a[k, n1] + pj * gravData.c[L, n1];
           %             b[k, n1] = b[k, n1] + pj * gravData.s[L, n1];
           %             xp2(j) = xp(j);
           %             ip2(j) = ip(j);
           %         end
           % 
           %         for (m = L - 2; m >= 0; m -= 1)
           % 
           %             gpeven(jmax, L, m, xp2, xp1, out xp0, ip2, ip1, out ip0);
           % 
           %             for (j = 0; j <= jmax; j++)
           % 
           %                 k = 2 * j;
           %                 pj = x2f(xp0(j), ip0(j));
           %                 %         pja[k, m] = pj;
           %                 a[k, m] = a[k, m] + pj * gravData.c[L, m];
           %                 b[k, m] = b[k, m] + pj * gravData.s[L, m];
           %                 xp2(j) = xp1(j);
           %                 ip2(j) = ip1(j);
           %                 xp1(j) = xp0(j);
           %                 ip1(j) = ip0(j);
           %             end
           %         end
           %     end
           % 
           %     for (L = 1; L <= 3; L += 2)
           % 
           %         jmax = str2num((L - 1) * 0.5);
           % 
           %         for (m = 0; m <= L; m++)
           % 
           %             pinit(L, m, ref xp);
           % 
           %             for (j = 0; j <= jmax; j++)
           % 
           %                 k = 2 * j + 1;
           %                 a[k, m] = a[k, m] + xp(j) * gravData.c[L, m];
           %                 b[k, m] = b[k, m] + xp(j) * gravData.s[L, m];
           %             end
           %         end
           %     end
           % 
           %     % odd terms
           %     for (j = 0; j <= jmax; j++)
           % 
           %         ip(j) = 0;
           %     end
           % 
           %     for (L = 5; L <= nmax; L += 2)
           % 
           % 
           %         for (j = 0; j <= jmax; j++)
           % 
           %             xpold(j) = xp(j);
           %             ipold(j) = ip(j);
           %         end
           %         jmax = str2num((L - 1) * 0.5);
           %         n1 = L - 1;
           %         dpodd(L, xpold, out xp, out xp1, ipold, out ip, out ip1);
           % 
           %         for (j = 0; j <= jmax; j++)
           % 
           %             k = 2 * j + 1;
           %             pj = x2f(xp(j), ip(j));
           %             a[k, L] = a[k, L] + pj * gravData.c[L, L];
           %             b[k, L] = b[k, L] + pj * gravData.s[L, L];
           %             %   pja[k, n] = pj;
           %             pj = x2f(xp1(j), ip1(j));
           %             a[k, n1] = a[k, n1] + pj * gravData.c[L, n1];
           %             b[k, n1] = b[k, n1] + pj * gravData.s[L, n1];
           %             xp2(j) = xp(j);
           %             ip2(j) = ip(j);
           %         end
           % 
           %         for (m = L - 2; m >= 0; m -= 1)
           % 
           %             gpodd(jmax, L, m, xp2, xp1, out xp0, ip2, ip1, out ip0);
           % 
           %             for (j = 0; j <= jmax; j++)
           % 
           %                 k = 2 * j + 1;
           %                 pj = x2f(xp0(j), ip0(j));
           %                 %     pja[k, m] = pj;
           %                 a[k, m] = a[k, m] + pj * gravData.c[L, m];
           %                 b[k, m] = b[k, m] + pj * gravData.s[L, m];
           %                 xp2(j) = xp1(j);
           %                 ip2(j) = ip1(j);
           %                 xp1(j) = xp0(j);
           %                 ip1(j) = ip0(j);
           %             end
           %         end
           %     end
           % end  % xfsh2f
           % 
           % % from 2017 fukishima paper, ALFs
           % %        j n = 2j Pn(1=2) d Pn(1=2)
           % %1 2 􀀀2.7950849718747371205114670859140954E􀀀01 +4.46E􀀀16
           % %2 4 􀀀8.6718750000000000000000000000000019E􀀀01 +1.28E􀀀16
           % %3 8 􀀀3.0362102888840987874508856147660683E􀀀01 +4.89E􀀀16
           % %4 16 􀀀8.6085221363787000197086857609730105E􀀀01 􀀀2.47E􀀀16
           % %5 32 􀀀3.1119962497760147366174972709413071E􀀀01 􀀀1.14E􀀀15
           % %6 64 􀀀8.5832418550243444685033693028444350E􀀀01 􀀀1.97E􀀀16
           % %7 128 􀀀3.1316449107909472026965626279232704E􀀀01 􀀀2.13E􀀀15
           % %8 256 􀀀8.5762286823362136598509949562358587E􀀀01 +4.73E􀀀16
           % %9 512 􀀀3.1365884210353617421696699741828314E􀀀01 +1.94E􀀀15
           % %10 1024 􀀀8.5744308441362078316451963835453876E􀀀01 􀀀1.12E􀀀15
           % %11 2048 􀀀3.1378260213136415436393598399168296E􀀀01 +1.01E􀀀14


           % % GMAT Pines approach
           % %------------------------------------------------------------------------------
           % function FullGeopPines
           %     (
           %      jday,
           %     pos,
           %      latgc,
           %     Int32 nn, Int32 mm,
           %     AstroLib.gravityConst gravData,
           %     out acc
           %     %out [,] gradient
           %     )
           % 
           %     acc = new (4);
           % 
           %     % Int32 XS = fillgradient ? 2 : 1;
           %     % calculate vector components ----------------------------------
           %      magr = sqrt(pos(1) * pos(1) + pos(2) * pos(2) + pos(3) * pos(3));    % Naming scheme from ref (4)
           %      s = pos(1) / magr;
           %      t = pos(2) / magr;
           %      u = pos(3) / magr; % sin(phi), phi = geocentric latitude
           % 
           %     % Calculate values for A -----------------------------------------
           %     ord = 750;
           %     [,] A = new [ord, ord];
           %     [,] N1 = new [ord, ord];
           %     [,] N2 = new [ord, ord];
           %     [,] V = new [ord, ord];
           %     [,] VR01 = new [ord, ord];
           %     [,] VR11 = new [ord, ord];
           %     [,] VR02 = new [ord, ord];
           %     [,] VR12 = new [ord, ord];
           %     [,] VR22 = new [ord, ord];
           %     Re = new [ord];
           %     Im = new [ord];
           %      sqrt2 = sqrt(2.0);
           %     Int32 XS = 2;
           %     u = sin(latgc);
           % 
           %     % get leg poly unitalization numbers (do once)
           %     % all 0
           %     for (Int32 m = 0; m <= mm + 2; ++m)
           % 
           %         for (Int32 L = m + 2; L <= nn + 2; ++L)
           % 
           %             N1[L, m] = sqrt(((2.0 * L + 1) * (2 * L - 1)) / ((L - m) * (L + m)));  %  in denom??
           %             N2[L, m] = sqrt(((2.0 * L + 1) * (L - m - 1) * (L + m - 1)) / ((2.0 * L - 3) * (L + m) * (L - m)));
           %         end
           %     end
           % 
           %     % Nansr
           %     for (Int32 L = 0; L <= nn + 2; ++L)
           % 
           %         V[L, 0] = sqrt((2.0 * (2 * L + 1)));   % Temporary, to make following loop work
           %         for (Int32 m = 1; m <= L + 2 && m <= mm + 2; ++m)
           % 
           %             V[L, m] = V[L, m - 1] / sqrt(((L + m) * (L - m + 1)));  % need real on L-m?
           %         end
           %         V[L, 0] = sqrt((2.0 * L + 1));       % Now set true value
           %     end
           % 
           %     for (Int32 L = 0; L <= nn; ++L)
           %         for (Int32 m = 0; m <= L && m <= mm; ++m)
           % 
           %             % nn = L;
           %             VR01[L, m] = sqrt(((nn - m) * (nn + m + 1)));  % need real on L-m?
           %             VR11[L, m] = sqrt(((2.0 * nn + 1) * (nn + m + 2) * (nn + m + 1)) / ((2.0 * nn + 3)));
           %             VR02[L, m] = sqrt(((nn - m) * (nn - m - 1) * (nn + m + 1) * (nn + m + 2)));  % need real on L-m?
           %             VR12[L, m] = sqrt((2.0 * nn + 1) / (2.0 * nn + 3) * ((nn - m) * (nn + m + 1) * (nn + m + 2) * (nn + m + 3)));  % need real on L-m?
           %             VR22[L, m] = sqrt((2.0 * nn + 1) / (2.0 * nn + 5) * ((nn + m + 1.0) * (nn + m + 2) * (nn + m + 3) * (nn + m + 4)));
           %             if (m == 0)
           % 
           %                 VR01[L, m] /= sqrt2;
           %                 VR11[L, m] /= sqrt2;
           %                 VR02[L, m] /= sqrt2;
           %                 VR12[L, m] /= sqrt2;
           %                 VR22[L, m] /= sqrt2;
           %             end
           %         end
           % 
           %         % generate legendre polynomials - the off-diagonal elements
           %         A(2, 1) = u * sqrt(3.0);
           %         for (Int32 L = 1; L <= nn + XS; ++L)
           %             A[L + 1, L] = u * sqrt(2.0 * L + 3) * A[L, L];
           % 
           %             % apply column-fill recursion formula (Table 2, Row I, Ref.(2))
           %             for (Int32 m = 0; m <= mm + XS; ++m)
           % 
           %                 for (Int32 L = m + 2; L <= nn + XS; ++L)
           %                     A[L, m] = u * N1[L, m] * A[L - 1, m] - N2[L, m] * A[L - 2, m];  % uses anm bnm from fukushima eq 6, 7
           %                     % Ref.(4), Eq.(24)
           %                     if (m == 0)
           %                         Re[m] = 1;
           %                     else
           %                         Re[m] = s * Re[m - 1] - t * Im[m - 1]; % real part of (s + i*t)^m
           %                         if (m == 0)
           %                             Im[m] = 0;
           %                         else
           %                             Im[m] = s * Im[m - 1] + t * Re[m - 1]; % imaginary part of (s + i*t)^m
           %                         end
           % 
           %                         % Now do summation ------------------------------------------------
           %                         % initialize recursion
           %                          FieldRadius = re;
           %                          rho = FieldRadius / magr;
           %                          Factor = mu;
           %                          rho_np1 = -Factor / magr * rho;   % rho(0) ,Ref(4), Eq 26 , factor = mu for gravity
           %                          rho_np2 = rho_np1 * rho;
           %                          a1 = 0;
           %                          a2 = 0;
           %                          a3 = 0;
           %                          a4 = 0;
           %                         for (Int32 L = 1; L <= nn; ++L)
           % 
           %                             rho_np1 *= rho;
           %                             rho_np2 *= rho;
           %                              sum1 = 0;
           %                              sum2 = 0;
           %                              sum3 = 0;
           %                              sum4 = 0;
           % 
           %                             for (Int32 m = 0; m <= L && m <= mm; ++m) % wcs - removed 'm<=L'
           % 
           %                                  Cval = gravData.c[L, m]; % Cnm(jday, L, m);
           %                                  Sval = gravData.s[L, m]; % Snm(jday, L, m);
           %                                 % Pines Equation 27 (Part of)
           %                                  D = (Cval * Re[m] + Sval * Im[m]) * sqrt2;
           %                                  E, F;
           %                                 if (m == 0)
           %                                     E = 0;
           %                                 else
           %                                     E = (Cval * Re[m - 1] + Sval * Im[m - 1]) * sqrt2;
           %                                     if (m == 0)
           %                                         F = 0;
           %                                     else
           %                                         F = (Sval * Re[m - 1] - Cval * Im[m - 1]) * sqrt2;
           % 
           %                                         % Correct for unitalization
           %                                          Avv00 = A[L, m];
           %                                          Avv01 = VR01[L, m] * A[L, m + 1];
           %                                          Avv11 = VR11[L, m] * A[L + 1, m + 1];
           %                                         % Pines Equation 30 and 30b (Part of)
           %                                         sum1 += m * Avv00 * E;
           %                                         sum2 += m * Avv00 * F;
           %                                         sum3 += Avv01 * D;
           %                                         sum4 += Avv11 * D;
           % 
           %                                         % Truncate the gradient at GRADIENT_MAX x GRADIENT_MAX
           %                                         %if (fillgradient)
           %                                         %
           %                                         %    if ((m <= gradientlimit) && (L <= gradientlimit))
           %                                         %
           %                                         %        % Pines Equation 27 (Part of)
           %                                         %        % 2015.09.18 GMT-5295 m<=2  -> m<=1
           %                                         %         G = m <= 1 ? 0 : (Cval * [m - 2] + Sval * Im[m - 2]) * sqrt2;
           %                                         %         H = m <= 1 ? 0 : (Sval * [m - 2] - Cval * Im[m - 2]) * sqrt2;
           %                                         %        % Correct for unitalization
           %                                         %         Avv02 = VR02[L][m] * A[L][m + 2];
           %                                         %         Avv12 = VR12[L][m] * A[L + 1][m + 2];
           %                                         %         Avv22 = VR22[L][m] * A[L + 2][m + 2];
           %                                         %        if (GmatMathUtil::IsNaN(Avv02) || GmatMathUtil::IsInf(Avv02))
           %                                         %            Avv02 = 0.0;  % ************** wcs added ****
           % 
           %                                         %        % Pines Equation 36 (Part of)
           %                                         %        sum11 += m * (m - 1) * Avv00 * G;
           %                                         %        sum12 += m * (m - 1) * Avv00 * H;
           %                                         %        sum13 += m * Avv01 * E;
           %                                         %        sum14 += m * Avv11 * E;
           %                                         %        sum23 += m * Avv01 * F;
           %                                         %        sum24 += m * Avv11 * F;
           %                                         %        sum33 += Avv02 * D;
           %                                         %        sum34 += Avv12 * D;
           %                                         %        sum44 += Avv22 * D;
           %                                         %    end
           %                                         %    else
           %                                         %
           %                                         %        if (matrixTruncationWasPosted == false)
           %                                         %
           %                                         %            MessageInterface::ShowMessage('*** WARNING *** Gradient data '
           % 
           %                                         %                  'for the state transrition matrix and A-matrix '
           % 
           %                                         %                  'computations are truncated at degree and order '
           % 
           %                                         %                  '<= %d.\L', gradientlimit);
           %                                         %            matrixTruncationWasPosted = true;
           %                                         %        end
           %                                         %    end
           %                                         %end
           %                                     end
           %                                     % Pines Equation 30 and 30b (Part of)
           %                                      rr = rho_np1 / FieldRadius;
           %                                     a1 += rr * sum1;
           %                                     a2 += rr * sum2;
           %                                     a3 += rr * sum3;
           %                                     a4 -= rr * sum4;
           %                                     %if (fillgradient)
           %                                     %
           %                                     %    % Pines Equation 36 (Part of)
           %                                     %    a11 += rho_np2 / FieldRadius / FieldRadius * sum11;
           %                                     %    a12 += rho_np2 / FieldRadius / FieldRadius * sum12;
           %                                     %    a13 += rho_np2 / FieldRadius / FieldRadius * sum13;
           %                                     %    a14 -= rho_np2 / FieldRadius / FieldRadius * sum14;
           %                                     %    a23 += rho_np2 / FieldRadius / FieldRadius * sum23;
           %                                     %    a24 -= rho_np2 / FieldRadius / FieldRadius * sum24;
           %                                     %    a33 += rho_np2 / FieldRadius / FieldRadius * sum33;
           %                                     %    a34 -= rho_np2 / FieldRadius / FieldRadius * sum34;
           %                                     %    a44 += rho_np2 / FieldRadius / FieldRadius * sum44;
           %                                     %end
           %                                 end
           % 
           %                                 % Pines Equation 31
           %                                 acc(1) = a1 + a4 * s;
           %                                 acc(2) = a2 + a4 * t;
           %                                 acc(3) = a3 + a4 * u;
           %                                 %if (fillgradient)
           %                                 %
           %                                 %    % Pines Equation 37
           %                                 %    gradient(1, 1) = a11 + s * s * a44 + a4 / magr + 2 * s * a14;
           %                                 %    gradient(2, 2) = -a11 + t * t * a44 + a4 / magr + 2 * t * a24;
           %                                 %    gradient(3, 3) = a33 + u * u * a44 + a4 / magr + 2 * u * a34;
           %                                 %    gradient(1, 2) = gradient(2, 1) = a12 + s * t * a44 + s * a24 + t * a14;
           %                                 %    gradient(1, 3) = gradient(3, 1) = a13 + s * u * a44 + s * a34 + u * a14;
           %                                 %    gradient(2, 3) = gradient(3, 2) = a23 + t * u * a44 + u * a24 + t * a34;
           %                                 %end
           %                             end  % FullGeopPines

% 
% % -----------------------------------------------------------------------------------------------\
% % Gottlieb approach for acceleration
% % gotpot in his nomenclature
% %
% % -----------------------------------------------------------------------------------------------\
% 
% function FullGeopGot(recef, unitArr, order, out legarrGot,out G,out string straccum )
% constastro;
%    %unitArr = new [order + 2, order + 2, 7];
% 
%    magr = mag(recef);
%    Ri = 1.0 / magr;
%    Xovr = recef(1) * Ri;
%    Yovr = recef(2) * Ri;
%    Zovr = recef(3) * Ri;
%    sinlat = Zovr;
%    coslat = cos(asin(sinlat));
%    Reor = re * Ri;
%    Reorn = Reor;
%    muor = mu * Ri;
%    muor2 = muor * Ri;
% 
%    % include two-body or not
%    %if (Want_Central_force == true)
%    %    Sum_Init = 1;
%    %else
%    % note, 1 makes the two body pretty close, except for the 1st component
%    Sum_Init = 0;
% 
%    % initial values
%    % ctrigArr(1) = 1.0;
%    ctrigArr(2) = Xovr;
%    %  strigArr(1) = 0.0;
%    strigArr(2) = Yovr;
%    Sumh = 0.0;
%    Sumj = 0.0;
%    Sumk = 0.0;
%    Sumgam = Sum_Init;
% 
%    % unitArr(L, m, 0) xi Gottlieb eta
%    % unitArr(L, m, 1) eta Gottlieb zeta
%    % unitArr(L, m, 2) alpha Gottlieb alpha
%    % unitArr(L, m, 3) beta Gottlieb beta
%    % unitArr(L, m, 5) delta Gottlieb zn
%    legarrGot(1, 1) = 1.0;
%    legarrGot(1, 2) = 0.0;
%    legarrGot(2, 1) = sqrt(3) * sinlat;
%    legarrGot(2, 2) = sqrt(3); % * coslat;
% 
%    for n = 2:order
%         % get the power for each n
%        Reorn = Reorn * Reor;
%        %pn = legPoly(n, 1);
%        cn = gravData.cNor(n, 1);
%        %sn = gravData.sNor(n, 1);
% 
%        nm1 = n - 1;
%        nm2 = n - 2;
% 
%        % eq 3-17, eq 7-14  alpha(n) beta(n)
%        legarrGot(n,1) = sinlat * unitArr[n, 0, 2] * legarrGot[nm1, 0] - unitArr[n, 0, 3] * legarrGot[nm2, 0];
%        % inner diagonal eq 7-16
%        % n-1,n-2, 6, not 5, no nm1
%        legarrGot[n, nm1] = unitArr[n - 1, nm2, 6] * sinlat * legarrGot[n, n];
%        %      legPoly[n, nm1] = unitArr[n, nm1, 7] * sinlat * legPoly[n, n];
%        % diagonal eq 7-8
%        legarrGot[n, n] = unitArr[n, n, 4] * coslat * legarrGot[nm1, nm1];
% 
%        Sumh_N = unitArr[1, 0, 6] * legarrGot(n, 1) * cn;  % 0 by 2016 paper
%        Sumgam_N = legarrGot(n, 1) * cn * (n + 1);  %
% 
%        if (order > 0)
% 
%            for m = 1:= nm2
% 
%                % eq 3-18, eq 7-12   xin(m) eta(m)
%                legarrGot[n, m] = unitArr[n, m, 0] * sinlat * legarrGot[nm1, m] - unitArr[n, m, 1] * legarrGot[nm2, m];
%            end
%            % got all the Legendre functions now
% 
%            Sumj_N = 0.0;
%            Sumk_N = 0.0;
%            ctrigArr(n) = ctrigArr(2) * ctrigArr[nm1] - strigArr(2) * strigArr[nm1]; % mm1????
%            strigArr(n) = strigArr(2) * ctrigArr[nm1] + ctrigArr(2) * strigArr[nm1];
% 
%            if (n < order)
%                Lim = n;
%            else
%                Lim = order;
% 
%                for (m = 1; m <= Lim; m++)
% 
%                    mm1 = m - 1;
%                    mp1 = m + 1;
%                    npmp1 = (n + mp1);  %
%                    pnm = legarrGot[n, m];
%                    cnm = gravData.cNor[n, m];
%                    snm = gravData.sNor[n, m];
%                    %ctmm1 = gravData.cNor[n, mm1];
%                    %stmm1 = gravData.sNor[n, mm1];
% 
%                    Mxpnm = m * pnm;  %
%                    BnmVal = cnm * ctrigArr[m] + snm * strigArr[m];
%                    Sumh_N = Sumh_N + legarrGot[n, mp1] * BnmVal * unitArr[n, m, 6];  % zn(m)
%                    Sumgam_N = Sumgam_N + npmp1 * pnm * BnmVal;
%                    Sumj_N = Sumj_N + Mxpnm * (cnm * ctrigArr[m] + snm * strigArr[m]);
%                    Sumk_N = Sumk_N - Mxpnm * (cnm * strigArr[m] - snm * ctrigArr[m]);
%                end   % for through m
% 
%                Sumj = Sumj + Reorn * Sumj_N;
%                Sumk = Sumk + Reorn * Sumk_N;
%            end  % if order > 0
% 
%            % ---- Sums bleow here have values when N m = 0
%            Sumh = Sumh + Reorn * Sumh_N;
%            Sumgam = Sumgam + Reorn * Sumgam_N;
%        end  % loop
% 
%        Lambda = Sumgam + sinlat * Sumh;
%        G(1) = -muor2 * (Lambda * Xovr - Sumj);
%        G(2) = -muor2 * (Lambda * Yovr - Sumk);
%        G(3) = -muor2 * (Lambda * Zovr - Sumh);
% 
%        % if (show == 'y')
% 
%        straccum = straccum + 'Gottlieb case nonspherical, no two-body ---------- ' + '\n';
%        straccum = straccum + 'legarrGot 4 0   ' + legarrGot(5, 1), '  4 1   '
%        + legarrGot(5, 2), '  4 4   ' + legarrGot(5, 5), '\n';
%        straccum = straccum + 'legarrGot 5 0   ' + legarrGot(6, 1), '  5 1   '
%        + ctrigArr(3), '  Tan   ' + strigArr(3), '\n';
%        straccum = straccum + 'legarrGot' + order + ' 0   ' + legarrGot[order, 0], '  ' + order + ' 1   '
%        + legarrGot[order, 1], ' + order + ' + legarrGot[order, order], '\n';
%        %straccum = straccum + 'trigarr ' + order + ' Sin  ' + trigArr[order, 0], '  Cos   '
%        %    + trigArr[order, 1], '  Tan   ' + trigArr[order, 3], '\n';
%        straccum = straccum + 'apertGot ecef ' + order, order, G(1), '     '
%        + G(2), '     ' + G(3), '\n';
%    end
% 
% end  % FullGeopGot;
% 
% 
% 
% function testproporbit(fid)
%     % 152 is arbitrary
%     constastro;
%     rad = 180.0 / pi;              % deg to rad
%     conv = pi / (180.0 * 3600.0);  % ' to rad
%     fileLoc = 'D:\Codes\LIBRARY\Matlab\';
%     fid = fopen(append(fileLoc, 'testallmat.out'), 'wt');
%     fid1 = fopen(append(fileLoc, 'testallmatsum.out'), 'wt');
% 
%     % ------------------ BOOK EXAMPLE ----------------------------------------------
%     % ------------------------------- initial state -------------------------------
%     reci = [ -605.79079600; -5870.23042200; 3493.05191600 ];
%     veci = [ -1.568251000; -3.702348000; -6.479485000 ];
%     % prout initial conditions
%     fprintf(fida,'reci %15.11f %15.11f %15.11f %15.11f %15.11f %15.11f\n', reci(1), reci(2), reci(3), veci(1), veci(2), veci(3));
%     cd = 2.2;
%     cr = 1.2;
%     area = 40.0;     % m^2
%     mass = 1000.0;   % kg
% 
%     % ------------------------------- establish time parameters -------------------------------
%     fileLoc = 'D:\Codes\LIBRARY\DataLib\';
%     [iau80arr] = iau80in(fileLoc);
% 
%     year = 2020;
%     mon = 2;
%     day = 18;
%     hr = 15;
%     minute = 8;
%     second = 47.23847;
%     [jdutc, jdFutc] = jday(year, mon, day, hr, minute, second);  % utc
% 
%     % these vaues are not consistent wih the EOP files - are they interpolated already????
%     % no. change to use actual at 0000 values, rerun probelm.
%     %2020 02 18 58897  0.030655  0.336009 -0.1990725  0.0000639 -0.108041 -0.007459  0.000315 -0.000055  37
%     %2020 02 19 58898  0.030313  0.337617 - 0.1991016  0.0000235 - 0.107939 - 0.007476  0.000324 - 0.000036  37
% 
%        xp = 0.030655 * conv;
%        yp = 0.336009 * conv;
%        lod = 0.0000639;
%        ddpsi = -0.108041 * conv;  % ' to rad
%        ddeps = -0.007459 * conv;
%        ddx = 0.000315 * conv;     % ' to rad
%        ddy = 0.000055 * conv;
%        dut1 = -0.1990725;   % second
%        dat = 37;            % second
%        opt = '80';
%        eqeterms = 2;
%        jdtt = jdutc;
%        jdftt = jdFutc + (dat + 32.184)/86400.0;
% 
%        % method to do calculations in
%        unitialized = 'y';
%        fprintf(fida,'unitialized = %c\n', unitialized);
% 
%      fprintf(fida,'%4i  %3i  %3i  %2i:%2i:%6.4f\n', year, mon, day, hr, minute, second);
%        fprintf(fida,'dat  %3i lod %8.6f\n', dat, lod);
%     fprintf(fida,'jdutc  %16.12f\n',jdutc + jdFutc );
%        fprintf(fida,'xp %15.11f  yp %15.11f arcsec\n', (xp / conv), (yp / conv));
%        fprintf(fida,'dpsi %15.11f  deps %15.11f arcsec\n', (ddpsi / conv), (ddeps / conv));
%        fprintf(fida,'ddx %15.11f  ddy %15.11f arcsec\n', ddx, ddy);
%        jdut1 = jdutc + jdFutc + dut1 / 86400.0;
%        fprintf(fida,'jdut1 %16.12f\n', jdut1);
% 
%        % watch if getting tdb that j2000 is also tdb
%        ttt = (jdutc + jdFutc + (dat + 32.184) / 86400.0 - 2451545.0) / 36525.0;
%        fprintf(fida,'jdttt ' + (jdutc + jdFutc + (dat + 32.184) / 86400.0), ' ttt ' + ttt, '\n');
% 
%        fundarg(ttt, opt, out fArgs);
% 
%        tmptdb = (dat + 32.184 + 0.001657 * sin(628.3076 * ttt + 6.2401)
%        + 0.000022 * sin(575.3385 * ttt + 4.2970)
%        + 0.000014 * sin(1256.6152 * ttt + 6.1969)
%        + 0.000005 * sin(606.9777 * ttt + 4.0212)
%        + 0.000005 * sin(52.9691 * ttt + 0.4444)
%        + 0.000002 * sin(21.3299 * ttt + 5.5431)
%        + 0.000010 * ttt * sin(628.3076 * ttt + 4.2490)) / 86400.0;  % USNO circ(14)
%        jday(year, mon, day, hr, minute, second + tmptdb, out jdtdb, out jdFtdb);
%        ttdb = (jdtdb + jdFtdb - 2451545.0) / 36525.0;
%        fprintf(fida,'jdttb ' + (jdtdb + jdFtdb), ' ttdb ' + ttdb);
% 
%        % get reduction matrices
%        deltapsi = 0.0;
%        meaneps = 0.0;
% 
%        omegaearth(1) = 0.0;
%        omegaearth(2) = 0.0;
%        omegaearth(3) = earthrot * (1.0 - lod / 86400.0);
% 
%        prec = precess(ttt, opt, out psia, out wa, out epsa, out chia);
%        nut = nutation(ttt, ddpsi, ddeps, iau80arr, opt, fArgs, out deltapsi, out deltaeps, out trueeps, out meaneps);
%        st = sidereal(jdut1, deltapsi, meaneps, fArgs, lod, eqeterms, opt);
%        pm = polarm(xp, yp, ttt, opt);
% 
%        %% ---- perform transrformations eci to ecef
%        pmp = mattransr(pm, 3);
%        stp = mattransr(st, 3);
%        nutp = mattransr(nut, 3);
%        precp = mattransr(prec, 3);
%        temp = matmult(pmp, stp, 3, 3, 3);
%        temp1 = matmult(temp, nutp, 3, 3, 3);
%        transreci2ecef = matmult(temp1, precp, 3, 3, 3);
%        recef = matvecmult(transreci2ecef, reci, 3);
%        fprintf(fida,'recef  ' + recef(1), recef(2), recef(3),
%        'v  ' + vecef(1), vecef(2), vecef(3));
% 
%        %----perform transrformations ecef to eci
%        % note the rotations occur only for velocity so the full transrformation is fine here
%        transrecef2eci = mattransr(transreci2ecef, 3);
%        reci = matvecmult(transrecef2eci, recef, 3);
%        fprintf(fida,'reci  ' + reci(1), reci(2), reci(3),
%        'v  ' + veci(1), veci(2), veci(3));
% 
%        eci_ecef(ref reci, ref veci, MathTimeLib.Edirection.eto, ref recef, ref vecef,
%        AstroLib.EOpt.e80, iau80arr, iau06arr,
%        jdtt, jdftt, jdut1, jdxysstart, lod, xp, yp, ddpsi, ddeps, ddx, ddy);
% 
%        % prout initial conditions
%        fprintf(fida,'reci  ' + reci(1), reci(2), reci(3),
%        'v  ' + veci(1), veci(2), veci(3));
%        fprintf(fida,'recef  ' + recef(1), recef(2), recef(3),
%        'v  ' + vecef(1), vecef(2), vecef(3));
% 
%        [latgc, latgd, lon, hellp] = ecef2ll(recef);
%        % or
%        latgc = 52.0 / rad;  % 52 and 34
%        lon = 5.0 / rad;
%        alt = 6880.0;
%        site(latgd, lon, alt, out recef, out vecef);
%        fprintf(fida,'new site loc');
%        fprintf(fida,'recef  ' + recef(1), recef(2), recef(3),
%        'v  ' + vecef(1), vecef(2), vecef(3));
% 
%        % ---------------------------------------------------------------------------------------------
%        % ------------------------------------ GRAVITY FIELD --------------------------------------
%        % get past text in each file
%        %if (fname.Contains('GEM'))    % GEM10bununit36.grv, GEMT3unit50.grv
%        %    startKtr = 17;
%        %if (fname.Contains('EGM-96')) % EGM-96unit70.grv
%        %    startKtr = 73;
%        %if (fname.Contains('EGM-08')) % EGM-08unit100.grv
%        %startKtr = 83;  % or 21 for the larger file... which has gfc in the first col too
%        % fully unitalized, 4415, .1363, order 100
%        %string fname = 'D:/Dataorig/Gravity/EGM-08unit100.grv';  % 83
%        % fully unitalized, 4415, .1363, order 2190
%        %startKtr = 21;  % or 21 for the larger file... which has gfc in the first col too
%        fname = 'D:/Dataorig/Gravity/EGM2008_to2190_TideFree.txt';
%        % fully unitalized, 4415, .1363, order 360
%        %string fname = 'D:/Dataorig/Gravity/GGM03C-Data.txt';
% 
%        % EGM2008_to2190_TideFree.txt
%        % L    m               c                      s                 csigma              ssigma
%        % 2    0   -0.484165143790815D-03    0.000000000000000D+00    0.7481239490D-11    0.0000000000D+00
%        % 2    1   -0.206615509074176D-09    0.138441389137979D-08    0.7063781502D-11    0.7348347201D-11
%        % 2    2    0.243938357328313D-05   -0.140027370385934D-05    0.7230231722D-11    0.7425816951D-11
%        [gravarr] = readgravityfield(fname, 'y');
%        fprintf(fida,'\nread in gravity field %s ---------------\n', fname);
%        fprintf(fida,'\ncoefficents --------------- ');
%        fprintf(fida,'c  2  0  ' + gravData.c(3, 1), ' s ' + gravData.s(3, 1));
%        fprintf(fida,'c  4  0  ' + gravData.c(5, 1), ' s ' + gravData.s(5, 1));
%        fprintf(fida,'c  4  4  ' + gravData.c(5, 5), ' s ' + gravData.s(5, 5));
%        fprintf(fida,'c 21  1 ' + gravData.c[21, 1], ' s ' + gravData.s[21, 1]);
%        fprintf(fida,'\nunitalized coefficents --------------- ');
%        fprintf(fida,'c  2  0  ' + gravData.cNor(3, 1), ' s ' + gravData.sNor(3, 1));
%        fprintf(fida,'c  4  0  ' + gravData.cNor(5, 1), ' s ' + gravData.sNor(5, 1));
%        fprintf(fida,'c  4  4  ' + gravData.cNor(5, 5), ' s ' + gravData.sNor(5, 5));
%        fprintf(fida,'c 21  1 ' + gravData.cNor[21, 1], ' s ' + gravData.sNor[21, 1]);
%        fprintf(fida,'c 500  1 ' + gravData.cNor[500, 1], ' s ' + gravData.sNor[500, 1]);
% 
%        % --------------------------------------------------------------------------------------------------
%        % calculate legendre polynomials
% 
%        fprintf(fida,'\nCalculate Legendre polynomial recursions Unn and Nor  --------------- ');
%        degree = 500;
%        order = 500;
%        % GTDS version
%        % does with  ununitalized elements, then unitalized from there. But ununitalized only go to about 170
%        LegPolyG(latgc, order, unitalized, convArr, unitArr, out LegArrGU, out LegArrGN);
% 
%        % Gottlieb version
%        LegPolyGot(latgc, order, unitalized, convArr, unitArr, out LegArrGotU, out LegArrGotN);
% 
%        % Montenbruck version
%        LegPolyM(latgc, order, unitalized, convArr, out LegArrMU, out LegArrMN);
% 
%        % Geodyn version
%        geodynlegp(latgc, degree, order, out LegArrOU, out LegArrON);
% 
%        % Exact values
%        LegPolyEx(latgc, order, out LegArrEx);
% 
%        % Fukushima approach do as 1-d arrays for now
%        LegPolyF(latgc, order, 'y', unitArr, out LegArrF);
%        %pmm = new (9);
%        %psm = new (9);
%        %Int32[] ipsm = new Int32(9);
%        % get the values in X-numbers
%        %alfsx(cos(latgc), 6, unitArr, out psm, out ipsm);
%        %alfmx(sin(latgc), 3, 6, unitArr, psm(4), ipsm(4), out pmm);
% 
% 
%        sumdr1 = 0.0;
%        sumdr2 = 0.0;
%        sumdr3 = 0.0;
%        fprintf(fida,'\nwrite out unitalized Legendre polynomials --------------- ');
% 
%        % order xxxxxxxxxxxxxxxxxx
%        for (L = 0; L <= 130; L++)
% 
%            string tempstr1 = 'MN  ';  % montenbruck
%            string tempstr2 = 'GN  ';  % gtds
%            string tempstr3 = 'MU  ';
%            string tempstr3a = 'LU  ';  % exact
%            string tempstr4 = 'OU  ';  % geodyn\
%            string tempstr5 = 'GtN ';  % gottlieb\
%            string tempstr6 = 'GtU ';  % gottlieb
%            string tempstr7 = 'FN  ';  % Fukushima, test ones
%            stopL = L;
%            for (m = 0; m <= stopL; m++)
% 
%                tempstr1 = tempstr1, L, '  ' + m, '   ' + LegArrMN[L, m];
%                tempstr2 = tempstr2, L, '  ' + m, '   ' + LegArrGN[L, m];
%                tempstr5 = tempstr5, L, '  ' + m, '   ' + LegArrGotN[L, m];
%                tempstr7 = tempstr7, L, '  ' + m, '   ' + LegArrF[L, m];
%                tempstr3 = tempstr3, L, '  ' + m, '   ' + LegArrMU[L, m];
%                tempstr3a = tempstr3a, L, '  ' + m, '   ' + LegArrEx[L, m];
%                tempstr6 = tempstr6, L, '  ' + m, '   ' + LegArrGotU[L, m];
%                tempstr4 = tempstr4, L, '  ' + m, '   ' + LegArrOU[L + 1, m + 1];
%                % check error values
%                dr1 = 100.0 * (LegArrF[L, m] - LegArrGotN[L, m]) / LegArrF[L, m];
%                dr2 = 100.0 * (LegArrF[L, m] - LegArrGN[L, m]) / LegArrF[L, m];
%                dr3 = 100.0 * (LegArrF[L, m] - LegArrMN[L, m]) / LegArrF[L, m];
%                sumdr1 = sumdr1 + dr1;
%                sumdr2 = sumdr2 + dr2;
%                sumdr3 = sumdr3 + dr3;
%                errstr = errstr + '\n' + L, '  ' + m, '   ' + dr1
%                , dr2, dr3;
%            end
%            % unitalized ones
%            fprintf(fida,tempstr2);
%            fprintf(fida,tempstr1);
%            fprintf(fida,tempstr5);
%            fprintf(fida,tempstr7 + '\n');
%            % ununitalized ones
%            fprintf(fid,tempstr3);
%            fprintf(fid,tempstr3a);
%            fprintf(fid,tempstr6);
%            fprintf(fid,tempstr4 + '\n');
%        end
%        strbuildplot.AppendLine(errstr);
% 
%        % -------------------- now accelerations -----------------------------------------------------
%        fprintf(fida,'\naccelerations --------------- ');
%        string straccum = '';
% 
%        %order = 4;
%        order = 120; % 10;
%        % GTDS acceleration for non-spherical portion
%        FullGeopG(recef, order, unitalized, convArr, unitArr, gravData, out aPertG, 'y', out straccum);
%        fprintf(fida,straccum);
%        aeci = matvecmult(transrecef2eci, aPertG, 3);
%        straccum = straccum + 'apertG eci  ' + order, order, aeci(1), '     '
%        + aeci(2), '     ' + aeci(3), '\n';
%        fprintf(fida,straccum);
% 
%        FullGeopG(recef, order, unitalized, convArr, unitArr, gravData, out aPertG, 'n', out straccum);
%        fprintf(fida,straccum);
% 
%        % Montenbruck acceleration
%        FullGeopM(recef, order, unitalized, convArr, gravData, out aPertM, 'y', out straccum);
%        fprintf(fida,straccum);
%        aeci = matvecmult(transrecef2eci, aPertM, 3);
%        straccum = straccum + 'apertM eci  ' + order, order, aeci(1), '     '
%        + aeci(2), '     ' + aeci(3), '\n';
%        fprintf(fida,straccum);
% 
%        FullGeopM(recef, order, unitalized, convArr, gravData, out aPertM, 'n', out straccum);
%        fprintf(fida,straccum);
% 
%        % Montenbruck code acceleration
%        FullGeopMC(recef, order, unitalized, convArr, gravData, out aPertM1, 'y', out straccum);
%        fprintf(fida,straccum);
%        aeci = matvecmult(transrecef2eci, aPertM1, 3);
%        straccum = straccum + 'apertM1 eci ' + order, order, aeci(1), '     '
%        + aeci(2), '     ' + aeci(3), '\n';
%        fprintf(fida,straccum);
% 
%        FullGeopMC(recef, order, unitalized, convArr, gravData, out aPertM1, 'n', out straccum);
%        fprintf(fida,straccum);
% 
%        % Gottlieb acceleration
%        fprintf(fida,'Gottlieb acceleration ');
%        G = new (4);
%        aPertGt = new (4);
%        FullGeopGot(gravData, recef, unitArr, order, out LegArrGott, out G, out straccum);
%        fprintf(fida,straccum);
% 
%        % Fukushima acceleration
%        fprintf(fida,'Fukushima acceleration ');
%        %LegPolyFF(recef, latgc, order, 'y', unitArr, gravData, out LegArrF);
%        [,] a = new [360, 360];
%        [,] b = new [360, 360];
%        xfsh2f(80, gravData, out a, out b);
%        fprintf(fida,a(3, 1));
%        fprintf(fida,'a  2  0  ' + a(3, 1), ' b ' + b(3, 1));
%        fprintf(fida,'a  2  1  ' + a(3, 2), ' b ' + b(3, 2));
%        fprintf(fida,'a  4  0  ' + a(5, 1), ' b ' + b(5, 1));
%        fprintf(fida,'a  4  1  ' + a(5, 1), ' b ' + b(5, 2));
%        fprintf(fida,'a  4  4  ' + a(5, 5), ' b ' + b(5, 5));
%        fprintf(fida,'a 10 10  ' + a[10, 0], ' b ' + b[10, 0]);
% fprintf(fida,'a 21  1 ' + a[21, 1], ' b ' + b[21, 1]);
% 
% % Pines approach
% fprintf(fida,'Pines acceleration ');
% FullGeopPines(jdutc, recef, latgc, order, order, gravData, out aeci);
% fprintf(fida,'apertP    4 4   ' + aeci(1), '     ' + aeci(2), '     ' + aeci(3));
% 
% fprintf(fida,straccum);
% fprintf(fida,'\ngravity field ' + fname, order, ' --------------- ');
% fprintf(fida,' summary accelerations ----------------------------------------------- ');
% fprintf(fida,'apertG bf  ' + order, order, aPertG(1), '     ' + aPertG(2), '     ' + aPertG(3));
% fprintf(fida,'apertM bf  ' + order, order, aPertM(1), '     ' + aPertM(2), '     ' + aPertM(3));
% fprintf(fida,'apertMC bf ' + order, order, aPertM1(1), '     ' + aPertM1(2), '     ' + aPertM1(3));
% fprintf(fida,'apertGt bf ' + order, order, G(1), '     ' + G(2), '     ' + G(3));
% 
% aPertG = matvecmult(transrecef2eci, aPertG, 3);
% aPertM = matvecmult(transrecef2eci, aPertM, 3);
% aPertM1 = matvecmult(transrecef2eci, aPertM1, 3);
% aPertGt = matvecmult(transrecef2eci, G, 3);
% fprintf(fida,'apertG  eci ' + order, order, aPertG(1), '     ' + aPertG(2), '     ' + aPertG(3));
% fprintf(fida,'apertM  eci ' + order, order, aPertM(1), '     ' + aPertM(2), '     ' + aPertM(3));
% fprintf(fida,'apertMC eci ' + order, order, aPertM1(1), '     ' + aPertM1(2), '     ' + aPertM1(3));
% fprintf(fida,'apertGt eci ' + order, order, aPertGt(1), '     ' + aPertGt(2), '     ' + aPertGt(3));
% 
% fprintf(fida,'STK ansr 4x4         -0.0000003723020	-0.0000031362090   	-0.0000102647170\n');  % no 2-body
% 
% 
% % -------------------------- add in two body term since full geop is only disturbing part
% eci_ecef(ref reci, ref veci, MathTimeLib.Edirection.efrom, ref recef, ref vecef,
% AstroLib.EOpt.e80, iau80arr, iau06arr,
% jdtt, jdftt, jdut1, jdxysstart, lod, xp, yp, ddpsi, ddeps, ddx, ddy);
% aeci2(1) = -mu * reci(1) / (Math.Pow(mag(reci), 3));
% aeci2(2) = -mu * reci(2) / (Math.Pow(mag(reci), 3));
% aeci2(3) = -mu * reci(3) / (Math.Pow(mag(reci), 3));
% fprintf(fida,'a2body      ' + aeci2(1), '     ' + aeci2(2), '     ' + aeci2(3));
% 
% aPertG(1) = aPertG(1) + aeci2(1);
% aPertG(2) = aPertG(2) + aeci2(2);
% aPertG(3) = aPertG(3) + aeci2(3);
% 
% temm = new (4);
% temm(1) = aPertG(1);
% temm(2) = aPertG(2);
% temm(3) = aPertG(3);
% 
% aPertM(1) = aPertM(1) + aeci2(1);
% aPertM(2) = aPertM(2) + aeci2(2);
% aPertM(3) = aPertM(3) + aeci2(3);
% 
% aPertM1(1) = aPertM1(1) + aeci2(1);
% aPertM1(2) = aPertM1(2) + aeci2(2);
% aPertM1(3) = aPertM1(3) + aeci2(3);
% 
% fprintf(fida,' now with two body included');
% fprintf(fida,'apertG ' + order, order, aPertG(1), '     ' + aPertG(2), '     ' + aPertG(3));
% fprintf(fida,'apertM ' + order, order, aPertM(1), '     ' + aPertM(2), '     ' + aPertM(3));
% fprintf(fida,'apertMC ' + order, order, aPertM1(1), '     ' + aPertM1(2), '     ' + aPertM1(3));
% fprintf(fida,'STK ansr 4x4 w2   0.0007483593980          0.0072522125910         -0.0043275195170\n');  % no 2-body
% %                    4x4 j2000   0.00074835849281         0.00725221243453        -0.00432751993509
% %                    4x4 icrf    0.00074835939828         0.00725221259059        -0.00432751951698
% %                    all j2000   0.00074845403274         0.00725223127396        -0.00432750265312
% %                    all icrf    0.00074845493821         0.00725223143002        -0.00432750223499
% 
% 
% fprintf(fida,'------------------ find drag acceleration');
%  density = 1.5e-12;  % kg / m3
%  magv = mag(vecef);
% vrel(1) = vecef(1); % vecef unital is veci to tod, then - wxr
% vrel(2) = vecef(2);
% vrel(3) = vecef(3);
% fprintf(fida,' vrel ' + vrel(1), vrel(2), vrel(3));
% %                 kg / m3        m2  /  kg     km / s  km / s
% adrag(1) = -0.5 * density * cd * area / mass * magv * vrel(1) * 1000.0;  % simplify vel, get units to km/s2
% adrag(2) = -0.5 * density * cd * area / mass * magv * vrel(2) * 1000.0;  % simplify vel, get units to km/s2
% adrag(3) = -0.5 * density * cd * area / mass * magv * vrel(3) * 1000.0;  % simplify vel, get units to km/s2
% 
% fprintf(fida,' adrag ecef' + adrag(1), adrag(2), adrag(3));
% 
% fprintf(fida,' agrav + drag ecef' +(temm(1)+ adrag(1)),
% (temm(2) + adrag(2)), (temm(3) + adrag(3)));
% 
% 
% transrecef2eci = matmult(temp1, pm, 3, 3, 3);
% aeci = matvecmult(transrecef2eci, adrag, 3);
% fprintf(fida,' adrag eci ' + aeci(1), aeci(2), aeci(3));
% fprintf(fida,'ansr drag JR spline      0.0000000001040	0.0000000002090	0.0000000003550\n');
% fprintf(fida,'ansr drag JR daily       0.0000000000840	0.0000000001720	0.0000000002900\n');
% fprintf(fida,'ansr drag MSIS daily     0.0000000000730	0.0000000001510	0.0000000002530\n');
% 
% temmm = new (4);
% temmm(1) = temm(1) + adrag(1);
% temmm(2) = temm(2) + adrag(2);
% temmm(3) = temm(3) + adrag(3);
% aeci = matvecmult(transrecef2eci, temmm, 3);
% fprintf(fida,' agrav+drag eci ' + aeci(1), aeci(2), aeci(3));
% 
% fprintf(fida,' ------------------ find third body acceleration');
% AstroLib.jpldedataClass[] jpldearr = jpldearr;
%  musun, mumoon, rsmag, rmmag;
% musun = 1.32712428e11;    % km3 / s2
% mumoon = 4902.799;        % km3 / s2
% infilename = append('D:\Codes\LIBRARY\DataLib\', 'sunmooneph_430t.txt');
% [jpldearr, jdjpldestart, jdjpldestartFrac] = initjplde(infilename);
% 
% % sun
% findjpldeparam(jdtdb, jdFtdb, 's', jpldearr, jdtdbjplstart, out rsun, out rsmag, out rmoon, out rmmag);
% % stk value (chk that tdb is argument)
% rsuns = [ 126916355.384390; -69567131.339884; -30163629.424510 ];
% % JPL ansr  2020  2 18  M          0.6306
% rmoonj = [ 14462.2967; -357096.9762; -151599.3021 ];
% %JPL ansr  2020  2 18 15:08:47.23847 S       0.6306
% rsunj = [ 126921698.4134; -69564121.8695; -30156263.9220 ];
% 
% addvec(1.0, rsuns, -1.0, rsun, out tempvec1);
% fprintf(fida,' diff rsun stk-mine ' + tempvec1(1), tempvec1(2),
% tempvec1(3), mag(tempvec1));
% addvec(1.0, rsunj, -1.0, rsun, out tempvec1);
% fprintf(fida,' diff rsun jpl-mine ' + tempvec1(1), tempvec1(2),
% tempvec1(3), mag(tempvec1));
% addvec(1.0, rsuns, -1.0, rsunj, out tempvec1);
% fprintf(fida,' diff rsun stk-jpl  ' + tempvec1(1), tempvec1(2),
% tempvec1(3), mag(tempvec1));
% addvec(1.0, rmoonj, -1.0, rmoon, out tempvec1);
% fprintf(fida,' diff rmoon jpl-mine ' + tempvec1(1), tempvec1(2),
% tempvec1(3), mag(tempvec1));
% fprintf(fida,' rsun  ' + rsun(1), rsun(2), rsun(3));
% fprintf(fida,' rmoon ' + rmoon(1), rmoon(2), rmoon(3));
% 
%  mu3 = musun;
% rsat3(1) = rsun(1) - reci(1);
% rsat3(2) = rsun(2) - reci(2);
% rsat3(3) = rsun(3) - reci(3);
%  magrsat3 = mag(rsat3);
% rearth3(1) = rsun(1);
% rearth3(2) = rsun(2);
% rearth3(3) = rsun(3);
%  magrearth3 = mag(rearth3);
% athirdbody(1) = mu3 * (rsat3(1) / Math.Pow(magrsat3, 3) - rearth3(1) / Math.Pow(magrearth3, 3));
% athirdbody(2) = mu3 * (rsat3(2) / Math.Pow(magrsat3, 3) - rearth3(2) / Math.Pow(magrearth3, 3));
% athirdbody(3) = mu3 * (rsat3(3) / Math.Pow(magrsat3, 3) - rearth3(3) / Math.Pow(magrearth3, 3));
% fprintf(fida,' a3bodyS  eci ' + athirdbody(1), athirdbody(2), athirdbody(3));
% athirdbody2(1) = -mu3 / Math.Pow(magrearth3, 3) * (rearth3(1) - 3.0 * rearth3(1) * (dot(reci, rearth3) / Math.Pow(magrearth3, 2))
% - 7.5 * rearth3(1) * Math.Pow(((dot(reci, rearth3) / Math.Pow(magrearth3, 2))), 2));
% athirdbody2(2) = -mu3 / Math.Pow(magrearth3, 3) * (rearth3(2) - 3.0 * rearth3(2) * (dot(reci, rearth3) / Math.Pow(magrearth3, 2))
% - 7.5 * rearth3(2) * Math.Pow(((dot(reci, rearth3) / Math.Pow(magrearth3, 2))), 2));
% athirdbody2(3) = -mu3 / Math.Pow(magrearth3, 3) * (rearth3(3) - 3.0 * rearth3(3) * (dot(reci, rearth3) / Math.Pow(magrearth3, 2))
% - 7.5 * rearth3(3) * Math.Pow(((dot(reci, rearth3) / Math.Pow(magrearth3, 2))), 2));
% fprintf(fida,' a3bodyS2 eci' + athirdbody2(1), athirdbody2(2), athirdbody2(3));
% q = (Math.Pow(mag(reci), 2) + 2.0 * dot(reci, rsat3)) *
% (Math.Pow(magrearth3, 2) + magrearth3 * magrsat3 + Math.Pow(magrsat3, 2)) /
% (Math.Pow(magrearth3, 3) * Math.Pow(magrsat3, 3) * (magrearth3 + magrsat3));
% athirdbody1(1) = mu3 * (rsat3(1) * q - reci(1) / Math.Pow(magrearth3, 3));
% athirdbody1(2) = mu3 * (rsat3(2) * q - reci(2) / Math.Pow(magrearth3, 3));
% athirdbody1(3) = mu3 * (rsat3(3) * q - reci(3) / Math.Pow(magrearth3, 3));
% fprintf(fida,' a3bodyS1 eci' + athirdbody1(1), athirdbody1(2), athirdbody1(3));
% fprintf(fida,'ansr sun        0.0000000001820	0.0000000001620	-0.0000000001800\n');
% a3body(1) = athirdbody1(1);
% a3body(2) = athirdbody1(2);
% a3body(3) = athirdbody1(3);
% 
% % moon
% mu3 = mumoon;
% rsat3(1) = rmoon(1) - reci(1);
% rsat3(2) = rmoon(2) - reci(2);
% rsat3(3) = rmoon(3) - reci(3);
% magrsat3 = mag(rsat3);
% rearth3(1) = rmoon(1);
% rearth3(2) = rmoon(2);
% rearth3(3) = rmoon(3);
% magrearth3 = mag(rearth3);
% athirdbody(1) = mu3 * (rsat3(1) / Math.Pow(magrsat3, 3) - rearth3(1) / Math.Pow(magrearth3, 3));
% athirdbody(2) = mu3 * (rsat3(2) / Math.Pow(magrsat3, 3) - rearth3(2) / Math.Pow(magrearth3, 3));
% athirdbody(3) = mu3 * (rsat3(3) / Math.Pow(magrsat3, 3) - rearth3(3) / Math.Pow(magrearth3, 3));
% fprintf(fida,' a3bodyM  eci ' + athirdbody(1), athirdbody(2), athirdbody(3));
% athirdbody2(1) = -mu3 / Math.Pow(magrearth3, 3) * (rearth3(1) - 3.0 * rearth3(1) * (dot(reci, rearth3) / Math.Pow(magrearth3, 2))
% - 7.5 * rearth3(1) * Math.Pow(((dot(reci, rearth3) / Math.Pow(magrearth3, 2))), 2));
% athirdbody2(2) = -mu3 / Math.Pow(magrearth3, 3) * (rearth3(2) - 3.0 * rearth3(2) * (dot(reci, rearth3) / Math.Pow(magrearth3, 2))
% - 7.5 * rearth3(2) * Math.Pow(((dot(reci, rearth3) / Math.Pow(magrearth3, 2))), 2));
% athirdbody2(3) = -mu3 / Math.Pow(magrearth3, 3) * (rearth3(3) - 3.0 * rearth3(3) * (dot(reci, rearth3) / Math.Pow(magrearth3, 2))
% - 7.5 * rearth3(3) * Math.Pow(((dot(reci, rearth3) / Math.Pow(magrearth3, 2))), 2));
% fprintf(fida,' a3bodyM2 eci' + athirdbody2(1), athirdbody2(2), athirdbody2(3));
% q = (Math.Pow(mag(reci), 2) + 2.0 * dot(reci, rsat3)) *
% (Math.Pow(magrearth3, 2) + magrearth3 * magrsat3 + Math.Pow(magrsat3, 2)) /
% (Math.Pow(magrearth3, 3) * Math.Pow(magrsat3, 3) * (magrearth3 + magrsat3));
% athirdbody1(1) = mu3 * (rsat3(1) * q - reci(1) / Math.Pow(magrearth3, 3));
% athirdbody1(2) = mu3 * (rsat3(2) * q - reci(2) / Math.Pow(magrearth3, 3));
% athirdbody1(3) = mu3 * (rsat3(3) * q - reci(3) / Math.Pow(magrearth3, 3));
% fprintf(fida,' a3bodyM1 eci' + athirdbody1(1), athirdbody1(2), athirdbody1(3));
% fprintf(fida,'ansr moon        0.0000000000860	-0.0000000004210	-0.0000000006980\n');
% a3body(1) = a3body(1) + athirdbody1(1);
% a3body(2) = a3body(2) + athirdbody1(2);
% a3body(3) = a3body(3) + athirdbody1(3);
% fprintf(fida,'ansr sun/moon    0.0000000002730	-0.0000000002680	-0.0000000008800\n');
% 
% 
% fprintf(fida,' ------------------ find srp acceleration\n');
%  psrp = 4.56e-6;  % N/m2 = kgm/s2 / m2 = kg/ms2
% rsatsun(1) = rsun(1) - reci(1);
% rsatsun(2) = rsun(2) - reci(2);
% rsatsun(3) = rsun(3) - reci(3);
%  magrsatsun = mag(rsatsun);
% %           kg/ms2      m2      kg      km            km
% asrp(1) = -(psrp * cr * area / mass * rsatsun(1) / magrsatsun) / 1000.0;  % result in km/s
% asrp(2) = -(psrp * cr * area / mass * rsatsun(2) / magrsatsun) / 1000.0;
% asrp(3) = -(psrp * cr * area / mass * rsatsun(3) / magrsatsun) / 1000.0;
% fprintf(fida,' asrp eci ' + asrp(1), asrp(2), asrp(3));
% fprintf(fida,'ansr srp        -0.0000000001970	0.0000000001150	0.0000000000480\n');
% 
% fprintf(fida,' ------------------ add perturbing accelerations\n');
% aecef(1) = adrag(1);  % plus gravity xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% aecef(2) = adrag(2);
% aecef(3) = adrag(3);
% fprintf(fida,' aecef ' + aecef(1), aecef(2), aecef(3));
% 
% % ---- move acceleration from earth fixed coordinates to eci
% % there are no cross products here as unital
% aeci = matvecmult(transrecef2eci, aecef, 3);
% fprintf(fida,' aeci ' + aeci(1), aeci(2), aeci(3));
% 
% % find two body component of eci acceleration
% aeci2(1) = -mu * reci(1) / (Math.Pow(mag(reci), 3));
% aeci2(2) = -mu * reci(2) / (Math.Pow(mag(reci), 3));
% aeci2(3) = -mu * reci(3) / (Math.Pow(mag(reci), 3));
% fprintf(fida,' aeci2body ' + aeci2(1), aeci2(2), aeci2(3));
% 
% % totla acceleration
% aeci(1) = aeci2(1) + a3body(1) + asrp(1) + aeci(1);
% aeci(2) = aeci2(2) + a3body(2) + asrp(2) + aeci(2);
% aeci(3) = aeci2(3) + a3body(3) + asrp(3) + aeci(3);
% fprintf(fida,'total aeci ' + aeci(1), aeci(2), aeci(3));
% 



% % ------------------------------------------- timing comparisons
% fprintf(fida,'\n ===================================== Timing Comparisons =====================================');
% % timing of routines
% var watch = System.Diagnostics.Stopwatch.StartNew(fid);
% 
% for (i = 0; i < 500; i++)
% 
% straccum = '';
% order = 50;
% % unitalized calcs, show
% FullGeopM(recef, order, 'y', convArr, gravData, out aPertM, 'n', out straccum);
% end
% %  stop timer
% watch.Stop(fid);
% var elapsedMs = watch.ElapsedMilliseconds;
% fprintf(fida,'Done with Montenbruck calcs ' + (watch.ElapsedMilliseconds * 0.001), ' second  ');
% 
% watch = System.Diagnostics.Stopwatch.StartNew(fid);
% for (i = 0; i < 500; i++)
% 
% straccum = '';
% order = 50;
% % unitalized calcs, show
% FullGeopG(recef, order, 'y', convArr, unitArr, gravData, out aPertG, 'n', out straccum);
% end
% %  stop timer
% watch.Stop(fid);
% elapsedMs = watch.ElapsedMilliseconds;
% fprintf(fida,'Done with GTDS calcs ' + (watch.ElapsedMilliseconds * 0.001), ' second  ');
% 
% 
% watch = System.Diagnostics.Stopwatch.StartNew(fid);
% for (i = 0; i < 500; i++)
% 
% straccum = '';
% order = 100;
% % unitalized calcs, show
% % GTDS version
% LegPolyG(latgc, order, 'y', convArr, unitArr, out LegArrGU, out LegArrGN);
% end
% %  stop timer
% watch.Stop(fid);
% elapsedMs = watch.ElapsedMilliseconds;
% fprintf(fida,'Done with GTDS ALF calcs ' + (watch.ElapsedMilliseconds * 0.001), ' second  ');
% 
% 
% watch = System.Diagnostics.Stopwatch.StartNew(fid);
% for (i = 0; i < 500; i++)
% 
% straccum = '';
% order = 100;
% % unitalized calcs, show
% % Gottlieb version
% LegPolyGot(latgc, order, 'y', convArr, unitArr, out LegArrGotU, out LegArrGotN);
% end
% %  stop timer
% watch.Stop(fid);
% elapsedMs = watch.ElapsedMilliseconds;
% fprintf(fida,'Done with Gott ALF calcs ' + (watch.ElapsedMilliseconds * 0.001), ' second  ');
% 
% watch = System.Diagnostics.Stopwatch.StartNew(fid);
% for (i = 0; i < 500; i++)
% 
% straccum = '';
% order = 100;
% % unitalized calcs, show
% % Montenbruck version
% LegPolyM(latgc, order, 'y', convArr, out LegArrMU, out LegArrMN);
% end
% %  stop timer
% watch.Stop(fid);
% elapsedMs = watch.ElapsedMilliseconds;
% fprintf(fida,'Done with Mont ALF calcs ' + (watch.ElapsedMilliseconds * 0.001), ' second  ');
% 
% 
% watch = System.Diagnostics.Stopwatch.StartNew(fid);
% for (i = 0; i < 500; i++)
% 
% straccum = '';
% order = 100;
% % unitalized calcs, show
% % Fukushima version
% LegPolyF(latgc, order, 'y', unitArr, out LegArrF);
% end
% %  stop timer
% watch.Stop(fid);
% elapsedMs = watch.ElapsedMilliseconds;
% fprintf(fida,'Done with Fukushima ALF calcs ' + (watch.ElapsedMilliseconds * 0.001), ' second  ');
% 


% ------------------------------------------- pole test case comparisons
% fprintf(fida,'\n ===================================== Pole Test Comparisons =====================================');
% 
% rad = 180.0 / pi;
% for i = 0:500
% 
% lon = 154.0 / rad;
% latgc = (89.9 + (i / 1000.0)) / rad;
%  magr = 7378.382745;
% 
% recef(1) = (magr * cos(latgc) * cos(lon));
% recef(2) = (magr * cos(latgc) * sin(lon));
% recef(3) = (magr * sin(latgc));
% 
% straccum = '';
% order = 50;
% % unitalized calcs, show
% [straccum] = FullGeopG(recef, order, 'y', convArr, unitArr, gravData, out aPertG, 'n');
% 
% fprintf(fida,'test pole ' + (latgc * rad), (lon * rad), aPertM(1), '     ' + aPertM(2), '     ' + aPertM(3));
% end


% available files:
% GEM10Bununit36.grv
% GEMT1unit36.grv
% GEMT2unit36.grv
% GEMT3unit36.grv
% GEMT3unit50.grv
% EGM-96unit70.grv
% EGM-96unit254.grv
% EGM-08unit100.grv
% EGM-96unit70.grv
% EGM-96unit70.grv
% EGM-96unit70.grv
% GGM01Cunit90.grv
% GGM02Cunit90.grv
% GGM03Cunit100.grv
% JGM2unit70.grv
% JGM3unit70.grv
% WGS-84_EGM96unit70.grv
% WGS-84Eunit180.grv
% WGS-84unit70.grv
%
%string fname = 'D:\Dataorig\Gravity\EGM96A.TXT';       % unit
%string fname = 'D:\Dataorig\Gravity\egm2008_gfc.txt';  % unit

% % --------------------gottlieb 1993 test
% fprintf(fida,'===================================== Gottlieb 1993 test case ===================================== ');
% fprintf(fida,'GEM-10B ununitalized 36x36 ');
% % get past text in each file
% %if (fname.Contains('GEM'))    % GEM10bununit36.grv, GEMT3unit50.grv
% %    startKtr = 17;
% %if (fname.Contains('EGM-96')) % EGM-96unit70.grv
% %    startKtr = 73;
% %if (fname.Contains('EGM-08')) % EGM-08unit100.grv
% %    startKtr = 83;  % or 21 for the larger file... which has gfc in the first col too
% fname = 'D:/Dataorig/Gravity/GEM10Bununit36.grv';
% unital = 'n';
% % latgc;
% %Int32 degree, order;
% %[,] LegArr;  % montenbruck
% %[,] LegArrN;
% %[,] LegArrG;  % gtds
% %[,] LegArrGN;
% %%  [,] LegArrEx;
% %[,] LegArr1;  % geodyn
% 
% %AstroLib.gravityModelData gravData;
% 
% recef = [ 5489.1500, 802.2220, 3140.9160 ];  % km
% fprintf(fida,'recef = ' + recef(1), recef(2), recef(3));
% % these are from the vector
% latgc = Math.Asin(recef(3) / mag(recef));
%  templ = sqrt(recef(1) * recef(1) + recef(2) * recef(2));
%  rtasc;
% if (abs(templ) < 0.0000001)
% rtasc = Math.Sign(recef(3)) * pi * 0.5;
% else
% rtasc = atan2(recef(2), recef(1));
% lon = rtasc;
% fprintf(fida,'latgc lon ' + (latgc * rad), (lon * rad));
% 
% 
% [gravarr] = readgravityfield(fname, normal);
% fprintf(fida,'\ncoefficents --------------- ');
% fprintf(fida,'c  2  0  ' + gravData.c(3, 1), ' s ' + gravData.s(3, 1));
% fprintf(fida,'c  4  0  ' + gravData.c(5, 1), ' s ' + gravData.s(5, 1));
% fprintf(fida,'c  4  4  ' + gravData.c(5, 5), ' s ' + gravData.s(5, 5));
% fprintf(fida,'c 21  1 ' + gravData.c[21, 1], ' s ' + gravData.s[21, 1]);
% fprintf(fida,'\nunitalized coefficents --------------- ');
% fprintf(fida,'c  2  0  ' + gravData.cNor(3, 1), ' s ' + gravData.sNor(3, 1));
% fprintf(fida,'c  4  0  ' + gravData.cNor(5, 1), ' s ' + gravData.sNor(5, 1));
% fprintf(fida,'c  4  4  ' + gravData.cNor(5, 5), ' s ' + gravData.sNor(5, 5));
% fprintf(fida,'c 21  1 ' + gravData.cNor[21, 1], ' s ' + gravData.sNor[21, 1]);
% 
% 
% degree = 36;  % 36
% order = 36;
% LegPolyG(latgc, order, unitalized, convArr, unitArr, out LegArrGU, out LegArrGN);
% LegPolyM(latgc, order, unitalized, convArr, out LegArrMU, out LegArrMN);
% % get geodyn version
% geodynlegp(latgc, degree, order, out LegArrOU, out LegArrON);
% % get exact values
% % LegPolyEx(latgc, order, out LegArrEx);
% 
% errstr = ' ';
% sumdr1 = 0.0;
% sumdr2 = 0.0;
% fprintf(fida,'\nLegendre polynomials --------------- ');
% for (L = 1; L <= 6; L++)  % order xxxxxxxxxxxxxxxxxx
% 
% string tempstr1 = 'MN ';  % montenbruck
% string tempstr2 = 'GN ';  % gtds
% string tempstr3 = 'MU ';
% string tempstr4 = 'OU ';  % geodyn
% for (m = 0; m <= L; m++)
% 
% tempstr1 = tempstr1, L, '  ' + m, '   ' + LegArrMN[L, m];
% tempstr2 = tempstr2, L, '  ' + m, '   ' + LegArrGN[L, m];
% tempstr3 = tempstr3, L, '  ' + m, '   ' + LegArrMU[L, m];
% tempstr4 = tempstr4, L, '  ' + m, '   ' + LegArrOU[L + 1, m + 1];
% %dr1 = 100.0 * (LegArr[L, m] - LegArrEx[L, m]) / LegArrEx[L, m];
% %dr2 = 100.0 * (LegArr1[L, m] - LegArrEx[L, m]) / LegArrEx[L, m];
% %sumdr1 = sumdr1 + dr1;
% %sumdr2 = sumdr2 + dr2;
% %errstr = errstr + '\n' + L, '  ' + m, '   ' + dr1
% %   , dr2;
% end
% fprintf(fida,tempstr1);
% fprintf(fida,tempstr2);
% %  fprintf(fid,tempstr3);
% fprintf(fida,tempstr4 + '\n');
% end
% fprintf(fida,'totals gtds ' + sumdr1, ' montenbruck ' + sumdr2);
% strbuildplot.AppendLine(errstr);
% 
% fprintf(fida,'\naccelerations --------------- ');
% jdutc = 2451573.0;
% jdF = 0.1;
% straccum = '';
% order = 4;
% % unitalized calcs, show
% FullGeopM(recef, order, 'y', convArr, gravData, out aPertM, 'y', out straccum);
% % add in two body term since full geop is only disturbing part
% jdut1 = jdutc + jdF;
% %eci_ecef(ref reci, ref veci, iau80arr, MathTimeLib.Edirection.efrom, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps, AstroLib.EOpt.e80, ref recef, ref vecef);
% % time is not given, so let ecef and eci be =
% reci(1) = recef(1);
% reci(2) = recef(2);
% reci(3) = recef(3);
% 
% aeci2(1) = -398600.47 * reci(1) / (Math.Pow(mag(reci), 3));
% aeci2(2) = -398600.47 * reci(2) / (Math.Pow(mag(reci), 3));
% aeci2(3) = -398600.47 * reci(3) / (Math.Pow(mag(reci), 3));
% %aPertG(1) = aPertG(1) + aeci2(1);
% %aPertG(2) = aPertG(2) + aeci2(2);
% %aPertG(3) = aPertG(3) + aeci2(3);
% aPertM(1) = aPertM(1) + aeci2(1);
% aPertM(2) = aPertM(2) + aeci2(2);
% aPertM(3) = aPertM(3) + aeci2(3);
% 
% 
% fprintf(fida,straccum);
% %   fprintf(fid,'apertG 4 4   ' + aPertG(1), '     ' + aPertG(2), '     ' + aPertG(3));
% fprintf(fida,'apertM 4 4   ' + aPertM(1), '     ' + aPertM(2), '     ' + aPertM(3));
% fprintf(fida,'ansr          -0.00844269212018857E+00 -0.00123393633785485E+00 -0.00484659352346614E+00  km/s2\n');
% 
% straccum = '';
% order = 5;
% % unitalized calcs, show
% FullGeopG(recef, order, 'y', convArr, unitArr, gravData, out aPertG, 'y', out straccum);
% fprintf(fida,straccum);
% fprintf(fida,'apertG 5 5   ' + aPertG(1), '     ' + aPertG(2), '     ' + aPertG(3));
% % fprintf(fid,'apertM 5 5   ' + aPertM(1), '     ' + aPertM(2), '     ' + aPertM(3));
% fprintf(fida,'ansr          -0.00844260633555472E+00 -0.00123393243051834E+00 -0.00484652486332608E+00  km/s2\n');
% 
% 
% 
% % --------------------fonte 1993 test
% % fprintf(fid,'\n ===================================== Fonte 1993 test case =====================================');
% % fprintf(fid,'GEM-10B ununitalized 36x36 ');
% % fname = 'D:/Dataorig/Gravity/GEM10Bununit36.grv';
% % unital = 'n';
% % recef = [ 180.295260378399, -1145.13224944286, -6990.09446227757 ]; % km
% % fprintf(fid,'recef = ' + recef(1), recef(2), recef(3));
% % latgc = -1.40645188850273;
% % lon = -4.09449590512370;
% % fprintf(fid,'latgc lon ' + (latgc * rad), (lon * rad));
% 
% 
% % % Un-unitalized Polynomial Validation GEM10B
% % % GTDS vs Lundberg Truth GTDS (21x21 GEM10B)
% % fprintf(fid,'\ncoefficients --------------- ');
% % [gravarr] = readgravityfield(fname, normal);
% % fprintf(fid,'c  4  0    ' + gravData.c(5, 1), ' s ' + gravData.s(5, 1));
% % fprintf(fid,'c 21  0   ' + gravData.c[21, 0], ' s ' + gravData.s[21, 0]);
% % fprintf(fid,'c 21  5    ' + gravData.c[21, 5], ' s ' + gravData.s[21, 5]);
% % fprintf(fid,'c 21 20   ' + gravData.c[21, 20], ' s ' + gravData.s[21, 20]);
% % fprintf(fid,'c 21 21    ' + gravData.c[21, 21], ' s ' + gravData.s[21, 21]);
% 
% % % GTDS Emulation vs Lundberg Truth (21x21 GEM10B)
% % degree = 21;
% % order = 21;
% % LegPoly(latgc, order, out LegArr, out LegArrG, out LegArrN, out LegArrGN);
% % % get geodyn version
% % geodynlegp(latgc, degree, order, out LegArr1);
% % % get exact values
% % %   LegPolyEx(latgc, order, out LegArrEx);
% 
% % dr1 = 0.0;
% % dr2 = 0.0;
% % sumdr1 = 0.0;
% % sumdr2 = 0.0;
% % fprintf(fid,'\nLegendre polynomials --------------- ');
% % for (L = 1; L <= 6; L++)  % order
% %
% %     string tempstr1 = 'M ';
% %     string tempstr2 = 'G ';
% %     string tempstr3 = 'E ';
% %     string tempstr4 = 'O ';
% %     for (m = 0; m <= L; m++)
% %
% %         tempstr1 = tempstr1, L, '  ' + m, '   ' + LegArrN[L, m];
% %         tempstr2 = tempstr2, L, '  ' + m, '   ' + LegArrGN[L, m];
% %         % tempstr3 = tempstr3, L, '  ' + m, '   ' + LegArrEx[L, m];
% %         tempstr4 = tempstr4, L, '  ' + m, '   ' + LegArr1[L + 1, m + 1];
% %         %    dr1 = 100.0 * (LegArr[L, m] - LegArrEx[L, m]) / LegArrEx[L, m];
% %         %    dr2 = 100.0 * (LegArr1[L, m] - LegArrEx[L, m]) / LegArrEx[L, m];
% %         %sumdr1 = sumdr1 + dr1;
% %         %sumdr2 = sumdr2 + dr2;
% %         %errstr = errstr + '\n' + L, '  ' + m, '   ' + dr1
% %         %   , dr2;
% %     end
% %     fprintf(fid,tempstr1);
% %     fprintf(fid,tempstr2);
% %    % fprintf(fid,tempstr3);
% %     fprintf(fid,tempstr4 + '\n');
% % end
% %% fprintf(fid,'totals gtds ' + sumdr1, ' montenbruck ' + sumdr2);
% % fprintf(fid,'ansr 21  0 0.385389365005720                                                                      21  5   354542.107743601  354542.1077435970657340');
% % fprintf(fid,'ansr 21 20         -2442182686.11423  -2442182686.11409981594');
% % fprintf(fid,'ansr 21 21          405012060.632803  405012060.6327805324689' + '\n');
%
% % fprintf(fid,'\naccelerations --------------- ');
% % FullGeop(recef, jd, jdF, order, gravData, out aPert, out aPert1);
%
% % fprintf(fid,'apertG 21 21   ' + aPert(1), aPert(2), aPert(3));
% % fprintf(fid,'apertM 21 21   ' + aPert1(1), '     ' + aPert1(2), '     ' + aPert1(3));
% % fprintf(fid,'ansr             8.653210294968294E-7  -6.515584998975128E-6  -1.931032474628621E-5 ');
% % fprintf(fid,'ansr             8.653210294968E-7     -6.5155849989750E-6    -1.931032474628616E-5');
%
% % % --------------------fonte 1993 test
% % fprintf(fid,'\n===================================== Fonte 1993 test case =====================================');
% % fprintf(fid,'GEM-T3 unitalized 50x50 ');
% % fname = 'D:/Dataorig/Gravity/GEMT3unit50.grv';          % unit only released as 36x36 though...
% % unital = 'y';
%
%
% % [gravarr] = readgravityfield(fname, normal);
% % fprintf(fid,'\ncoefficients --------------- ');
% % fprintf(fid,'c  4  0   ' + gravData.c(5, 1), gravData.s(5, 1));
% % fprintf(fid,'c 21 20   ' + gravData.c[21, 20], gravData.s[21, 20]);
% % fprintf(fid,'c 50  0   ' + gravData.c[50, 0], gravData.s[50, 0]);
% % fprintf(fid,'c 50 50   ' + gravData.c[50, 50], gravData.s[50, 50]);
% % fprintf(fid,'c 50  5   ' + gravData.c[50, 5], gravData.s[50, 5]);
%
% % fprintf(fid,'\nLegendre polynomials --------------- ');
% % % GTDS Emulation vs Lundberg Truth (21x21 GEM10B)
% % degree = 50;
% % order = 50;
%
% % LegPoly(latgc, order, out LegArr, out LegArrG, out LegArrN, out LegArrGN);
% % % get geodyn version
% % geodynlegp(latgc, degree, order, out LegArr1);
% % % get exact
% % %  LegPolyEx(latgc, order, out LegArrEx);
%
% % fprintf(fid,'legarr4    0   ' + LegArrN(5, 1), LegArrN(5, 2));
%
% % fprintf(fid,'50  0          ' + LegArrN[50, 0]);
% % fprintf(fid,'50  0 alt      ' + LegArrGN[50, 0]);
% % fprintf(fid,'ansr 50  0      0.09634780379822722     9.634780379823085162E-02');
% % fprintf(fid,'50  0 geody    ' + LegArr1[50, 0], '\n');
% % %   fprintf(fid,'50  0 exact    ' + LegArrEx[50, 0], '\n');
% % %    fprintf(fid,'50  0 exact    ' + LegArrEx[50, 0], '\n');
%
% % fprintf(fid,'50 21       ' + LegArrN[50, 21]);
% % fprintf(fid,'50 21 alt   ' + LegArrGN[50, 21]);
% % fprintf(fid,'ansr 50  21  -1.443200082785759E+28  -14432000827857661203015450149.6553');
% % fprintf(fid,'50 21 geody ' + LegArr1[50, 21], '\n');
% % %   fprintf(fid,'50 21 exact  ' + LegArrEx[50, 21], '\n');
%
% % fprintf(fid,'50 49       ' + LegArrN[50, 49]);
% % fprintf(fid,'50 49 alt   ' + LegArrGN[50, 49]);
% % fprintf(fid,'ansr 50  49  -8.047341511222794E+39  -8.047341511222872818E+39');
% % fprintf(fid,'50 49 geody ' + LegArr1[50, 49], '\n');
% % % fprintf(fid,'50 49 exact ' + ex5049, '\n');
%
% % fprintf(fid,'50 50       ' + LegArrN[50, 50]);
% % fprintf(fid,'50 50 alt   ' + LegArrGN[50, 50]);
% % fprintf(fid,'ansr 50 50      1.334572710963763E+39   1.334572710963775698E+39' + '\n');
% % fprintf(fid,'50 50 geody ' + LegArr1[50, 50], '\n');
% % %fprintf(fid,'50 50 exact ' + ex550, '\n');
%
% % fprintf(fid,'\naccelerations --------------- ');
% % unitalized calcs, show
% % FullGeop(recef, order, unitalized, gravData, out aPert, out aPert1);
% % fprintf(fid,'apert 50 50   ' + aPert(1), aPert(2), aPert(3));
% % fprintf(fid,'ansr           8.683465146150188E-007    -6.519678538340073E-006   -1.931876804829165E-005');
% % fprintf(fid,'ansr           8.68346514615019361E-07   -6.51967853834008023E-06  -1.93187680482916393E-05');
%
% % recef = [ 487.0696937, -5330.5022406, 4505.7372146 ];  % m
% % vecef = [ -2.101083975, 4.624581986, 5.688300377 ];
%
% % write out results
% string directory = @'d:\codes\library\matlab\';
% File.WriteAllText(directory + 'legpoly.txt', strbuildall);
%
% File.WriteAllText(directory + 'legendreAcc.txt', strbuildplot);
%
% end
%                                        end
%                                end
%                            end
%                        end
%                    end
%                end
%            end
%     end
% end
%


function testhill(fid)
    constastro;
    dts = 1400.0; % second

    % circular orbit
    alt = 590.0;
    r(1) = re + alt;
    r(2) = 0.0;
    r(3) = 0.0;
    v(1) = 0.0;
    v(2) = sqrt(mu / mag(r));
    v(3) = 0.0;

    rh(1) = 0.0;
    rh(2) = 0.0;
    rh(3) = 0.0;
    vh(1) = -0.1;
    vh(2) = -0.04;
    vh(3) = -0.02;

    for i = 1:50
        dts = i * 60.0;  % second
        [rint, vint] = hillsr(rh, vh, alt, dts);
        fprintf(fid,'rint %15.11f %15.11f %15.11f %15.11f %15.11f %15.11f %15.11f\n', dts, rint(1), rint(2), rint(3), vint(1), vint(2), vint(3));
    end

    [vint] = hillsv(r, alt, dts);

end  % test hill


% ----------------------------------------------------------------------------
%
%                                 function printcov
%
% this function prints a covariance matrix
%
% author        : david vallado                  719 - 573 - 2600   23 may 2003
%
% revisions
%
% inputs description range / units
%   covin     - 6x6 input covariance matrix
%   covtype   - type of covariance             'cl','ct','fl','sp','eq',
%   cu        - covariance units(deg or rad)  't' or 'm'
%   anom      - anomaly                        'mean' or 'true' or 'tau'
%
% outputs       :
%
% locals        :
%
% references    :
% none
%
% [strout] = printcov(covin, covtype, cu, anom)
% ----------------------------------------------------------------------------*/

function [strout] = printcov(covin, covtype, cu, anom, fid)

    strout = '';
    % if (strcmp(anom, 'truea') == 1) || (strcmp(anom, 'meana') == 1)
    %     semi = 'a m  ';
    % else
    %     if (strcmp(anom, 'truen') == 1) || (strcmp(anom, 'meann') == 1)
    %         semi = 'n rad';
    %     end
    % end
    % 
    % if (strcmp(covtype, 'ct') == 1)
    %     append (strout,'cartesian covariance\n');
    %     append(strout,'        x  m            y m             z  m           xdot  m/s       ydot  m/s       zdot  m/s\n');
    % end
    % 
    % if (strcmp(covtype, 'cl') == 1)
    % 
    %     append(strout, 'classical covariance\n');
    %     if (cu == 'm')
    % 
    %         append(strout, '          ' + semi + '          ecc           incl rad      raan rad         argp rad        ');
    %         if (contains(anom, 'mean') == 1) % 'meana' 'meann'
    %             append(strout, 'm rad\n');
    %         else     % 'truea' 'truen'
    %             append(strout, ' nu rad\n');
    %         end
    %     else
    % 
    %         append(strout, '          ' + semi + '           ecc           incl deg      raan deg         argp deg        ');
    %         if (contains(anom, 'mean') == 1) % 'meana' 'meann'
    %             append(strout, ' m deg\n');
    %         else     % 'truea' 'truen'
    %             append(strout, ' nu deg\n');
    %         end
    %     end
    % end
    % 
    % if (strcmp(covtype, 'eq') == 1)
    %
    %     append(strout, 'equinoctial covariance\n');
    %     %            if (cu == 'm')
    %     if (contains(anom, 'mean') == 1) % 'meana' 'meann'
    %         append(strout, '         ' + semi + '           af              ag           chi             psi         meanlonM rad\n');
    %     else     % 'truea' 'truen'
    %         append(strout, '         ' + semi + '           af              ag           chi             psi         meanlonNu rad\n');
    %     end
    % end
    % if (strcmp(covtype, 'fl') == 1)
    %     append(strout, 'flight covariance\n');
    %     append(strout, '       lon  rad      latgc rad        fpa rad         az rad           r  m           v  m/s\n');
    % end
    %
    % if (strcmp(covtype, 'sp') == 1)
    %     append(strout, 'spherical covariance\n');
    %     append(strout, '      rtasc deg       decl deg        fpa deg         az deg           r  m           v  m/s\n');
    % end

    % format strings to show signs 'and' to not round off if trailing 0!!
    for i=1:6
        %    append(strout, covin(i, 1), num2str(covin(i, 2)), num2str(covin(i, 3)), num2str(covin(i, 4)), num2str(covin(i, 5)), num2str(covin(i, 6)),'\n');
    end
    fprintf(fid,'%17.10g %17.10g %17.10g %17.10g %17.10g %17.10g\n', covin);

end  % printcov


%----------------------------------------------------------------------------
%
%                                  function printdiff
%
% this function prints a covariance matrix difference
%
% author        : david vallado                  719 - 573 - 2600   23 may 2003
%
% revisions
%
% inputs description range / units
%   strin    - title
%   mat1     - 6x6 input matrix
%   mat2     - 6x6 input matrix
%
% outputs       :
%
% locals        :
%
% ----------------------------------------------------------------------------*/

function [strout] = printdiff(strin, mat1, mat2, fid)

    small = 1e-18;

    % format strings to show signs 'and' to not round off if trailing 0!!

    strout = '';
    %strout = 'diff ' + strin + '\n';
    dr = mat1-mat2;
    %for i = 1: 6
        %    for j=1:6
        %        dr(i, j) = mat1(i, j) - mat2(i, j);
        %        strout = append(strout, dr(i, 1), dr(i, 2), dr(i, 3), dr(i, 4), dr(i, 5), dr(i, 6),'\n');
%        fprintf(fid,' diff %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f \,', dr(:, i));
 %   end

    fprintf(fid,'%17.10g %17.10g %17.10g %17.10g %17.10g %17.10g\n', dr);

    %    append(strout, 'pctdiff % ',strin , ' pct over 1e-18\n');
    % fprintf(fid, '%14.4f%14.4f%14.4f%14.4f%14.4f%14.4f\n', 100.0 * ((mat1' - mat2') / mat1'));
    % fprintf(fid, 'Check consistency of both approaches tmct2cl-inv(tmcl2ct) diff pct over 1e-18\n');
    % fprintf(fid, '-------- accuracy of tm comparing ct2cl and cl2ct ---------\n');
    %tm1 = mat1';
    %tm2 = mat2';
    for i=1:6
        for j=1:6
            if (abs(dr(i, j)) < small || abs(mat1(i, j)) < small)
                diffmm(i, j) = 0.0;
            else
                diffmm(i, j) = 100.0 * (dr(i, j) / mat1(i, j));
            end
            %           strout = append(strout, diffmm(i, 1), diffmm(i, 2), diffmm(i, 3), diffmm(i, 4), diffmm(i, 5), diffmm(i, 6),'\n');
        end
    end
    
    fprintf(fid,'%17.10g %17.10g %17.10g %17.10g %17.10g %17.10g\n', diffmm);
end  % printdiff


function testcovct2rsw(fid)
    anom = 'meana';  % truea/n, meana/n
    anomflt = 'latlon'; % latlon  radec


    reci = [ -605.79221660; -5870.22951108; 3493.05319896 ];
    veci = [ -1.56825429; -3.70234891; -6.47948395 ];
    aeci = [ 0.001; 0.002; 0.003 ];

    % StringBuilder strbuild = new StringBuilder(fid);
    % strbuild.Clear(fid);

    fileLoc = 'D:\Codes\LIBRARY\DataLib\';
    [iau80arr] = iau80in(fileLoc);

    year = 2000;
    mon = 12;
    day = 15;
    hr = 16;
    minute = 58;
    second = 50.208;
    dut1 = 0.10597;
    dat = 32;
    timezone = 0;
    xp = 0.0;
    yp = 0.0;
    lod = 0.0;
    terms = 2;
    timezone = 0;
    ddpsi = 0.0;
    ddeps = 0.0;
    ddx = 0.0;
    ddy = 0.0;
    eqeterms = 2;


    [ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb, jdtdb, jdtdbfrac] ...
        = convtime ( year, mon, day, hr, minute, second, timezone, dut1, dat );

    year = 2004;
    mon  =   5;
    day  =  14;
    hr   =  10;
    minute  =  43;
    second  =   0.0;
    dut1 = -0.463326;
    dat  = 32;
    xp   =  0.0;
    yp   =  0.0;
    lod  =  0.0;
    timezone= 6;

    % -------- convtime    - convert time from utc to all the others
    %, tcg, jdtcg,jdtcgfrac, tcb, jdtcb,jdtcbfrac
    [ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb, jdtdb, jdtdbfrac] ...
        = convtime ( year, mon, day, hr, minute, second, timezone, dut1, dat );

    fprintf(fid,'ut1 %8.6f tut1 %16.12f jdut1 %18.11f\n',ut1,tut1,jdut1+jdut1frac );
    fprintf(fid,'utc %8.6f\n',utc );
    fprintf(fid,'tai %8.6f\n',tai );
    fprintf(fid,'tt  %8.6f ttt  %16.12f jdtt  %18.11f\n',tt,ttt,jdtt + jdttfrac );
    fprintf(fid,'tdb %8.6f ttdb %16.12f jdtdb %18.11f\n',tdb,ttdb,jdtdb + jdtdbfrac );


    % ---convert the eci state into the various other state formats(classical, equinoctial, etc)
    cartcov = [ ...
        100.0, 1.0e-2, 1.0e-2, 1.0e-4, 1.0e-4, 1.0e-4 ;...
        1.0e-2, 100.0,  1.0e-2, 1.0e-4,   1.0e-4,   1.0e-4;...
        1.0e-2, 1.0e-2, 100.0,  1.0e-4,   1.0e-4,   1.0e-4;...
        1.0e-4, 1.0e-4, 1.0e-4, 0.0001,   1.0e-6,   1.0e-6;...
        1.0e-4, 1.0e-4, 1.0e-4, 1.0e-6,   0.0001,   1.0e-6;...
        1.0e-4, 1.0e-4, 1.0e-4, 1.0e-6,   1.0e-6,   0.0001];
    cartstate = [ reci(1); reci(2); reci(3); veci(1); veci(2); veci(3) ];  % in km


    % test position and velocity going back
    avec = [ 0.0; 0.0; 0.0 ];

    [recef, vecef, aecef] = eci2ecef(reci, veci, aeci, iau80arr, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps );
    fprintf(fid,'==================== do the sensitivity tests\n');
    fprintf(fid,'1.  Cartesian Covariance\n');
    [strout ] = printcov(cartcov, 'ct', 'm', anom, fid);
    %%fprintf(fid,strout);

    fprintf(fid,'2.  RSW Covariance from Cartesian #1 above  -------------------\n');
    [cartcovrsw, tmct2cl] = covct2rsw(cartcov, cartstate);
    [strout ] = printcov(cartcovrsw, 'ct', 'm', anom, fid);
    %%fprintf(fid,strout);

    fprintf(fid,'2.  NTW Covariance from Cartesian #1 above  -------------------\n');
    [cartcovntw, tmct2cl] = covct2ntw(cartcov, cartstate);
    [strout ] = printcov(cartcovntw, 'ct', 'm', anom, fid);
    %fprintf(fid,strout);
end


function testcovct2ntw(fid)
    % done above
end


% test eci_ecef too
function testcovct2clmean(fid)
    anom = 'meana';  % truea/n, meana/n
    anomflt = 'latlon'; % latlon  radec

    reci = [ -605.79221660; -5870.22951108; 3493.05319896 ];
    veci = [ -1.56825429; -3.70234891; -6.47948395 ];
    aeci = [ 0.001; 0.002; 0.003 ];

    % StringBuilder strbuild = new StringBuilder(fid);
    % strbuild.Clear(fid);

    fileLoc = 'D:\Codes\LIBRARY\DataLib\';
    [iau80arr] = iau80in(fileLoc);

    year = 2000;
    mon = 12;
    day = 15;
    hr = 16;
    minute = 58;
    second = 50.208;
    dut1 = 0.10597;
    dat = 32;
    timezone = 0;
    xp = 0.0;
    yp = 0.0;
    lod = 0.0;
    terms = 2;
    timezone = 0;
    ddpsi = 0.0;
    ddeps = 0.0;
    ddx = 0.0;
    ddy = 0.0;

    [ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb, jdtdb, jdtdbfrac] ...
        = convtime ( year, mon, day, hr, minute, second, timezone, dut1, dat );

    % ---convert the eci state into the various other state formats(classical, equinoctial, etc)
    cartcov = [
        100.0, 1.0e-2, 1.0e-2, 1.0e-4, 1.0e-4, 1.0e-4 ; ...
        1.0e-2, 100.0,  1.0e-2, 1.0e-4,   1.0e-4,   1.0e-4; ...
        1.0e-2, 1.0e-2, 100.0,  1.0e-4,   1.0e-4,   1.0e-4; ...
        1.0e-4, 1.0e-4, 1.0e-4, 0.0001,   1.0e-6,   1.0e-6; ...
        1.0e-4, 1.0e-4, 1.0e-4, 1.0e-6,   0.0001,   1.0e-6; ...
        1.0e-4, 1.0e-4, 1.0e-4, 1.0e-6,   1.0e-6,   0.0001 ];
    cartstate = [ reci(1); reci(2); reci(3); veci(1); veci(2); veci(3) ];  % in km

    % --------convert to a classical orbit state
    [p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper ] = rv2coe (reci, veci);
    classstate(1) = a;   % km
    classstate(2) = ecc;
    classstate(3) = incl;
    classstate(4) = raan;
    classstate(5) = argp;
    if (contains(anom, 'mean')) % meann or meana
        classstate(6) = m;
    else  % truea or truen
        classstate(6) = nu;
    end
    % -------- convert to an equinoctial orbit state
    [a, n, af, ag, chi, psi, meanlonM, meanlonNu, fr] = rv2eq(reci, veci);
    if (strcmp(anom,'meana')) == 1 || (strcmp(anom,'truea') == 1 )
        eqstate(1) = a;  % km
    else % meann or truen
        eqstate(1) = n;
    end
    eqstate(2) = af;
    eqstate(3) = ag;
    eqstate(4) = chi;
    eqstate(5) = psi;
    if (contains(anom, 'mean')==1) %  meana or meann
        eqstate(6) = meanlonM;
    else % truea or truen
        eqstate(6) = meanlonNu;
    end
    % --------convert to a flight orbit state
[lon, latgc, rtasc, decl, fpa, az, magr, magv] = rv2flt ...
        ( reci, veci, iau80arr, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps )
    if (strcmp(anomflt, 'radec')==1)

        fltstate(1) = rtasc;
        fltstate(2) = decl;
    else
        if (strcmp(anomflt, 'latlon') == 1)

            fltstate(1) = lon;
            fltstate(2) = latgc;
        end
    end
    fltstate(3) = fpa;
    fltstate(4) = az;
    fltstate(5) = magr;  % km
    fltstate(6) = magv;

    % test position and velocity going back
    aeci = [ 0.0, 0.0, 0.0 ];

    [recef, vecef, aecef] = eci2ecef(reci', veci', aeci, iau80arr, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps );
    %vx = magv* ( -cos(lon)*sin(latgc)*cos(az)*cos(fpa) - sin(lon)*sin(az)*cos(fpa) + cos(lon)*cos(latgc)*sin(fpa) );
    %vy = magv* ( -sin(lon)*sin(latgc)*cos(az)*cos(fpa) + cos(lon)*sin(az)*cos(fpa) + sin(lon)*cos(latgc)*sin(fpa) );
    %vz = magv* (sin(latgc) * sin(fpa) + cos(latgc)*cos(az)*cos(fpa) );
    %% correct:
    %ve1 = magv* ( -cos(rtasc)*sin(decl)*cos(az)*cos(fpa) - sin(rtasc)*sin(az)*cos(fpa) + cos(rtasc)*cos(decl)*sin(fpa) ); % m/s
    % ve2 = magv* (-sin(rtasc) * sin(decl) * cos(az) * cos(fpa) + cos(rtasc) * sin(az) * cos(fpa) + sin(rtasc) * cos(decl) * sin(fpa));
    %ve3 = magv* (sin(decl) * sin(fpa) + cos(decl)*cos(az)*cos(fpa) );

    fprintf(fid,'==================== do the sensitivity tests\n');

    fprintf(fid,'1.  Cartesian Covariance\n');
    [strout ] = printcov(cartcov, 'ct', 'm', anom,fid);

    fprintf(fid,'2.  Classical Covariance from Cartesian #1 above (%s) -------------------\n', anom);

    [classcovmeana, tmct2cl] = covct2cl(cartcov, cartstate, anom);
    [strout ] = printcov(classcovmeana, 'cl', 'm', anom,fid);

    fprintf(fid,'  Cartesian Covariance from Classical #2 above\n');
    [cartcovmeanarev, tmcl2ct] = covcl2ct(classcovmeana, classstate, anom);
    [strout ] = printcov(cartcovmeanarev, 'ct', 'm', anom,fid);
    fprintf(fid,'\n');

    printdiff(' cartcov - cartcovmeanarev\n', cartcov, cartcovmeanarev);

    %coveci_ecef(ref cartcov, cartstate, MathTimeLib.Edirection.eto,  ref ecefcartcov, out tm, iau80arr,
    %            ttt, jdut1, lod, xp, yp, 2, ddpsi, ddeps, AstroLib.EOpt.e80);
    %printcov(cartcovmeanarev, 'ct', 'm', anom, out strout);
    %fprintf(fid,strout);
    fprintf(fid,'\n');

end  % testcovct2clmean


function testcovct2cltrue(fid)
    anom = 'truea';  % truea/n, meana/n
    anomflt = 'latlon'; % latlon  radec

    reci = [ -605.79221660; -5870.22951108; 3493.05319896 ];
    veci = [ -1.56825429; -3.70234891; -6.47948395 ];
    aeci = [ 0.001; 0.002; 0.003 ];

    %StringBuilder strbuild = new StringBuilder(fid);
    %strbuild.Clear(fid);

    fileLoc = 'D:\Codes\LIBRARY\DataLib\';
    [iau80arr] = iau80in(fileLoc);

    year = 2000;
    mon = 12;
    day = 15;
    hr = 16;
    minute = 58;
    second = 50.208;
    dut1 = 0.10597;
    dat = 32;
    timezone = 0;
    xp = 0.0;
    yp = 0.0;
    lod = 0.0;
    terms = 2;
    timezone = 0;
    ddpsi = 0.0;
    ddeps = 0.0;
    ddx = 0.0;
    ddy = 0.0;

    [ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb, jdtdb, jdtdbfrac] ...
        = convtime ( year, mon, day, hr, minute, second, timezone, dut1, dat );

    % ---convert the eci state into the various other state formats(classical, equinoctial, etc)
    cartcov = [ ...
        100.0, 1.0e-2, 1.0e-2, 1.0e-4, 1.0e-4, 1.0e-4; ...
        1.0e-2, 100.0,  1.0e-2, 1.0e-4,   1.0e-4,   1.0e-4; ...
        1.0e-2, 1.0e-2, 100.0,  1.0e-4,   1.0e-4,   1.0e-4; ...
        1.0e-4, 1.0e-4, 1.0e-4, 0.0001,   1.0e-6,   1.0e-6; ...
        1.0e-4, 1.0e-4, 1.0e-4, 1.0e-6,   0.0001,   1.0e-6; ...
        1.0e-4, 1.0e-4, 1.0e-4, 1.0e-6,   1.0e-6,   0.0001 ];
    cartstate = [ reci(1); reci(2); reci(3); veci(1); veci(2); veci(3) ];  % in km

    % --------convert to a classical orbit state
    [p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper ] = rv2coe (reci, veci);
    classstate(1) = a;  % in km
    classstate(2) = ecc;
    classstate(3) = incl;
    classstate(4) = raan;
    classstate(5) = argp;
    if (contains(anom, 'mean')) % meann or meana
        classstate(6) = m;
    else  % truea or truen
        classstate(6) = nu;
    end
    % -------- convert to an equinoctial orbit state
    [a, n, af, ag, chi, psi, meanlonM, meanlonNu, fr] = rv2eq(reci, veci);
    if (strcmp(anom, 'meana') || strcmp(anom, 'truea'))
        eqstate(1) = a;  % km
    else % meann or truen
        eqstate(1) = n;
    end
    eqstate(2) = af;
    eqstate(3) = ag;
    eqstate(4) = chi;
    eqstate(5) = psi;
    if (contains(anom, 'mean')) %  meana or meann
        eqstate(6) = meanlonM;
    else % truea or truen
        eqstate(6) = meanlonNu;
    end
    % --------convert to a flight orbit state
[lon, latgc, rtasc, decl, fpa, az, magr, magv] = rv2flt ...
        ( reci, veci, iau80arr, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps )
    if (anomflt.Equals('radec'))

        fltstate(1) = rtasc;
        fltstate(2) = decl;
    else
        if (anomflt.Equals('latlon'))

            fltstate(1) = lon;
            fltstate(2) = latgc;
        end
    end
    fltstate(3) = fpa;
    fltstate(4) = az;
    fltstate(5) = magr;  % in km
    fltstate(6) = magv;

    % test position and velocity going back
    avec = [ 0.0, 0.0, 0.0 ];

    [recef, vecef, aecef] = eci2ecef(reci, veci, aeci, iau80arr, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps );
    %vx = magv* ( -cos(lon)*sin(latgc)*cos(az)*cos(fpa) - sin(lon)*sin(az)*cos(fpa) + cos(lon)*cos(latgc)*sin(fpa) );
    %vy = magv* ( -sin(lon)*sin(latgc)*cos(az)*cos(fpa) + cos(lon)*sin(az)*cos(fpa) + sin(lon)*cos(latgc)*sin(fpa) );
    %vz = magv* (sin(latgc) * sin(fpa) + cos(latgc)*cos(az)*cos(fpa) );
    %% correct:
    %ve1 = magv* ( -cos(rtasc)*sin(decl)*cos(az)*cos(fpa) - sin(rtasc)*sin(az)*cos(fpa) + cos(rtasc)*cos(decl)*sin(fpa) ); % m/s
    % ve2 = magv* (-sin(rtasc) * sin(decl) * cos(az) * cos(fpa) + cos(rtasc) * sin(az) * cos(fpa) + sin(rtasc) * cos(decl) * sin(fpa));
    %ve3 = magv* (sin(decl) * sin(fpa) + cos(decl)*cos(az)*cos(fpa) );

    fprintf(fid,'==================== do the sensitivity tests\n');

    fprintf(fid,'1.  Cartesian Covariance\n');
    [strout ] = printcov(cartcov, 'ct', 'm', anom,fid);
    %fprintf(fid,strout);

    fprintf(fid,'2.  Classical Covariance from Cartesian #1 above (%s) -------------------\n', anom);

    [classcovtruea, tmct2cl] = covct2cl(cartcov, cartstate, anom);
    [strout ] = printcov(classcovtruea, 'cl', 'm', anom,fid);
    %fprintf(fid,strout);

    fprintf(fid,'  Cartesian Covariance from Classical #2 above\n');
    [cartcovtruearev, tmcl2ct] = covcl2ct(classcovtruea, classstate, anom);
    [strout ] = printcov(cartcovtruearev, 'ct', 'm', anom,fid);
    %fprintf(fid,strout);
    fprintf(fid,'\n');

    printdiff(' cartcov - cartcovtruearev\n', cartcov, cartcovtruearev);
    %fprintf(fid,strout);
end  % testcovct2cltrue



function testcovcl2eq(anom, fid)
    anomflt = 'latlon'; % latlon  radec

    reci = [ -605.79221660; -5870.22951108; 3493.05319896 ];
    veci = [ -1.56825429; -3.70234891; -6.47948395 ];
    aeci = [ 0.001; 0.002; 0.003 ];

    % StringBuilder strbuild = new StringBuilder(fid);
    % strbuild.Clear(fid);

    fileLoc = 'D:\Codes\LIBRARY\DataLib\';
    [iau80arr] = iau80in(fileLoc);

    year = 2000;
    mon = 12;
    day = 15;
    hr = 16;
    minute = 58;
    second = 50.208;
    dut1 = 0.10597;
    dat = 32;
    timezone = 0;
    xp = 0.0;
    yp = 0.0;
    lod = 0.0;
    terms = 2;
    timezone = 0;
    ddpsi = 0.0;
    ddeps = 0.0;
    ddx = 0.0;
    ddy = 0.0;
    eqeterms = 2;

    [ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb, jdtdb, jdtdbfrac] ...
        = convtime ( year, mon, day, hr, minute, second, timezone, dut1, dat );

    % ---convert the eci state into the various other state formats(classical, equinoctial, etc)
    cartcov = [ ...
        100.0, 1.0e-2, 1.0e-2, 1.0e-4, 1.0e-4, 1.0e-4; ...
        1.0e-2, 100.0,  1.0e-2, 1.0e-4,   1.0e-4,   1.0e-4; ...
        1.0e-2, 1.0e-2, 100.0,  1.0e-4,   1.0e-4,   1.0e-4; ...
        1.0e-4, 1.0e-4, 1.0e-4, 0.0001,   1.0e-6,   1.0e-6; ...
        1.0e-4, 1.0e-4, 1.0e-4, 1.0e-6,   0.0001,   1.0e-6; ...
        1.0e-4, 1.0e-4, 1.0e-4, 1.0e-6,   1.0e-6,   0.0001 ];
    cartstate = [ reci(1); reci(2); reci(3); veci(1); veci(2); veci(3) ];  % in km

    % --------convert to a classical orbit state
    [p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper ] = rv2coe (reci, veci);
    classstate(1) = a;   % km
    classstate(2) = ecc;
    classstate(3) = incl;
    classstate(4) = raan;
    classstate(5) = argp;
    if (contains(anom, 'mean')) % meann or meana
        classstate(6) = m;
    else  % truea or truen
        classstate(6) = nu;
    end

    % -------- convert to an equinoctial orbit state
    [a, n, af, ag, chi, psi, meanlonM, meanlonNu, fr] = rv2eq(reci, veci);
    if (strcmp(anom, 'meana')==1) || (strcmp(anom, 'truea')==1)
        eqstate(1) = a;  % km
    else % meann or truen
        eqstate(1) = n;
    end
    eqstate(2) = af;
    eqstate(3) = ag;
    eqstate(4) = chi;
    eqstate(5) = psi;
    if (contains(anom, 'mean')==1) %  meana or meann
        eqstate(6) = meanlonM;
    else % truea or truen
        eqstate(6) = meanlonNu;
    end
    % --------convert to a flight orbit state
   [lon, latgc, rtasc, decl, fpa, az, magr, magv] = rv2flt ...
        ( reci, veci, iau80arr, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps );

   if (strcmp(anomflt, 'radec')==1)
        fltstate(1) = rtasc;
        fltstate(2) = decl;
    else
        if (strcmp(anomflt, 'latlon')==1)

            fltstate(1) = lon;
            fltstate(2) = latgc;
        end
    end

    fltstate(3) = fpa;
    fltstate(4) = az;
    fltstate(5) = magr;  % km
    fltstate(6) = magv;

    % test position and velocity going back
    aeci = [ 0.0; 0.0; 0.0 ];

    [recef, vecef, aecef] = eci2ecef(reci, veci, aeci, iau80arr, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps );
    %vx = magv* ( -cos(lon)*sin(latgc)*cos(az)*cos(fpa) - sin(lon)*sin(az)*cos(fpa) + cos(lon)*cos(latgc)*sin(fpa) );
    %vy = magv* ( -sin(lon)*sin(latgc)*cos(az)*cos(fpa) + cos(lon)*sin(az)*cos(fpa) + sin(lon)*cos(latgc)*sin(fpa) );
    %vz = magv* (sin(latgc) * sin(fpa) + cos(latgc)*cos(az)*cos(fpa) );
    %% correct:
    %ve1 = magv* ( -cos(rtasc)*sin(decl)*cos(az)*cos(fpa) - sin(rtasc)*sin(az)*cos(fpa) + cos(rtasc)*cos(decl)*sin(fpa) ); % m/s
    % ve2 = magv* (-sin(rtasc) * sin(decl) * cos(az) * cos(fpa) + cos(rtasc) * sin(az) * cos(fpa) + sin(rtasc) * cos(decl) * sin(fpa));
    %ve3 = magv* (sin(decl) * sin(fpa) + cos(decl)*cos(az)*cos(fpa) );

    fprintf(fid,'==================== do the sensitivity tests\n');

    fprintf(fid,'1.  Cartesian Covariance\n');
    [strout] = printcov(cartcov, 'ct', 'm', anom,fid);
    %fprintf(fid,strout);

    fprintf(fid,'3.  Equinoctial Covariance from Classical (Cartesian) #1 above (%s) -------------------\n', anom);
    [classcovmeana, tmct2cl] = covct2cl(cartcov, cartstate, anom);
    [eqcovmeana, tmcl2eq] = covcl2eq(classcovmeana, classstate, anom, fr);

    [strout] =printcov(eqcovmeana, 'eq', 'm', anom,fid);
    %fprintf(fid,strout);

    fprintf(fid,'  Cartesian Covariance from Classical #3 above\n');
    [classcovmeana, tmeq2cl] = coveq2cl(eqcovmeana, eqstate, anom, fr);
    [strout] =printcov(classcovmeana, 'cl', 'm', anom,fid);
    %fprintf(fid,strout);

    [cartcovmeanarev, tmcl2ct] = covcl2ct(classcovmeana, classstate, anom);
    [strout] =printcov(cartcovmeanarev, 'ct', 'm', anom,fid);
    %fprintf(fid,strout);
    fprintf(fid,'\n');

    [strout] = printdiff(append(' cartcov - cartcov', anom, 'rev\n'), cartcov, cartcovmeanarev, fid);
    %fprintf(fid,strout);
end  % testcovcl2eq

function testcovct2eq(anom, fid)
    anomflt = 'latlon'; % latlon  radec


    reci = [ -605.79221660; -5870.22951108; 3493.05319896 ];
    veci = [ -1.56825429; -3.70234891; -6.47948395 ];
    aeci = [ 0.001; 0.002; 0.003 ];

    % StringBuilder strbuild = new StringBuilder(fid);
    % strbuild.Clear(fid);

    fileLoc = 'D:\Codes\LIBRARY\DataLib\';
    [iau80arr] = iau80in(fileLoc);

    year = 2000;
    mon = 12;
    day = 15;
    hr = 16;
    minute = 58;
    second = 50.208;
    dut1 = 0.10597;
    dat = 32;
    timezone = 0;
    xp = 0.0;
    yp = 0.0;
    lod = 0.0;
    terms = 2;
    timezone = 0;
    ddpsi = 0.0;
    ddeps = 0.0;
    ddx = 0.0;
    ddy = 0.0;
    eqeterms = 2;

    [ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb, jdtdb, jdtdbfrac] ...
        = convtime ( year, mon, day, hr, minute, second, timezone, dut1, dat );

    % ---convert the eci state into the various other state formats(classical, equinoctial, etc)
    cartcov = [ ...
        100.0, 1.0e-2, 1.0e-2, 1.0e-4, 1.0e-4, 1.0e-4; ...
        1.0e-2, 100.0,  1.0e-2, 1.0e-4,   1.0e-4,   1.0e-4; ...
        1.0e-2, 1.0e-2, 100.0,  1.0e-4,   1.0e-4,   1.0e-4; ...
        1.0e-4, 1.0e-4, 1.0e-4, 0.0001,   1.0e-6,   1.0e-6; ...
        1.0e-4, 1.0e-4, 1.0e-4, 1.0e-6,   0.0001,   1.0e-6; ...
        1.0e-4, 1.0e-4, 1.0e-4, 1.0e-6,   1.0e-6,   0.0001 ];
    cartstate = [ reci(1); reci(2); reci(3); veci(1); veci(2); veci(3) ];  % in km

    % --------convert to a classical orbit state
    [p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper ] = rv2coe (reci, veci);
    classstate(1) = a;   % km
    classstate(2) = ecc;
    classstate(3) = incl;
    classstate(4) = raan;
    classstate(5) = argp;
    if (contains(anom, 'mean')==1) % meann or meana
        classstate(6) = m;
    else  % truea or truen
        classstate(6) = nu;
    end
    % -------- convert to an equinoctial orbit state
    [a, n, af, ag, chi, psi, meanlonM, meanlonNu, fr] = rv2eq(reci, veci);
    if (strcmp(anom, 'meana')==1) || (strcmp(anom, 'truea')==1)
        eqstate(1) = a;  % km
    else % meann or truen
        eqstate(1) = n;
    end
    eqstate(2) = af;
    eqstate(3) = ag;
    eqstate(4) = chi;
    eqstate(5) = psi;
    if (contains(anom, 'mean')==1) %  meana or meann
        eqstate(6) = meanlonM;
    else % truea or truen
        eqstate(6) = meanlonNu;
    end

    % --------convert to a flight orbit state
[lon, latgc, rtasc, decl, fpa, az, magr, magv] = rv2flt ...
        ( reci, veci, iau80arr, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps );
    if (strcmp(anomflt, 'radec')==1)

        fltstate(1) = rtasc;
        fltstate(2) = decl;
    else
        if (strcmp(anomflt, 'latlon')==1)

            fltstate(1) = lon;
            fltstate(2) = latgc;
        end
    end
    fltstate(3) = fpa;
    fltstate(4) = az;
    fltstate(5) = magr;  % km
    fltstate(6) = magv;

    % test position and velocity going back
    avec = [ 0.0; 0.0; 0.0 ];

    [recef, vecef, aecef] = eci2ecef(reci, veci, aeci, iau80arr, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps );
    %vx = magv* ( -cos(lon)*sin(latgc)*cos(az)*cos(fpa) - sin(lon)*sin(az)*cos(fpa) + cos(lon)*cos(latgc)*sin(fpa) );
    %vy = magv* ( -sin(lon)*sin(latgc)*cos(az)*cos(fpa) + cos(lon)*sin(az)*cos(fpa) + sin(lon)*cos(latgc)*sin(fpa) );
    %vz = magv* (sin(latgc) * sin(fpa) + cos(latgc)*cos(az)*cos(fpa) );
    %% correct:
    %ve1 = magv* ( -cos(rtasc)*sin(decl)*cos(az)*cos(fpa) - sin(rtasc)*sin(az)*cos(fpa) + cos(rtasc)*cos(decl)*sin(fpa) ); % m/s
    % ve2 = magv* (-sin(rtasc) * sin(decl) * cos(az) * cos(fpa) + cos(rtasc) * sin(az) * cos(fpa) + sin(rtasc) * cos(decl) * sin(fpa));
    %ve3 = magv* (sin(decl) * sin(fpa) + cos(decl)*cos(az)*cos(fpa) );

    fprintf(fid,'==================== do the sensitivity tests\n');

    fprintf(fid,'1.  Cartesian Covariance\n');
    [strout] =printcov(cartcov, 'ct', 'm', anom,fid);
    %fprintf(fid,strout);

    fprintf(fid,'3.  Equinoctial Covariance from Cartesian #1 above (%s) -------------------\n', anom);
    [eqcovmeana, tmct2eq] = covct2eq(cartcov, cartstate, anom, fr);

    [strout] =printcov(eqcovmeana, 'eq', 'm', anom, fid);
    %fprintf(fid,strout);

    fprintf(fid,'  Cartesian Covariance from Classical #3 above\n');
    [cartcovmeanarev, tmeq2ct] = coveq2ct(eqcovmeana, eqstate, anom, fr);

    [strout] =printcov(cartcovmeanarev, 'ct', 'm', anom, fid);
    %fprintf(fid,strout);
    fprintf(fid,'\n');

    [strout] = printdiff(append(' cartcov - cartcov', anom, 'rev\n'), cartcov, cartcovmeanarev, fid);
    %fprintf(fid,strout);
end  % testcoveq2clmeann


function testcovct2fl(anomflt, fid)
    anom = 'meann';

    reci = [ -605.79221660; -5870.22951108; 3493.05319896 ];
    veci = [ -1.56825429; -3.70234891; -6.47948395 ];
    aeci = [ 0.001; 0.002; 0.003 ];

    % StringBuilder strbuild = new StringBuilder(fid);
    % strbuild.Clear(fid);

    fileLoc = 'D:\Codes\LIBRARY\DataLib\';
    [iau80arr] = iau80in(fileLoc);

    year = 2000;
    mon = 12;
    day = 15;
    hr = 16;
    minute = 58;
    second = 50.208;
    dut1 = 0.10597;
    dat = 32;
    timezone = 0;
    xp = 0.0;
    yp = 0.0;
    lod = 0.0;
    terms = 2;
    timezone = 0;
    ddpsi = 0.0;
    ddeps = 0.0;
    ddx = 0.0;
    ddy = 0.0;
    eqeterms = 2;

    [ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb, jdtdb, jdtdbfrac] ...
        = convtime ( year, mon, day, hr, minute, second, timezone, dut1, dat );

    % ---convert the eci state into the various other state formats(classical, equinoctial, etc)
    cartcov = [ ...
        100.0, 1.0e-2, 1.0e-2, 1.0e-4, 1.0e-4, 1.0e-4; ...
        1.0e-2, 100.0,  1.0e-2, 1.0e-4,   1.0e-4,   1.0e-4; ...
        1.0e-2, 1.0e-2, 100.0,  1.0e-4,   1.0e-4,   1.0e-4; ...
        1.0e-4, 1.0e-4, 1.0e-4, 0.0001,   1.0e-6,   1.0e-6; ...
        1.0e-4, 1.0e-4, 1.0e-4, 1.0e-6,   0.0001,   1.0e-6; ...
        1.0e-4, 1.0e-4, 1.0e-4, 1.0e-6,   1.0e-6,   0.0001 ];
    cartstate = [ reci(1); reci(2); reci(3); veci(1); veci(2); veci(3) ];  % in km

    % --------convert to a classical orbit state
    [p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper ] = rv2coe (reci, veci);
    classstate(1) = a;   % km
    classstate(2) = ecc;
    classstate(3) = incl;
    classstate(4) = raan;
    classstate(5) = argp;
    if (contains(anom, 'mean')==1) % meann or meana
        classstate(6) = m;
    else  % truea or truen
        classstate(6) = nu;
    end
    % -------- convert to an equinoctial orbit state
    [a, n, af, ag, chi, psi, meanlonM, meanlonNu, fr] = rv2eq(reci, veci);
    if (strcmp(anom, 'meana')==1) || (strcmp(anom, 'truea')==1)
        eqstate(1) = a;  % km
    else % meann or truen
        eqstate(1) = n;
    end
    eqstate(2) = af;
    eqstate(3) = ag;
    eqstate(4) = chi;
    eqstate(5) = psi;
    if (contains(anom, 'mean')==1) %  meana or meann
        eqstate(6) = meanlonM;
    else % truea or truen
        eqstate(6) = meanlonNu;
    end
    % --------convert to a flight orbit state
[lon, latgc, rtasc, decl, fpa, az, magr, magv] = rv2flt ...
        ( reci, veci, iau80arr, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps );
    if (strcmp(anomflt,'radec')==1)

        fltstate(1) = rtasc;
        fltstate(2) = decl;
    else
        if (strcmp(anomflt, 'latlon')==1)

            fltstate(1) = lon;
            fltstate(2) = latgc;
        end
    end
    fltstate(3) = fpa;
    fltstate(4) = az;
    fltstate(5) = magr;  % km
    fltstate(6) = magv;

    % test position and velocity going back
    avec = [ 0.0; 0.0; 0.0 ];

    [recef, vecef, aecef] = eci2ecef(reci, veci, aeci, iau80arr, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps );
    %vx = magv* ( -cos(lon)*sin(latgc)*cos(az)*cos(fpa) - sin(lon)*sin(az)*cos(fpa) + cos(lon)*cos(latgc)*sin(fpa) );
    %vy = magv* ( -sin(lon)*sin(latgc)*cos(az)*cos(fpa) + cos(lon)*sin(az)*cos(fpa) + sin(lon)*cos(latgc)*sin(fpa) );
    %vz = magv* (sin(latgc) * sin(fpa) + cos(latgc)*cos(az)*cos(fpa) );
    %% correct:
    %ve1 = magv* ( -cos(rtasc)*sin(decl)*cos(az)*cos(fpa) - sin(rtasc)*sin(az)*cos(fpa) + cos(rtasc)*cos(decl)*sin(fpa) ); % m/s
    % ve2 = magv* (-sin(rtasc) * sin(decl) * cos(az) * cos(fpa) + cos(rtasc) * sin(az) * cos(fpa) + sin(rtasc) * cos(decl) * sin(fpa));
    %ve3 = magv* (sin(decl) * sin(fpa) + cos(decl)*cos(az)*cos(fpa) );

    fprintf(fid,'==================== do the sensitivity tests\n');

    fprintf(fid,'1.  Cartesian Covariance\n');
    [strout] =printcov(cartcov, 'ct', 'm', anomflt, fid);
    %fprintf(fid,strout);

    fprintf(fid,'7.  Flight Covariance from Cartesian #1 above (%s) -------------------\n', anomflt);
    [fltcovmeana, tmct2fl] = covct2fl( cartcov, cartstate, anomflt, iau80arr, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps);
    
    if (strcmp(anomflt,'latlon')==1)
        [strout] =printcov(fltcovmeana, 'fl', 'm', anomflt, fid);
    else
        [strout] =printcov(fltcovmeana, 'sp', 'm', anomflt, fid);
    end
    %fprintf(fid,strout);

    fprintf(fid,'  Cartesian Covariance from Flight #7 above\n');
    [cartcovmeanarev, tmfl2ct] = covfl2ct( fltcovmeana, fltstate, anomflt, iau80arr, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps);

    [strout] = printcov(cartcovmeanarev, 'ct', 'm', anomflt, fid);
    %fprintf(fid,strout);
    fprintf(fid,'\n');

    [strout] = printdiff(append(' cartcov - cartcov', anomflt, 'rev\n'), cartcov, cartcovmeanarev, fid);
    %fprintf(fid,strout);
end  % testcovct2fl




