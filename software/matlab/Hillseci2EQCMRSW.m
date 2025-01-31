% ------------------------------------------------------------------------------
%
%                           function Hillseci2EQCMRSW
%
%  this function finds the relative pos/vel vectors in the Modified
%  Equidistant Cylindrical (EQCM) frame given the ECI target and
%  interceptor vectors with RTN (=RSW) ordered components.
%
%  Routine IS dependent on km for distance unit due to mu
%  all vectors are column vectors
%
%  units are in meters, seconds, and radians
%  all vectors are column vectors
%
%  author        : sal alfano              719-573-2600   25 sep 2012
%
%  [rintEQCM, vintEQCM] = Hillseci2EQCMRSW(rtgtECI, vtgtECI, rintECI,
%  vintECI);
% ------------------------------------------------------------------------------

function [rintEQCM, vintEQCM] = Hillseci2EQCMRSW(rtgtECI, vtgtECI, rintECI, vintECI)

    %  find rotation matrix from ECI to RTN1 frame for target
    %  convert target and interceptor, compute vector magnitudes
    rotECItoRTN1 = fECItoRTNsal(rtgtECI, vtgtECI);
    rtgtRTN1 = rotECItoRTN1*rtgtECI;
    vtgtRTN1 = rotECItoRTN1*vtgtECI;
    rintRTN1 = rotECItoRTN1*rintECI;
    vintRTN1 = rotECItoRTN1*vintECI;
    magrtgt = norm(rtgtRTN1);
    magrint = norm(rintRTN1);

    %  find lamda/phi (long/lat) rotation angles
    %  to go from target to interceptor (lambdatgt will be 0)
    sinphi = rintRTN1(3,1) / magrint;
    phi = asin(sinphi);
    cosphi = cos(phi);
    lambda = atan2(rintRTN1(2,1), rintRTN1(1,1));
    coslambda = cos(lambda);
    sinlambda = sin(lambda);

    %  find necessary orbital elements of target at present (nu1)
    %  and future (nu2) locations
    mum =  398600.4415; %  for km, else *10^9 for m;
    hvectgt =  cross(rtgtRTN1, vtgtRTN1);
    ptgt =  dot(hvectgt,hvectgt)/mum;
    eccvectgt =  cross(vtgtRTN1, hvectgt)/mum - rtgtRTN1/magrtgt;
    ecctgt =  norm(eccvectgt);
    atgt =  ptgt/(1 - ecctgt*ecctgt);
    % the check on ecc is sensitive. can result in diffs at certain
    % times
    if ecctgt > 0.0000001
        perigeeunit =  eccvectgt/ecctgt;
    else
        perigeeunit =  rtgtRTN1/magrtgt;
    end
    if ecctgt > 0.99
        rtgtRTN1
        vtgtRTN1
        rintRTN1
        vintRTN1
        dbstop;
    end
    lambdaperigee =  atan2(perigeeunit(2,1), perigeeunit(1,1));
    nu1 =  -lambdaperigee;
    nu2 =  lambda-lambdaperigee;

    %fprintf(1,'angles %11.6f %11.6f %11.6f %11.6f \n', phi*57.295, lambda*57.295, nu1*57.295, nu2*57.295);

    %  find future position and velocity of target
    r2tgt =  ptgt / (1.0 + ecctgt*cos(nu2));
    Pvec =  perigeeunit;
    Qvec =  cross([0 0 1]',Pvec);
    r2vectgt =  r2tgt*(cos(nu2)*Pvec + sin(nu2)*Qvec);
    v2vectgt =  sqrt(mum/ptgt)*(-sin(nu2)*Pvec + (ecctgt+cos(nu2))*Qvec);

    %  rotate all to future target (RTN2) frame & adjust for phi
    rotRTN1toRTN2 =  fECItoRTNsal(r2vectgt, v2vectgt);
    rtgtRTN2 =  rotRTN1toRTN2*r2vectgt;
    vtgtRTN2 =  rotRTN1toRTN2*v2vectgt;

    %  find interceptor SEZ components
    rottoSEZ =  zeros(3,3);
    rottoSEZ(1,1) = sinphi*coslambda;
    rottoSEZ(1,2) = sinphi*sinlambda;
    rottoSEZ(1,3) = -cosphi;
    rottoSEZ(2,1) = -sinlambda;
    rottoSEZ(2,2) = coslambda;
    rottoSEZ(2,3) = 0.0;
    rottoSEZ(3,1) = cosphi*coslambda;
    rottoSEZ(3,2) = cosphi*sinlambda;
    rottoSEZ(3,3) = sinphi;
    rintSEZ = rottoSEZ*rintRTN1;
    vintSEZ = rottoSEZ*vintRTN1;

    %  find position component positions
    rintEQCM(1,1) = rintSEZ(3,1) - rtgtRTN2(1,1);

    % try function for elliptic integral instead
    [ea0,m] = newtonnu ( ecctgt,nu1 );
    [ea1,m] = newtonnu ( ecctgt,nu2 );

    % fix quadrants for special cases
    if abs(ea1-ea0) > pi
        if ea0 < 0.0
            ea0 = 2.0*pi + ea0;
        else
            ea0 = 2.0*pi - ea0;
        end
    end

    [F0, E0] = elliptic12( ea0, ecctgt^2 );
    [F1, E1] = elliptic12( ea1, ecctgt^2 );  % with implied tol of 1e-22, no error. Shows up at about 1e-5

    fixit = 0.0;
    if E1-E0 < 0.0
        fixit = pi;
    end

    % arc length value
    rintEQCM(2,1) = (atgt * (E1 - E0+fixit));
    %        fprintf(1,'elliptic integral, %14.6f, %14.8f, %14.6f, %12.5f, %12.5f, ANS %12.6f  %12.6f \n', ...
    %                      atgt, ecctgt, b, m0*rad, m1*rad, arclength, rintEQCM(2,1) );

    rintEQCM(3,1) = phi*r2tgt;
    %fprintf( 1,'%f,%f, ea0, %f, ea1, %f \n',E1-E0, ea1-ea0, ea0, ea1);
    %  find velocity component positions
    lamdadot =  vintSEZ(2,1) / (magrint*cosphi);
    vintEQCM(1,1) = vintSEZ(3,1) - vtgtRTN2(1,1);
    vintEQCM(2,1) = lamdadot*r2tgt - (vtgtRTN1(2,1) / magrtgt) * magrtgt;
    vintEQCM(3,1) = (-vintSEZ(1,1) / magrint)*r2tgt;

end



