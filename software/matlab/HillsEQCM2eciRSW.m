% ------------------------------------------------------------------------------
%
%                           function HillsEQCM2eciRSW
%
%  this function finds the interceptor's ECI (RTN) pos/vel vectors
%  given the ECI target and Modified Equidistant Cylindrical (EQCM)
%  interceptor vectors.
%
%  Routine IS dependent on km for distance unit
%  all vectors are column vectors
%
%  units are in km, seconds, and radians
%  all vectors are column vectors
%
%  author        : sal alfano         719-573-2600   25 sep 2012
%  
%  [rintECI, vintECI] = HillsEQCM2eciRSW(rtgtECI, vtgtECI, rintEQCM, vintEQCM );
% ------------------------------------------------------------------------------

function [rintECI, vintECI] = HillsEQCM2eciRSW(rtgtECI, vtgtECI, rintEQCM, vintEQCM )

    %  find rotation matrix from ECI to RTN frame
    %  convert target and compute position vector magnitude
    %  In RTN frame lambdatgt will be 0
    rotECItoRTN1 = fECItoRTNsal(rtgtECI,vtgtECI);
    rtgtRTN1 = rotECItoRTN1*rtgtECI;
    vtgtRTN1 = rotECItoRTN1*vtgtECI;
    magrtgt  = norm(rtgtRTN1);

    %  find necessary orbital elements of target at present (nu1) location
    mum = 398600.4415; % use km!
    hvectgt = cross(rtgtRTN1,vtgtRTN1);
    ptgt = dot(hvectgt,hvectgt)/mum;
    eccvectgt = cross(vtgtRTN1,hvectgt)/mum-rtgtRTN1/magrtgt;
    ecctgt = norm(eccvectgt);
    atgt= ptgt/(1-ecctgt*ecctgt);
    if ecctgt > 0.00001
        perigeeunit = eccvectgt/ecctgt;
    else
        perigeeunit = rtgtRTN1/magrtgt;
    end;
    lambdaperigee = atan2(perigeeunit(2,1),perigeeunit(1,1));
    nu1 = -lambdaperigee;

    %  find nu2 and lambda from orbit arc
    arclength = rintEQCM(2,1);
    [ea1, m] = newtonnu ( ecctgt,nu1 );
    if abs(arclength) > 0.001   % 1 m tolerance on the arclength
        DE = arclength / atgt;  % arclength
        Deltaea = inverselliptic2(DE,ecctgt^2 );  % this will be an eccentric anomaly
        [F1, E1] = elliptic12( ea1, ecctgt^2 );
        ea2e = ea1 + Deltaea;

        % refine range if ea0 is non-zero because Deltaea is not the same at
        % different points in the orbit
        ii = 1;
        arclength1a = arclength + 10.0;  % initial value
        while (ii < 10) && (abs(arclength1a - arclength) > 0.001)  % 1 m
            [F2, E2] = elliptic12( ea2e, ecctgt^2 );
            arclength1a = atgt*(E2 - E1);
            corr = arclength / (ea2e-ea1);   % km/rad
            ea2e = ea2e - (arclength1a-arclength)/corr;
            ii = ii + 1;
        end
        [m,nu2] = newtone ( ecctgt,ea2e );  % convert eccentric back to true
    else
        nu2 = nu1;
    end   % if arclength is nearly 0.0

    lambda = nu2 - nu1;
    sinlambda = sin(lambda);
    coslambda = cos(lambda);

    %  find future position and velocity of target using nu2
    r2tgt = ptgt / (1.0 + ecctgt*cos(nu2));
    Pvec = perigeeunit;
    Qvec = cross([0 0 1]',Pvec);
    r2vectgt = r2tgt*(cos(nu2)*Pvec+sin(nu2)*Qvec);
    v2vectgt = sqrt(mum/ptgt)*(-sin(nu2)*Pvec + (ecctgt+cos(nu2))*Qvec);

    %  rotate to future target (RTN2) frame
    rotRTN1toRTN2 = fECItoRTNsal(r2vectgt,v2vectgt);
    rtgtRTN2 = rotRTN1toRTN2*r2vectgt;
    vtgtRTN2 = rotRTN1toRTN2*v2vectgt;

    %  find phi
    phi = rintEQCM(3,1)/r2tgt;
    cosphi = cos(phi);
    sinphi = sin(phi);

    %  find interceptor position unit vector in RTN1 frame
    rintunitRTN1(3,1) = sinphi;
    rintunitRTN1(2,1) = sinlambda*cosphi;
    rintunitRTN1(1,1) = coslambda*cosphi;

    %  find SEZ conversion
    rottoSEZ = zeros(3,3);
    rottoSEZ(1,1) = sinphi*coslambda;
    rottoSEZ(1,2) = sinphi*sinlambda;
    rottoSEZ(1,3) = -cosphi;
    rottoSEZ(2,1) = -sinlambda;
    rottoSEZ(2,2) = coslambda;
    rottoSEZ(2,3) = 0.0;
    rottoSEZ(3,1) = cosphi*coslambda;
    rottoSEZ(3,2) = cosphi*sinlambda;
    rottoSEZ(3,3) = sinphi;
    rintunitSEZ = rottoSEZ*rintunitRTN1;

    %  determine proper scaling from Z component of interecptor
    rintSEZ(3,1) = rintEQCM(1,1)+rtgtRTN2(1,1);
    magrint = rintSEZ(3,1)/rintunitSEZ(3,1);
    rintRTN1 = magrint*rintunitRTN1;

    %  find velocity component positions in RTN1 frame
    lamdadot = (vintEQCM(2,1) + (vtgtRTN1(2,1)/magrtgt)*magrtgt)/r2tgt;
    vintSEZ(1,1) = (-vintEQCM(3,1)/r2tgt)*magrint;
    vintSEZ(2,1) = lamdadot*magrint*cosphi;
    vintSEZ(3,1) = vintEQCM(1,1) + vtgtRTN2(1,1);
    vintRTN1 = rottoSEZ'*vintSEZ;

    %  now rotate all into original frame
    rintECI = rotECItoRTN1'*rintRTN1;
    vintECI = rotECItoRTN1'*vintRTN1;
end
