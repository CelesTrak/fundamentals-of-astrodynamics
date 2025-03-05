% ------------------------------------------------------------------------------
%
%                           function Hillseci2EQCMNTW
%
%  this function finds the relative pos/vel vectors in the Equidistant
%  Cylindrical frame given the ECI target and interceptor vectors.
%
%  units are in meters and seconds
%  all vectors are column vectors
%
%  author        : sal alfano              719-573-2600   26 aug 2010
%
%  [rintEQCM, vintEQCM]  =  Hillseci2EQCMNTW(rtgtECI, vtgtECI, rintECI, vintECI);
% ------------------------------------------------------------------------------

function [rintEQCM, vintEQCM]  =  Hillseci2EQCMNTW(rtgtECI, vtgtECI, rintECI, vintECI)
    constastro;

    %  find rotation matrix from ECI to NTW1 frame for target
    %  convert target and interceptor, compute vector magnitudes
    rotECItoNTW1 = feci2NTW(rtgtECI,vtgtECI);
    rtgtNTW1 = rotECItoNTW1*rtgtECI;
    vtgtNTW1 = rotECItoNTW1*vtgtECI;
    rintNTW1 = rotECItoNTW1*rintECI;
    vintNTW1 = rotECItoNTW1*vintECI;
    magrtgt = norm(rtgtECI);
    magrint = norm(rintECI);

    %  find long/lat rotation angles (radians) to go from target to interceptor
    sinphiint = rintNTW1(3,1)/magrint;
    phiint = asin(sinphiint);
    cosphiint = cos(phiint);
    lambdatgt = atan2(rtgtNTW1(2,1),rtgtNTW1(1,1));
    lambdaint = atan2(rintNTW1(2,1),rintNTW1(1,1));

    %  find necessary orbital elements of target at present (nu1) and future (nu2) locations
    mum = 398600.4415*10^9;
    hvectgt = cross(rtgtNTW1,vtgtNTW1);
    ptgt = dot(hvectgt,hvectgt)/mum;
    eccvectgt = cross(vtgtNTW1,hvectgt)/mum-rtgtNTW1/magrtgt;
    ecctgt = norm(eccvectgt);
    atgt = ptgt/(1-ecctgt*ecctgt);
    if ecctgt > 0.00001
        perigeeunit = eccvectgt/ecctgt;
    else
        perigeeunit = rtgtNTW1/magrtgt;
    end;
    lambdaperigee = atan2(perigeeunit(2,1),perigeeunit(1,1));
    nu1 = lambdatgt-lambdaperigee;
    nu2 = lambdaint-lambdaperigee;

    %  find future position and velocity of target
    r2tgt = ptgt/(1.0 + ecctgt*cos(nu2));
    Pvec = perigeeunit;
    Qvec = cross([0 0 1]',perigeeunit);
    %         r1vectgt = magrtgt*(cos(nu1)*Pvec + sin(nu1)*Qvec); % sanity check target r&v
    %         v1vectgt = sqrt(mum/ptgt)*(-sin(nu1)*Pvec  +  (ecctgt + cos(nu1))*Qvec);
    r2vectgt = r2tgt*(cos(nu2)*Pvec + sin(nu2)*Qvec);
    v2vectgt = sqrt(mum/ptgt)*(-sin(nu2)*Pvec  +  (ecctgt + cos(nu2))*Qvec);

    %  rotate all to future target (NTW2) frame & adjust for phiint
    rotNTW1toNTW2 = feci2NTW(r2vectgt,v2vectgt);
    rtgtNTW2 = rotNTW1toNTW2*r2vectgt;
    vtgtNTW2 = rotNTW1toNTW2*v2vectgt;
    rotphi = [cosphiint 0 -sinphiint; 0 1 0; sinphiint 0 cosphiint];
    rottoNTW3 = rotphi*rotNTW1toNTW2;
    rintNTW3 = rottoNTW3*rintNTW1;
    vintNTW3 = rottoNTW3*vintNTW1;

    %  find position component positions
    rintEQCM(1,1) = rintNTW3(1,1)-rtgtNTW2(1,1);
    rintEQCM(2,1) = orbit_arc_sal(atgt,ecctgt,nu1,nu2);
    rintEQCM(3,1) = phiint*norm(rtgtNTW2);

    %  find velocity component positions
    vintEQCM(1,1) = vintNTW3(1,1)-vtgtNTW2(1,1);
    vintEQCM(2,1) = vintNTW3(2,1)-vtgtNTW1(2,1);
    vintEQCM(3,1) = vintNTW3(3,1);

end
