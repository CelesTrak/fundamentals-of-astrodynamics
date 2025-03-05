% ------------------------------------------------------------------------------
%
%                           function HillsEQCM2eciNTW
%
%  this function finds the interceptor's ECI pos/vel vectors
%  given the ECI target and Equidistant Cylindrical
%  EQ1 (relative) interceptor vectors.
%
%  Routine not dependent on km or m for distance unit
%  all vectors are column vectors
%
%  author        : sal alfano         719-573-2600   13 aug 2011
%
%  [rintECI, vintECI] = HillsEQCM2eciNTW(rtgtECI, vtgtECI, rintEQCM, vintEQCM);
% ------------------------------------------------------------------------------

function [rintECI, vintECI] = HillsEQCM2eciNTW(rtgtECI, vtgtECI, rintEQCM, vintEQCM)
    constastro;
    
    %  find rotation matrix from ECI to NTW frame
    %  convert target and interceptor, compute vector magnitude and lambda
    rotECItoNTW1 = feci2NTW(rtgtECI,vtgtECI);
    rtgtNTW1 = rotECItoNTW1*rtgtECI;
    vtgtNTW1 = rotECItoNTW1*vtgtECI;
    magrtgt = norm(rtgtECI);
    lambdatgt = atan2(rtgtNTW1(2,1),rtgtNTW1(1,1));

    %  find necessary orbital elements of target at present (nu1) location
    mum = 398600.4415;  % do in km
    hvectgt = cross(rtgtNTW1,vtgtNTW1);
    ptgt = dot(hvectgt,hvectgt)/mum;
    eccvectgt = cross(vtgtNTW1,hvectgt)/mum-rtgtNTW1/magrtgt;
    ecctgt = norm(eccvectgt);
    atgt = ptgt/(1-ecctgt*ecctgt);
    if ecctgt > 0.00001
        perigeeunit = eccvectgt/ecctgt;
    else
        perigeeunit = rtgtNTW1/magrtgt;
    end
    lambdaperigee = atan2(perigeeunit(2,1),perigeeunit(1,1));
    nu1 = lambdatgt-lambdaperigee;

    %  reverse engineer position components
    %       find nu2 from orbit arc

    nu2 = find_nu_from_orbit_arc_sal(atgt, ecctgt, nu1, rintEQCM(2,1));
    lambdaint = nu2+lambdaperigee;
    %       find future position and velocity of target using nu2
    r2tgt = ptgt/(1.0+ecctgt*cos(nu2));
    Pvec = perigeeunit;
    Qvec = cross([0 0 1]',perigeeunit);
    r2vectgt = r2tgt*(cos(nu2)*Pvec+sin(nu2)*Qvec);
    v2vectgt = sqrt(mum/ptgt)*(-sin(nu2)*Pvec + (ecctgt+cos(nu2))*Qvec);

    %       rotate to future target (NTW2) frame & adjust for phiint
    rotNTW1toNTW2 = feci2NTW(r2vectgt,v2vectgt);
    rtgtNTW2 = rotNTW1toNTW2*r2vectgt;
    vtgtNTW2 = rotNTW1toNTW2*v2vectgt;

    %       find phi and unit vector for rintNTW1
    phiint=rintEQCM(3,1)/norm(rtgtNTW2);
    cosphiint=cos(phiint);
    sinphiint=sin(phiint);
    rintNTW1u(2,1)=cosphiint*sin(lambdaint);
    rintNTW1u(3,1)=sinphiint;
    rintNTW1u(1,1)=cosphiint*cos(lambdaint);

    %  reverse engineer velocity components
    rotphi=[cosphiint 0 -sinphiint; 0 1 0; sinphiint 0 cosphiint];
    rottoNTW3=rotphi*rotNTW1toNTW2;
    vintNTW3(2,1)=vintEQCM(2,1)+vtgtNTW1(2,1);
    vintNTW3(3,1)=vintEQCM(3,1);
    vintNTW3(1,1)=vintEQCM(1,1)+vtgtNTW2(1,1);
    vintNTW1=rottoNTW3'*vintNTW3;
    vintECI=rotECItoNTW1'*vintNTW1;

    % now finish finding position components
    rintNTW3u=rottoNTW3*rintNTW1u;
    rintNTW3(1,1)=rintEQCM(1,1)+rtgtNTW2(1,1);
    rintscale=rintNTW3(1,1)/rintNTW3u(1,1);
    rintNTW1=rintscale*rintNTW1u;
    rintECI=rotECItoNTW1'*rintNTW1;

end





