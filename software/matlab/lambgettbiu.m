%
% form the minimum time and universal variable matrix
%
% [tbiu, tbilu] = lambgettbiu(r1, r2, 3)

function [tof, kbi] = lambgettbiu( r1, r2, order )

    tof = zeros(order, 2);
    %tbi = [0 0; 0 0; 0 0; 0 0; 0 0];
    for i = 1: order
        [kbil, tofl] = lambertumins( r1, r2, i, 'S' );
        tof(i,1) = kbil;
        tof(i,2) = tofl;
    end

    kbi = zeros(order, 2);
    for i = 1: order
        [kbil, tofl] = lambertumins( r1, r2, i, 'L' );
        kbi(i,1) = kbil;
        kbi(i,2) = tofl;
    end

end

