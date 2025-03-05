% ------------------------------------------------------------------------------
%
%                           function orbit_arc_sal
%
%
% Computes the length of an arc on a closed orbit
% between true anomaly angles nu1 and nu2 (in radians).
% Arc length is determined by dividing the arc in small straight segments.
%
%  author        : sal alfano         719-573-2600   13 aug 2011
%
%  [arc] = orbit_arc_sal(a,ecc,nu1,nu2);
% ------------------------------------------------------------------------------

function [arc] = orbit_arc_sal(a,ecc,nu1,nu2)

    % determine eccentric anomalies from true
    E1 = 2*atan(sqrt((1.0-ecc)/(1.0+ecc))*tan(nu1/2));
    E2 = 2*atan(sqrt((1.0-ecc)/(1.0+ecc))*tan(nu2/2));

    % correct direction by changing E2 to be within half rev of E1
    E_diff = E2-E1;
    if E_diff < -pi
        E_diff = E_diff+2*pi;
    end
    if E_diff > pi
        E_diff = E_diff-2*pi;
    end
    E2 = E1+E_diff;

    % set angle step size (adjust ang_step as needed)
    ang_step = 0.001;
    num_steps = ceil(abs(E_diff)/ang_step)+1;
    dE = ang_step*sign(E_diff);

    % loop through and sum parts
    arc_unit = 0;
    E = E1;
    for stp = 1:num_steps
        E = E+dE;
        ecc_cos = ecc*cos(E-dE/2);
        d_arc = sqrt(1-ecc_cos*ecc_cos);
        arc_unit = arc_unit+d_arc;
        if abs(E-E1) >=  abs(E_diff)
            arc_unit = arc_unit-d_arc*(E-E1-E_diff)/dE;
            break;
        end
    end
    arc = arc_unit*a*dE;

end
