function [N_received_curve] = generate_return(N_sent_curve, time, dt, altitudes, constituent_density, T_a, T_c, sigma_eff)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


atmospheric_absorption_factor = T_a*T_c*sigma_eff*T_c*T_a; %scaling factor to account for atmospheric absorptions

%compute N received photons at each time bin
for i=1:size(time,2)

    t = i*dt;

    %compute effects at each altitude - will be bounded by overall passage
    %of time & ToF consideratons
    for j=i:-1:1
        z = j*dz; %altitude bin - (considering factor of 2 ToF)

        %determine index of current altitude bin in "altitudes"
        [~,ind] = min(abs(altitudes-z));

        %can ignore returns from above our maximum altitude
        if z > max(altitudes)
            continue
        end

        %does not yet include effects of PMT QE
        N_received_curve(i) = N_received_curve(i) + ...
            N_sent_curve(i-j+1)* atmospheric_absorption_factor * (A/(4*pi*z^2)) * constituent_density(ind);
    end

end


end