function [bool_recorded] = PMT_QE(arrival_times, t_d, linear_QE, paralyzable)
%Attenuates the number of received photons according to the specifications
%of the PMT/Discriminator system
% INPUTS:
%   - arrival_times: numeric array listing the arrival times of photons in
%   seconds
%   - t_d: dead time of pmt in seconds
%   - linear_QE: quantum efficiency of PMT assuming no saturation effects
%   - paralyzable: boolean, 1=paralyzable, 0=nonparalyzable
% OUTPUTS:
%   - bool_recorded: boolean array that flags the received photons that were
% actually recorded/registered/measured


%initialize tracking variables
bool_recorded = zeros(size(arrival_times));
prev_received_time = -100000;
prev_recorded_time = -100000;

%% main loop - loop over arrival times and determine which photons will be counted
for i=1:max(size(arrival_times))

    if paralyzable

        if prev_received_time + t_d < arrival_times(i)

            %poll random number to simulate QE effects. Even if meets
            %deadtime criteria, might not be recorded because QE < 1
            bool_recorded(i) = rand(1) < linear_QE;

        end

    else

        if prev_recorded_time + t_d < arrival_times(i)

            %poll random number to simulate QE effects. Even if meets
            %deadtime criteria, might not be recorded because QE < 1
            bool_recorded(i) = rand(1) < linear_QE;
            prev_recorded_time = arrival_times(i);
        end

    end

    %always update most recent received time
    prev_received_time = arrival_times(i);

end



end