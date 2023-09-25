function [event_times] = disperse_photons(photon_count_bins,time_bins,upscale_factor)
%disperse_photons: Converts a vector of the number of photons sent or
%received in a given time bin to a list of timestamped photon send/receive
%events.
%   This function generates a vector of timestampped photon events (either
%   sent or received), given a vector that contains the total number of
%   photons sent or received in a given bin and the start and end times of
%   each bin. First, the original vector containing the number of photons
%   arriving in a given bin is (re-sampled? not sure what the most precise
%   verb for this process is) to a higher resolution that is a factor of
%   *upscale_factor* higher than the original time resolution/bin widths.
%   This enables a more precise distribution of arriving photons within a
%   bin, although it can add significantly to total simulation runtime if
%   the *upscale_factor* is too large, 10 seems to be a good starting
%   point. The re-sampled vector is then integrated from the left, and when
%   the cummulative sum reaches an integer value, we choose a sub-bin in
%   which to generate a photon from all the sub-bins included in that
%   integration segment, weighted by each sub-bin's expected photon count
%   (which is typically 0<x<1). The photon arrival timestamp is chosen
%   randomly from the times included in the chosen sub-bin. We repeat this
%   process until we have fully integrated the original signal and
%   dispersed all the photon events.
%   INPUTS:
%       - photon_count_bins: original vector of # of photon counts per bin
%       - time_bins: vector containing the right edges of each original bin in
%            photon_count_bins. First bin's left edge is assumed to be 0
%       - upscale_factor: re-sampling resolution increase factor for more
%            precise distribution of photons within a bin. Must be integer
%            and >=1, 10 is a good starting point
%   OUTPUTS:
%       - event_times: vector of photon event timestamps

if ~mod(upscale_factor,1) == 0
    error("upscale_factor must be integer valued and >= 1")
end


%want to convert photon_count_bins into a list of photon arrival
%timestamps.
event_times = [];
prev_reached_bin = [1,1]; %need to track both bin & sub-bin of most recent arrival
running_sum = 0;
dt = mean(diff(time_bins));
N_exp = zeros(max(size(photon_count_bins)),upscale_factor);

for i=1:max(size(photon_count_bins))

    if i==1
        t_lower = 0 ;
    else
        t_lower = (i-1)*dt;
    end
    t_upper = time_bins(i);
    t_interp = t_lower:dt*(1/upscale_factor):t_upper;

    if i==1
        prev_N = 0;
    else
        prev_N = photon_count_bins(i-1);
    end

    %for debug
%     if photon_count_bins(i) >0
%         pause(0.1)
%     end

    %linear interpolation at upscaled resolution.
    for j=1:size(t_interp,2)-1
        % Number of photons expected within jth subbin
        N_exp(i,j) = (prev_N+(j*(1/upscale_factor))*(photon_count_bins(i)-prev_N))/upscale_factor;
        %we want to track the cummulative number of photons we expect to
        %see through our current bin.
        running_sum = running_sum + N_exp(i,j);

        %we should have seen a photon by this point (expected N photons > 1)
        % Choose one of the bins in the range from our last received photon to put this photon in.
        %Choice should be weighted by the expected number of photons
        %received in that bin.
        %If we happen to get multiple photons expected in a single step,
        %distribute them all in the same manner.
        while running_sum >= 1

            %need to make sure we account for any subbins from the previous
            %main bins, ensure they are included as valid bins for photon assignment
            if prev_reached_bin(1) ~= i
                %we've moved on to a new main bin since our last photon event.
                %When we choose where to assign this photon event, we need to
                %be sure to include those subbins from the previous main
                %bin(s) that contributed to the running sum. First choose the main bin
                %based on the total expected events that were contributed from
                %each main bin

                %iterate over each previous contributing main bin &
                %calculate weights
                weights = zeros(1+i-prev_reached_bin(1),1);
                for k=prev_reached_bin(1):i

                    %only the earliest and latest main bins may have a fraction of its
                    %subbins contributing
                    shifted_ind = 1+k-prev_reached_bin(1);
                    if k==1
                        weights(shifted_ind) = sum(N_exp(k,prev_reached_bin(2):end));
                        subweights = N_exp(k,prev_reached_bin(2):end);
                        subbins = prev_reached_bin(2):upscale_factor;
                    elseif k==i
                        weights(shifted_ind) = sum(N_exp(k,1:j));
                        subweights = N_exp(k,1:j);
                        subbins = 1:j;
                        %all other main bins are fully represented
                    else
                        weights(shifted_ind) = sum(N_exp(k,:));
                        subweights = N_exp(k,:);
                        subbins = 1:upscale_factor;
                    end

                end
                %choose main bin
                mainbin = randsample(prev_reached_bin(1):i, 1, true, weights);
                %if only the current main bin is represented, just use that
            else
                mainbin = i;
                subweights = N_exp(i,prev_reached_bin(2):j);
                subbins = prev_reached_bin(2):j;
            end


            %Now that we know which main bin we're using,
            %select the subbin in which we'll log photon, using the expected
            %photon count in each subbin as weight
            %Randsample acts differently if the first argument is a single
            %number instead of a vector. Only use it if we have multiple
            %subbins included, otherwise just pick the one subbin
            if prev_reached_bin(2) ~= j
                subbin = randsample(subbins,1,true,subweights);
            else
                subbin = j;
            end

            %we now know the main and sub bin number that we want to place
            %our photon in. Calculate the time bounds of this chosen bin
            %and randomly pick a time included in the bounds
            if mainbin == 1
                main_lower = 0;
            else
                main_lower = time_bins(mainbin-1);
            end
            main_upper = time_bins(mainbin);
            subbin_lower = main_lower + (dt/upscale_factor)*(subbin-1);
            subbin_upper = main_lower + (dt/upscale_factor)*(subbin);

            event_times(end+1) = subbin_lower + (subbin_upper-subbin_lower)*rand;
            running_sum = running_sum-1;
            prev_reached_bin = [i,j];

        end
    end

end





end