function [photon_count_bins] = bin_photons(event_times,time_bins)
%bin_photons Bins a list of event times into a specified time bin vector.
%Effectively the opposite of disperse_photons.
%   Assumes that time_bins is strictly increasing!

photon_count_bins = zeros(size(time_bins));

if ~any(time_bins == sort(time_bins))
    error("Time bins are not strictly increasing!!!")
end

if isempty(event_times)
    return
end

%bin all events
for i=1:max(size(event_times))
    
    bin = find(event_times(i)<time_bins, 1,'first');
    if isempty(bin)
        warning("Event occurring later than latest specified bin. Check provided time_bins and verify event_times are accurate")
    else
        photon_count_bins(bin) = photon_count_bins(bin)+1;
        if photon_count_bins(bin) > 1
            pause(0.0001)
        end
    end
    

end