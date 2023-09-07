function [photon_count_bins] = bin_photons(event_times,time_bins)
%bin_photons Bins a list of event times into a specified time bin vector.
%Effectively the opposite of disperse_photons.
%   Detailed explanation goes here

photon_count_bins = zeros(size(time_bins));

%bin all events
for i=1:max(size(event_times))
    
    bin = find(event_times(i)<time_bins, 1,'first');
    if isempty(bin)
        warning("Event occurring later than latest specified bin. Check provided time_bins and verify event_times are accurate")
    else
        photon_count_bins(bin) = photon_count_bins(bin)+1;
    end
    

end