function [received_curves,detected_curves,arrival_times,detected_timestamps] = simple_return(N_sent_curve,dt,time,t_travel,t_d,linear_QE,f_n)
%UNTITLED17 Summary of this function goes here
%   Detailed explanation goes here


%shift by travel time
for i=1:size(N_sent_curve,2)
    N_received_curve_signal( mod(i-1+t_travel,size(N_sent_curve,2))+1 ) = N_sent_curve(i);
end

%normalize received signal photon counts to N_s (expected number of
%arriving signal photons in an illumination). Then, use that expected
%number to determine whether a received photon is placed in a given bin
N_s_exp = 1;
N_received_curve_signal = N_received_curve_signal/sum(N_received_curve_signal)*N_s_exp;

tmp = rand(size(N_received_curve_signal));
N_received_curve_signal = double(N_received_curve_signal > tmp);


N_received_curve = N_received_curve_signal + poissrnd(f_n*dt,size(N_received_curve_signal));




%convert binned photon counts into a list of arrival times
%set upscale factor such that interpolated bin width is < 50% of the photon
%deadtime. Ensures sufficient resolution for photon assignment to reduce
%effects of stochasticity on pmt saturation
upscale_factor = 10^(ceil( log10(dt/t_d)));
upscale_factor = 1;

arrival_times = sort(disperse_photons(N_received_curve,time,upscale_factor));

%PMT saturation effects. Boolean list of whether a photon was detected
npar_detected = PMT_QE((arrival_times), t_d, linear_QE, 0);
par_detected = PMT_QE((arrival_times), t_d, linear_QE, 1);

%save only the detected timestamps:
npar_detected_timestamps = arrival_times(logical(npar_detected));
par_detected_timestamps = arrival_times(logical(par_detected));

%convert detected timestamps back to counts per bin
npar_detected_binned = bin_photons(npar_detected_timestamps,time);
par_detected_binned = bin_photons(par_detected_timestamps,time);

received_curves = {N_received_curve, N_received_curve_signal};
detected_curves = {par_detected_binned, npar_detected_binned};

detected_timestamps = {par_detected, npar_detected};


end