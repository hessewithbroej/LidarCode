function [received_curves,detected_curves,arrival_times,detected_timestamps] = simulate_lidar(constants, altitudes, time, constituent_density, N_sent_curve, atmospheric_parameters, hardware_parameters)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%access information from inputted vectors
lambda = constants(1);
c = constants(2);
h = constants(3);
E_photon = h*c/lambda;

T_a = atmospheric_parameters(1);
T_c = atmospheric_parameters(2);
sigma_eff = atmospheric_parameters(3);
f_n = atmospheric_parameters(4);

A = hardware_parameters(1);
linear_QE = hardware_parameters(2);
t_d = hardware_parameters(3);


dz = altitudes(2)-altitudes(1);
dt = 2*dz/c;

%generate received signal curve (continuous)
N_received_curve_signal = generate_return(N_sent_curve, time, dt, altitudes, constituent_density, T_a, T_c, sigma_eff,A);

%generate noise photons from a poisson distribution with parameter=f_n*dt
noise = poissrnd(f_n*dt, size(time));
% N_total_received = N_received_curve + noise;
N_received_curve = N_received_curve_signal+noise;

%convert binned photon counts into a list of arrival times
%set upscale factor such that interpolated bin width is < 50% of the photon
%deadtime. Ensures sufficient resolution for photon assignment to reduce
%effects of stochasticity on pmt saturation
upscale_factor = 10^(ceil( log10(dt/t_d)));

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