%Lidar return simulation via lidar eqn

close all
clearvars
clc

%fundamental constants
lambda = 374e-9; %laser wl in m
c = 299792498; %speed of light
h = 6.626e-34; %planck constant j*s
E_photon = h*c/lambda; %photon energy


%generate mock constituent density function of altitude
z_bounds = [1e+3,120e+3];
dz = 1000;
altitudes = z_bounds(1):dz:z_bounds(2); %altitude in m
constituent_density = [zeros([1,75]), linspace(0,100,10), repmat(100,[1,10]), linspace(100,1,25)]; %unitless density of constituent species
% constituent_density = [zeros(1,75), ones(1,45)];
constituent_density = constituent_density/sum(constituent_density); %density/concentration normalized such that full series sums to 1

dt = 2*dz/c; %time step for data collection. Gives altitude resolution of dz (accounting for there and back travel)


figure
plot(altitudes,constituent_density)
title("Constituent Density vs Altitude")
xlabel("Altitude (m)")
ylabel("Constituent Density (normalized)")
%
%Square laser power output
% t_start = 0; %laser turn on time (s)
% t_stop = t_start+dt; %laser turn off time (s)
% p_on = 1; %laser power in W when on
%
% %total time window where events may be happening (accounting for ToF for last photons emitted by laser)
% t_bounds = [0,t_stop+2*z_bounds(2)/c];
% time = t_bounds(1):dt:t_bounds(2); %time axis of power curve
%
% %use this for delta & square function power curve
% power_curve = [zeros(1,nearestDiv(t_start-t_bounds(1),dt)), p_on * ones(1,nearestDiv(t_stop-t_start,dt)), zeros(1,nearestDiv(t_bounds(2)+dt-t_stop,dt))];


%for gaussian-like curves
stdev = 3*dt; %standard deviation of laser power curve is 3 dt bins
cutoff_num_stds = 5; %truncate guassian after n stds

t_start = 0; %laser turn on time (s)
t_stop = t_start+(stdev/dt)*cutoff_num_stds*2*dt; %laser turn off time (s)
p_on = 1; %laser power in W when on



%total time window where events may be happening (accounting for ToF for last photons emitted by laser)
t_bounds = [dt,t_stop+2*z_bounds(2)/c+dt];
time = t_bounds(1):dt:t_bounds(2); %time axis of power curve
p_peak = 1; %peak power in W

%use this for a guassian-like power curve
power_curve = 1/(2*pi*stdev) * exp(-0.5*( ((t_start:dt:t_stop)-(stdev*cutoff_num_stds))/stdev ).^2);

%re-scale gaussian s.t. max(power_curve) = p_peak
power_curve = [power_curve/max(power_curve)*p_peak, zeros(1,size(time,2)-size(power_curve,2))];


% figure
% plot(time,power_curve);

%generate curve of N_sent photons over time
energy_curve = power_curve*dt; %total energy emitted by laser (j) at each timestep
N_sent_curve = floor(energy_curve/E_photon); %total number of photons emitted at each timestep
N_sent = sum(N_sent_curve); %total number of photons emitted during pulse

% figure
% plot(time, energy_curve);
figure
plot(time,N_sent_curve);
title("N Sent Photons vs Time")
xlabel("Time (sec)")
ylabel("N Sent Photons")

%attenuation parameters
T_a = 0.8;
T_c = 0.95;

%cross section parameters
%TODO: Add 1st-order temp and wind dependence --> will require vertical
%profiles of temp and vertical wind
sigma_eff = 1;

%receiver parameters
A = 10; %receiver area in m^2
linear_QE = 0.4; %QE in low-signal intensity region
t_d = 100e-9; %dead time in seconds
paralyzable = 1; %boolean for paralyzable detector




%% generate return signal

%return signal received by PMT (before QE, saturation etc)
N_received_curve = generate_return(N_sent_curve, time, dt, altitudes, constituent_density, T_a, T_c, sigma_eff);

%convert binned photon counts into a list of arrival times
%set upscale factor such that interpolated bin width is < 50% of the photon
%deadtime. Ensures sufficient resolution for photon assignment to reduce
%effects of stochasticity on pmt saturation
upscale_factor = 10^(1+ceil( log10(dt/t_d)));
arrival_times = (disperse_photons(N_received_curve,time,upscale_factor));

%PMT saturation effects. Boolean list of whether a photon was detected
npar_detected = PMT_QE(sort(arrival_times), t_d, linear_QE, 0);
par_detected = PMT_QE(sort(arrival_times), t_d, linear_QE, 1);

%save only the detected timestamps:
npar_detected_timestamps = arrival_times(logical(npar_detected));
par_detected_timestamps = arrival_times(logical(par_detected));

%convert detected timestamps back to counts per bin
npar_detected_binned = bin_photons(npar_detected_timestamps,time);
par_detected_binned = bin_photons(par_detected_timestamps,time);


%% begin including noise
f_n = 1000000; %avg freq of noise events in Hz
n_n = f_n*dt;
n_n = 0;

%generate noise photons from a poisson distribution with parameter=f_n*dt
noise = poissrnd(n_n, size(time));

N_total_received = N_received_curve + noise;

%% plots
f = figure;
ax1 = axes(f);
hold on
plot(ax1, time, N_received_curve/dt, 'k-')
plot(ax1, time, N_total_received/dt, 'g-')
plot(ax1, time, npar_detected_binned/dt, 'r-')
plot(ax1, time, par_detected_binned/dt, 'b-')
plot(ax1, [min(time),max(time)], [1/t_d,1/t_d], 'k--')
legend(["$f_{sig,received}$", "$f_{tot,received}$", "$f_{recorded,nonparalyzable}$","$f_{recorded,paralyzable}$", "$f_{dead}$"], 'interpreter', 'latex')
title("Frequency of Received Photons vs Time")
xlabel("Time (sec)")
ylabel("Freq. Received Photons (Hz)")

%estimate SNR
n_noise_tot = sum(noise);
n_sig_tot = sum(N_total_received);
SNR_est_simple = n_sig_tot/n_noise_tot;
text(0,max(N_total_received), sprintf("SNR: %0.3f", SNR_est_simple))





% figure
% xcorr_res = xcorr(constituent_density, power_curve);
% plot(xcorr_res)
% xlabel("Time (indexed, shifted)")
% ylabel("Cons. Density $\star$ Pow. Curve", 'Interpreter','latex')
% title("Cross correlation of power curve and density profile")
% N_sent
% N_received = sum(N_received_curve)
% ratio = N_received/N_sent




