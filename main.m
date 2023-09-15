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
%constituent_density = [zeros([1,75]), linspace(0,100,10), repmat(100,[1,10]), linspace(100,1,25)]; %unitless density of constituent species
%constituent_density = constituent_density/sum(constituent_density); %density/concentration normalized such that full series sums to 1

constituent_density = generate_gaussian(85000,3000, [z_bounds(1),z_bounds(2),dz], 1, 5);

dt = 2*dz/c; %time step for data collection. Gives altitude resolution of dz (accounting for there and back travel)


figure
plot(altitudes,constituent_density)
title("Constituent Density vs Altitude")
xlabel("Altitude (m)")
ylabel("Constituent Density (normalized)")



p_peak = 1; %laser peak power in W
t_rise = ceil((200e-6)/dt)*dt; %time from 10%->90% peak laser output

perc_10 = -2.1460; %z-score for 10% of max of gaussian. CONSTANT DO NOT CHANGE
perc_90 = -0.4590; %z-score for 90% of max of gaussian. CONSTANT DO NOT CHANGE
stdev = ceil((t_rise/(perc_90-perc_10))/dt)*dt; %corresponding stdev of gaussian given rise time

%generate power curve now that we know effective stdev of gaussian given
%the rise time between 10% and 90% peak output
power_curve = generate_gaussian(2*t_rise,stdev,[0,4*t_rise,dt],1,0);
power_curve = power_curve/(max(power_curve))*p_peak;

t_stop = 4*t_rise; %time when laser power is 0

% %total time window where events may be happening (accounting for ToF for last photons emitted by laser)
t_bounds = [dt,t_stop+2*z_bounds(2)/c+dt];
time = t_bounds(1):dt:t_bounds(2); %time axis of power curve

%append zeros on laser power curve so it's the same length as altitude, etc
power_curve = [power_curve, zeros([1,size(time,2)-size(power_curve,2)])];

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
A = 1; %receiver area in m^2
linear_QE = 0.4; %QE in low-signal intensity region
t_d = 100e-9; %dead time in seconds
paralyzable = 1; %boolean for paralyzable detector




%% generate return signal

%return signal received by PMT (before QE, saturation etc)
N_received_curve = generate_return(N_sent_curve, time, dt, altitudes, constituent_density, T_a, T_c, sigma_eff,A);

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
plot(ax1, time, N_received_curve, 'k-')
plot(ax1, time, N_total_received, 'g-')
plot(ax1, time, npar_detected_binned, 'r-')
plot(ax1, time, par_detected_binned, 'b-')
legend(["$N_{sig,received}$", "$N_{tot,received}$", "$N_{recorded,nonparalyzable}$","$N_{recorded,paralyzable}$"], 'interpreter', 'latex')
title("Number of Received Photons vs Time")
xlabel("Time (sec)")
ylabel("Num. Received Photons")



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




