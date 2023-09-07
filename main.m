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
constituent_density = [repmat(0,[1,75]), linspace(0,100,10), repmat(100,[1,10]), linspace(100,1,25)]; %unitless density of constituent species
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




%% main loop
N_received_curve = zeros(size(N_sent_curve));
N_recorded_curve_npar = N_received_curve;
N_recorded_curve_par = N_received_curve;

%atmospheric absorptions impact signal at all altitudes equally (when
%ignoring sigma_eff variation with alt/temp)
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
    %
    %     %want to convert n_received_curve into a list of photon arrival
    %     %timestamps.
    %
    %     upres_factor = 10;
    %     t_lower = (i-1)*dt;
    %     t_upper = t;
    %     t_interp = t_lower:dt*(1/upres_factor):t_upper;
    %     if i==1
    %         prev_n_received = 0;
    %     else
    %         prev_n_received = N_received_curve(i-1);
    %     end
    %
    %     N_exp = zeros(1,size(t_interp,2)-1);
    %     arrival_times = [];
    %     prev_reached_bin = 1;
    %     running_sum = 0;
    %     %linear interpolation at upscaled resolution.
    %     for j=1:size(t_interp,2)-1
    %         % Number of photons expected within higher resolution bin
    %         N_exp(j) = (prev_n_received+(j*(1/upres_factor))*(N_received_curve(i)-prev_n_received))/upres_factor;
    %         %we want to track the cummulative number of photons we expect to
    %         %see through our current bin.
    %         running_sum = running_sum + N_exp(j);
    %
    %         %we should have seen a photon by this point (expected N photons > 1)
    %         % Choose one of the bins in the range from our last received photon to put this photon in.
    %         %Choice should be weighted by the expected number of photons
    %         %received in that bin.
    %         %If we happen to get multiple photons expected in a single step,
    %         %distribute them all in the same manner.
    %         while running_sum >= 1
    %             %select the bin in which we'll log photon, using the expected
    %             %photon count in each bin as weight
    %             %Randsample acts differently if the first argument is a single
    %             %number instead of a vector, need both cases
    %             if prev_reached_bin ~= j
    %                 bin = randsample(prev_reached_bin:j,1,true, N_exp(prev_reached_bin:j));
    %             else
    %                 bin = j;
    %             end
    %             arrival_times(end+1) = t_lower + (t_upper-t_lower).*rand;
    %             running_sum = running_sum-1;
    %             prev_reached_bin = j;
    %         end
    %     end



end

upscale_factor = 10;
arrival_times = sort(disperse_photons(N_received_curve,time,upscale_factor));

%PMT saturation effects. Boolean of whether a photon was detected
npar_detected = PMT_QE(sort(arrival_times), t_d, linear_QE, 0);
par_detected = PMT_QE(sort(arrival_times), t_d, linear_QE, 1);

%pull only the detected timestamps:
npar_detected_timestamps = arrival_times(logical(npar_detected));
par_detected_timestamps = arrival_times(logical(par_detected));

%convert timestamps back to counts per bin
npar_detected_binned = bin_photons(npar_detected_timestamps,time);
par_detected_binned = bin_photons(par_detected_timestamps,time);

% N_recorded_curve_npar(i) = sum(npar_detected);
% N_recorded_curve_par(i) = sum(par_detected);

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




