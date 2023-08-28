%Lidar return simulation via lidar eqn

close all
clear all
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
%TODO: go from delta -> square -> curve
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
stdev = 3*dt;
cutoff_num_stds = 5; %truncate guassian after n stds

t_start = 0; %laser turn on time (s)
t_stop = t_start+(stdev/dt)*cutoff_num_stds*2*dt; %laser turn off time (s)
p_on = 1; %laser power in W when on



%total time window where events may be happening (accounting for ToF for last photons emitted by laser)
t_bounds = [0,t_stop+2*z_bounds(2)/c]; 
time = t_bounds(1):dt:t_bounds(2); %time axis of power curve
p_peak = 1; %peak power in W

%use this for a guassian-like power curve
power_curve = 1/(2*pi*stdev) * exp(-0.5*( ((t_start:dt:t_stop)-(stdev*cutoff_num_stds))/stdev ).^2);

%re-scale gaussian s.t. max(power_curve) = p_peak
power_curve = [power_curve/max(power_curve)*p_peak, zeros(1,size(time,2)-size(power_curve,2))];


% figure
% plot(time,power_curve);

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
paralyzable = 0; %boolean for paralyzable detector




%% main loop
N_received_curve = zeros(size(N_sent_curve));
N_recorded_curve_nonpar = N_received_curve;

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
        
        %
    end

    %PMT saturation effects
    N_recorded_curve_nonpar(i) = PMT_QE(N_received_curve(i), dt, t_d, linear_QE, 0);
    N_recorded_curve_par(i) = PMT_QE(N_received_curve(i), dt, t_d, linear_QE, 1);
end

f = figure;
ax1 = axes(f);
hold on
plot(ax1, time, N_received_curve, 'b-')
plot(ax1, time, N_recorded_curve_nonpar, 'r-')
plot(ax1, time, N_recorded_curve_par, 'g-')
% ax1.XColor = 'b';
% ax1.YColor = 'b';
% ax1.XAxisLocation = 'top';
% ax1.YAxisLocation = 'right';
title("N Received Photons vs Time")
xlabel("Time (sec)")
ylabel("N Received Photons")
% ax2 = axes(f);
% hold on
% plot(ax2, altitudes,constituent_density, 'r-')
% ax2.XColor = 'r';
% ax2.YColor = 'r';
% ax2.Color = 'none';
% ax2.XAxisLocation = 'bottom';
% ax2.YAxisLocation = 'left';
% ax1.Box = 'off';
% ax2.Box = 'off';


% plot(time, N_sent_curve)
% ylabel("N Sent Photons")

figure
xcorr_res = xcorr(constituent_density, power_curve);
plot(xcorr_res)
xlabel("Time (indexed, shifted)")
ylabel("Cons. Density $\star$ Pow. Curve", 'Interpreter','latex')
title("Cross correlation of power curve and density profile")
N_sent
N_received = sum(N_received_curve)
ratio = N_received/N_sent

