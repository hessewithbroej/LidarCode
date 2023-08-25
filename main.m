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
% constituent_density = [repmat(0,[1,75]), linspace(1,100,10), repmat(100,[1,10]), linspace(100,1,25)]; %unitless density of constituent species
constituent_density = [zeros(1,75), ones(1,45)];
constituent_density = constituent_density/sum(constituent_density); %density/concentration normalized such that full series sums to 1

dt = 2*dz/c; %time step for data collection. Gives altitude resolution of dz (accounting for there and back travel)


figure
plot(altitudes,constituent_density)

%TODO: go from delta -> square -> curve
t_start = 0; %laser turn on time (s)
t_stop = t_start+dt; %laser turn off time (s)
p_on = 0.1; %laser power in W when on



t_bounds = [0,t_stop+2*z_bounds(2)/c]; %total time window where events may be happening (accounting for ToF for last photons emitted by laser)


time = t_bounds(1):dt:t_bounds(2); %time axis of power curve
power_curve = [p_on * ones(1,nearestDiv(t_stop-t_start,dt)), zeros(1,nearestDiv(t_bounds(2)+dt-t_stop,dt))];

figure
plot(time,power_curve);

energy_curve = power_curve*dt; %total energy emitted by laser (j) at each timestep
N_sent_curve = floor(energy_curve/E_photon); %total number of photons emitted at each timestep
N_sent = sum(N_sent_curve); %total number of photons emitted during pulse

figure
plot(time, energy_curve);
figure
plot(time,N_sent_curve);

%attenuation parameters
T_a = 0.8;
T_c = 0.95;
transition_altitude = 76;

%cross section parameters
%TODO: Add 1st-order temp and wind dependence
sigma_eff = 1;

%receiver parameters
A = 1; %receiver area in m^2
QE = 0.4; %TODO: make signal-dependent

%% main loop
N_received_curve = zeros(size(N_sent_curve));


%atmospheric absorptions impact signal at all altitudes equally (when
%ignoring sigma_eff variation with alt/temp)
atmospheric_absorption_factor = T_a*T_c*sigma_eff*T_c*T_a; %scaling factor to account for atmospheric absorptions



%compute received photon at each time bin
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

        
        N_received_curve(i) = N_received_curve(i) + ...
            N_sent_curve(i-j+1) * atmospheric_absorption_factor * (A/(4*pi*z^2)) * constituent_density(ind);
        
        
        

    end

end

figure
plot(time, N_received_curve)

N_sent
N_received = sum(N_received_curve)
ratio = N_received/N_sent

