%Lidar return simulation via lidar eqn
%% freshen workspace
close all
clearvars
clc

%% populate lists of constants
%-------------------------------------------------
%fundamental constants
lambda = 374e-9; %laser wl in m
c = 299792498; %speed of light
h = 6.626e-34; %planck constant j*s
E_photon = h*c/lambda; %photon energy

constants = [lambda, c, h];
%-------------------------------------------------
%atmospheric 'constants'
T_a = 0.8;
T_c = 0.95;
%cross section parameters
%TODO: Add 1st-order temp and wind dependence --> will require vertical
%profiles of temp and vertical wind
sigma_eff = 1;
f_n = 1e+5;

atmospheric_parameters = [T_a, T_c, sigma_eff, f_n];
%-------------------------------------------------
%receiver parameters
A = 0.75; %receiver area in m^2
linear_QE = 0.4; %QE in low-signal intensity region
t_d = 100e-9; %dead time in seconds

receiver_parameters = [A,linear_QE,t_d];
%-------------------------------------------------
%laser parameters
t_rise = 200e-6; %10%-90% rise time
p_peak = 1; %peak power in watts

laser_parameters = [t_rise,p_peak];

%% generate constituent density function of altitude
z_bounds = [1e+3,120e+3];
dz = 1000;
altitudes = z_bounds(1):dz:z_bounds(2); %altitude in m
constituent_density = generate_gaussian(85000,3000, [z_bounds(1),z_bounds(2),dz], 1, 5);

figure
plot(altitudes,constituent_density)
title("Constituent Density vs Altitude")
xlabel("Altitude (m)")
ylabel("Constituent Density (normalized)")


%% generate laser power function

dt = 2*dz/c;
t_rise = ceil((t_rise)/dt)*dt; %time from 10%->90% peak laser output, in interger multiple of dt

perc_10 = -2.1460; %z-score for 10% of max of gaussian. CONSTANT DO NOT CHANGE
perc_90 = -0.4590; %z-score for 90% of max of gaussian. CONSTANT DO NOT CHANGE
stdev = ceil((t_rise/(perc_90-perc_10))/dt)*dt; %corresponding stdev of gaussian given rise time

%generate power curve now that we know effective stdev of gaussian given
%the rise time between 10% and 90% peak output
power_curve = generate_gaussian(2*t_rise,stdev,[0,4*t_rise,dt],1,0);

%scale power curve to p_peak
power_curve = power_curve/(max(power_curve))*p_peak;

%% generate appropriate time domain
t_stop = 4*t_rise; %time when laser power falls to 0

% %total time window where events may be happening (accounting for ToF for last photons emitted by laser)
t_bounds = [dt,t_stop+2*z_bounds(2)/c+dt];
time = t_bounds(1):dt:t_bounds(2); %time axis of power curve

%append zeros on laser power curve so it's the same length as altitude, etc
power_curve = [power_curve, zeros([1,size(time,2)-size(power_curve,2)])];

%% convert laser power to N_sent photons
%generate curve of N_sent photons over time
energy_curve = power_curve*dt; %total energy emitted by laser (j) at each timestep
N_sent_curve = floor(energy_curve/E_photon); %total number of photons emitted at each timestep
N_sent = sum(N_sent_curve); %total number of photons emitted during pulse

figure
plot(time,N_sent_curve);
title("N Sent Photons vs Time")
xlabel("Time (sec)")
ylabel("N Sent Photons")


%% generate return signals - received by PMT and detected by PMT
paralyzable_detections = {};
nonparalyzable_detections = {};
received_noisy = {};
received_signal = {};

arrival_times_cell = {};
par_detected_logicals = {};
npar_detected_logicals = {};

N_runs = 1000;
runtimes = zeros(1,N_runs);



for i=1:N_runs
    tic
    [received_curves, detected_curves, arrival_times, detected_timestamps] = simulate_lidar(constants, altitudes, time, constituent_density, N_sent_curve, atmospheric_parameters, receiver_parameters);
    paralyzable_detections{end+1} = detected_curves{1};
    nonparalyzable_detections{end+1} = detected_curves{2};
    received_noisy{end+1} = received_curves{1};
    received_signal{end+1} = received_curves{2};

    arrival_times_cell{end+1} = arrival_times;
    par_detected_logicals{end+1} = detected_timestamps{1};
    npar_detected_logicals{end+1} = detected_timestamps{2};

    runtimes(i) = toc;
    disp("Completed run: " + string(i) + " | Execution time: " + string(runtimes(i)) + " | Net Runtime: " + string(sum(runtimes)) + " | Estimated time remaining: " + string( (N_runs-i)*mean(runtimes(runtimes ~= 0)) ) )
end
disp("Completed all runs in: " + string(sum(runtimes)))
save(string(datetime(), "yyyyMMdd_HHmmss"), "constants", "atmospheric_parameters", "receiver_parameters", "laser_parameters", "paralyzable_detections", "nonparalyzable_detections", "received_noisy", "received_signal", "arrival_times_cell", "par_detected_logicals", "npar_detected_logicals")
%% process & plot


% figure
% hold on
% plot(time,received_curves{1},'g-');
% plot(time,received_curves{2},'k-');
% plot(time,detected_curves{2},'r-');
% plot(time,detected_curves{1},'b-');
% legend(["$N_{sig,received}$", "$N_{tot,received}$", "$N_{recorded,nonparalyzable}$","$N_{recorded,paralyzable}$"], 'interpreter', 'latex')
% title("Number of Received Photons vs Time")
% xlabel("Time (sec)")
% ylabel("Num. Received Photons")
%
% figure
% hold on
% plot(time,received_curves{1}/dt,'g-');
% plot(time,received_curves{2}/dt,'k-');
% plot(time,detected_curves{2}/dt,'r-');
% plot(time,detected_curves{1}/dt,'b-');
% plot([min(time),max(time)], [1/t_d,1/t_d], 'k--')
% legend(["$f_{sig,received}$", "$f_{tot,received}$", "$f_{recorded,nonparalyzable}$","$f_{recorded,paralyzable}$", "$f_{dead}$"], 'interpreter', 'latex')
% title("Frequency of Received Photons vs Time")
% xlabel("Time (sec)")
% ylabel("Freq. Received Photons (Hz)")


%% liu correction
tic


n_td = t_d/dt;
correction_rebinning_factor = 1;

lowest_constituent = find(constituent_density~=0, 1, 'first'); % time bins prior to this index will have only noise photons (based on ToF)
t_signal = time(lowest_constituent);
sub_times = time(1):correction_rebinning_factor*t_d:time(end);
[~,lowest_constituent_subbin] = min(abs(t_signal-sub_times));

n_td = 1/correction_rebinning_factor; %we intentionally upscale the resolution to make dt_subbins smaller
n_n = f_n * correction_rebinning_factor*t_d;

net_detected = zeros(size(sub_times));
H = zeros(size(sub_times));
for i=1:size(arrival_times_cell,2)

    par_detected_timestamps = arrival_times_cell{i}(logical(par_detected_logicals{i}));
    par_detected_subbinned = bin_photons(par_detected_timestamps, sub_times);
    net_detected = net_detected + par_detected_subbinned;
end
%avg number of photons detected in each bin
H = net_detected/size(arrival_times_cell,2); %this plus the preceding for loop is Liu's equation 11
% H = 1-exp(-net_detected);
H = movmean(H,10);
rt = toc
figure
hold on
plot(sub_times,H)
%%
[P_receive_no_n,P_detect_no_n] = calculate_noise_probability(lowest_constituent_subbin,n_td,H); %Liu's equation 12 inside here

par_detected_corrected_subbin = zeros(size(H));
par_detected_corrected_subbin_saturated = zeros(size(H));
par_detected_corrected_subbin_unsaturated = zeros(size(H));

P_receive_no_ns_list = [repmat(P_receive_no_n, [1,lowest_constituent_subbin]), zeros([1,size(par_detected_subbinned,2)-lowest_constituent_subbin])];
P_receive_no_ns_list_saturated = [repmat(P_receive_no_n, [1,lowest_constituent_subbin]), zeros([1,size(par_detected_subbinned,2)-lowest_constituent_subbin])];
P_receive_no_ns_list_unsaturated = [repmat(P_receive_no_n, [1,lowest_constituent_subbin]), zeros([1,size(par_detected_subbinned,2)-lowest_constituent_subbin])];


P_receive_no_ns = 1-H(lowest_constituent_subbin)/(P_receive_no_n^n_td); %liu eqn 13
P_receive_no_ns_list(lowest_constituent_subbin) = P_receive_no_ns;
P_receive_no_ns_list_saturated(lowest_constituent_subbin) = P_receive_no_ns;
P_receive_no_ns_list_unsaturated(lowest_constituent_subbin) = P_receive_no_ns;

par_detected_corrected_subbin(lowest_constituent_subbin) = max([0,-log(P_receive_no_ns)-n_n]); %liu eqn 14
par_detected_corrected_saturated_subbin(lowest_constituent_subbin) = max([0,-log(P_receive_no_ns)-n_n]);
par_detected_corrected_unsaturated_subbin(lowest_constituent_subbin) = max([0,-log(P_receive_no_ns)-n_n]);


for i=lowest_constituent_subbin+1:max(size(par_detected_subbinned))

    %applying liu eqn 14 to all bins i > m
    [P_receive_no_ns_saturated, P_receive_no_ns_unsaturated] = calculate_receival_probability(H(i), n_td); 
    P_receive_no_ns = 1-H(i)/P_receive_no_ns_list(i-1);

    P_receive_no_ns_list(i) = P_receive_no_ns;
    P_receive_no_ns_list_saturated(i) = P_receive_no_ns_saturated;
    P_receive_no_ns_list_unsaturated(i) = P_receive_no_ns_unsaturated;

    par_detected_corrected_subbin(i) = (max([0,-log(P_receive_no_ns)-n_n]));
    par_detected_corrected_subbin_saturated(i) = (max([0,-log(P_receive_no_ns_saturated)-n_n]));
    par_detected_corrected_subbin_unsaturated(i) = (max([0,-log(P_receive_no_ns_unsaturated)-n_n]));

    if ~isreal(par_detected_corrected_subbin(i)) || ~isreal(par_detected_corrected_subbin_saturated(i)) || ~isreal(par_detected_corrected_subbin_unsaturated(i))
        pause(0.01)
    end

end

par_detected_corrected = zeros(size(time));
par_detected_corrected_saturated = zeros(size(time));
par_detected_corrected_unsaturated = zeros(size(time));

% now we need to re-bin our corrected data into original bin sizes
t = sub_times(1);
ind = 1;
while t < sub_times(end)

    high_ind = find(sub_times < t+dt, 1, 'last');
    low_ind = find(sub_times >= t, 1, 'first');

    par_detected_corrected(ind) = sum(par_detected_corrected_subbin(low_ind:min([high_ind,size(par_detected_corrected_subbin,2)])));
    par_detected_corrected_saturated(ind) = sum(par_detected_corrected_subbin_saturated(low_ind:min([high_ind,size(par_detected_corrected_subbin,2)])));
    par_detected_corrected_unsaturated(ind) = sum(par_detected_corrected_subbin_unsaturated(low_ind:min([high_ind,size(par_detected_corrected_subbin,2)])));


    ind = ind+1;
    t = t+dt;
end





%%
% figure
% hold on
% plot(sub_times, movmean(H,20))
% plot(sub_times, movmean(1-P_receive_no_ns_list,20))
%
% figure
% plot(sub_times, movmean((1-P_receive_no_ns_list)-H,20))
%
% figure
% hold on
% plot(sub_times, movmean(H,1))
% plot(sub_times, movmean(H,10))
% plot(sub_times, movmean(H,20))
% plot(sub_times, movmean(H,50))

%%
%
%
mean_paralyzable_detections = mean(cell2mat(paralyzable_detections'));
mean_received_signal = mean(cell2mat(received_signal'));
mean_received_noisy = mean(cell2mat(received_noisy'));

%correction error metric
norm_area= sum( abs(par_detected_corrected/linear_QE - mean_received_signal)*dt )/sum(mean_received_signal*dt);

figure
hold on
plot(time, round(mean_received_signal), 'k-')
plot(time, round(mean_received_noisy), 'g-')
plot(time, round(mean_paralyzable_detections), 'r-')
% plot(time, round(par_detected_corrected), 'b-')
plot(time, round(par_detected_corrected/linear_QE),'m-')
plot(time, round(mean_paralyzable_detections/linear_QE), 'y-')
plot(time, round(par_detected_corrected_saturated/linear_QE),'c-')
plot(time, round(par_detected_corrected_unsaturated/linear_QE),'LineStyle','-','Color',[.3,0,.15])
xlabel("Time (sec)")
ylabel("Photon Counts")
legend(["$N_{r,sig}$", "$N_{r,sig+noise}$", "$N_{d}$", "$N_{d,corr (Liu)}$", "$N_{d,corr,  (naive)}$", "$N_{d,corr (sat)}$", "$N_{d,corr (unsat)}$"], "interpreter", "latex")

figure
hold on
plot(time, round(mean_received_signal)/dt, 'k-')
plot(time, round(mean_received_noisy)/dt, 'g-')
plot(time, round(mean_paralyzable_detections)/dt, 'r-')
% plot(time, round(par_detected_corrected), 'b-')
plot(time, round(par_detected_corrected/linear_QE)/dt,'m-')
plot(time, round(mean_paralyzable_detections/linear_QE)/dt, 'b-')
% plot(time, round(par_detected_corrected_saturated/linear_QE)/dt,'c-')
% plot(time, round(par_detected_corrected_unsaturated/linear_QE)/dt,'LineStyle','-','Color',[.3,0,.15])
plot([time(1),time(end)],[1/t_d,1/t_d], "k:")
xlabel("Time (sec)")
ylabel("Count Frequency (Hz)")
legend(["$f_{r,sig}$", "$f_{r,sig+noise}$", "$f_{d}$", "$f_{d,corr (Liu)}$", "$f_{d,corr,  (naive)}$", "$f_d$"], "interpreter", "latex")
% legend(["$f_{r,sig}$", "$f_{r,sig+foise}$", "$f_{d}$", "$f_{d,corr (Liu)}$", "$f_{d,corr,  (naive)}$", "$f_{d,corr (sat)}$", "$f_{d,corr,  (unsat)}$", "$f_d$"], "interpreter", "latex")
title("Simulated PMT Saturation Correction (A=" + string(receiver_parameters(1)) + ")")
ylim([0,1.2*(1/t_d)])
text(0,0.6*(1/t_d),"Error: " + string(norm_area))
