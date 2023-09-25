%Lidar return simulation via lidar eqn
%% freshen workspace
close all
clearvars
clc

%% populate lists of constants
%-------------------------------------------------
%fundamental constants
lambda = 532-9; %laser wl in m
c = 299792498; %speed of light
h = 6.626e-34; %planck constant j*s
E_photon = h*c/lambda; %photon energy

constants = [lambda, c, h];
%-------------------------------------------------

linear_QE = 0.4; %QE in low-signal intensity region
t_d = 6e-9; %dead time in seconds

receiver_parameters = [linear_QE,t_d];
%-------------------------------------------------
%laser parameters
dt = 200e-13;
FWHM = 1e-9;
stdev = FWHM/2.355; %corresponding stdev of gaussian given rise time

E_pulse = 1.8e-6;

n_stdevs = 10;
%generate power curve now that we know effective stdev of gaussian given
%the rise time between 10% and 90% peak output
power_curve = generate_gaussian(0.5*n_stdevs*stdev,stdev,[dt,n_stdevs*stdev,dt],1,0);

%% generate appropriate time domain
t_stop = n_stdevs*stdev; %time when laser power falls to 0

t_travel = 30; %# of time bins for travel to receiver
% %total time window where events may be happening (accounting for ToF for last photons emitted by laser)
t_bounds = [dt,t_stop+t_travel*dt];
time = t_bounds(1):dt:t_bounds(2); %time axis of power curve

f_n = 4e+6;

power_curve = [power_curve, zeros(1,t_travel)];

%% convert laser power to N_sent photons
%generate curve of N_sent photons over time
energy_curve = E_pulse*power_curve/sum(power_curve);

atten = 10^-19;

energy_curve = atten*energy_curve;
N_sent_curve = (energy_curve/E_photon); %total number of photons emitted at each timestep
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

N_runs = 20000;
runtimes = zeros(1,N_runs);



for i=1:N_runs
    tic
    [received_curves, detected_curves, arrival_times, detected_timestamps] = simple_return(N_sent_curve, dt, time, t_travel, t_d, linear_QE, f_n);
    paralyzable_detections{end+1} = detected_curves{1};
    nonparalyzable_detections{end+1} = detected_curves{2};
    received_noisy{end+1} = received_curves{2};
    received_signal{end+1} = received_curves{1};

    arrival_times_cell{end+1} = arrival_times;
    par_detected_logicals{end+1} = detected_timestamps{1};
    npar_detected_logicals{end+1} = detected_timestamps{2};

    runtimes(i) = toc;
    if mod(i,10) == 0
        disp("Completed run: " + string(i) + " | Execution time: " + string(runtimes(i)) + " | Net Runtime: " + string(sum(runtimes)) + " | Estimated time remaining: " + string( (N_runs-i)*mean(runtimes(runtimes ~= 0)) ) )
    end
end
disp("Completed all runs in: " + string(sum(runtimes)))
% save(string(datetime(), "yyyyMMdd_hhmmss"), "constants", "atmospheric_parameters", "receiver_parameters", "laser_parameters", "paralyzable_detections", "nonparalyzable_detections", "received_noisy", "received_signal", "arrival_times_cell", "par_detected_logicals", "npar_detected_logicals")
%% process (generate averages) & plot

mean_paralyzable_detections = mean(cell2mat(paralyzable_detections'));
mean_received_signal = mean(cell2mat(received_signal'));
mean_received_noisy = mean(cell2mat(received_noisy'));

f_final=figure
hold on
plot(time,mean_received_noisy, 'g-')
plot(time,mean_received_signal,'k-')
plot(time,mean_paralyzable_detections,'b-')

% 
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
% tic

n_td = t_d/dt;
correction_rebinning_factor = 1/n_td;

% lowest_constituent = find(constituent_density~=0, 1, 'first'); % time bins prior to this index will have only noise photons (based on ToF)
lowest_constituent = t_travel;
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
H = net_detected/size(arrival_times_cell,2);
% H = 1-exp(-net_detected);
% H = movmean(H,10);
% rt = toc
figure
hold on
plot(sub_times,H)
%%
[P_receive_no_n,P_detect_no_n] = calculate_noise_probability(lowest_constituent_subbin,n_td,H);

par_detected_corrected_subbin = zeros(size(H));

P_receive_no_ns_list = [repmat(P_receive_no_n, [1,lowest_constituent_subbin]), zeros([1,size(par_detected_subbinned,2)-lowest_constituent_subbin])];

P_receive_no_ns = 1-H(lowest_constituent_subbin)/(P_receive_no_n^n_td);
P_receive_no_ns_list(lowest_constituent_subbin) = P_receive_no_ns;

par_detected_corrected_subbin(lowest_constituent_subbin) = max([0,-log(P_receive_no_ns)-n_n]);


for i=lowest_constituent_subbin+1:max(size(par_detected_subbinned))


    P_receive_no_ns = 1-H(i)/P_receive_no_ns_list(i-1);
    P_receive_no_ns_list(i) = P_receive_no_ns;

    par_detected_corrected_subbin(i) = (max([0,-log(P_receive_no_ns)-n_n]));

    if ~isreal(par_detected_corrected_subbin(i))
        pause(0.01)
    end

end

par_detected_corrected = zeros(size(time));

% now we need to re-bin our corrected data into original bin sizes
t = sub_times(1);
ind = 1;
while t < sub_times(end)

    high_ind = find(sub_times < t+dt, 1, 'last');
    low_ind = find(sub_times >= t, 1, 'first');

    par_detected_corrected(ind) = sum(par_detected_corrected_subbin(low_ind:min([high_ind,size(par_detected_corrected_subbin,2)])));
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

figure
hold on
plot(time, mean_received_signal, 'g-')
plot(time, mean_received_noisy, 'k-')
plot(time, mean_paralyzable_detections, 'b-')
plot(time, (par_detected_corrected), 'r-')
plot(time, (par_detected_corrected/linear_QE),'m-')
% plot(time, (par_detected_binned/linear_QE), 'y-')

legend({"P_{s}", "P_{s+n}", "P_{d}", "P_{d,corrected}", "P_{d,corrected,2}"})