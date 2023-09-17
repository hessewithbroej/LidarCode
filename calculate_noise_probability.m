function [P_receive_no_n, P_detect_no_n] = calculate_noise_probability(lowest_constituent,n_td,binned_detections)


pre_signal_bins = binned_detections(1:lowest_constituent-1);
% tmp = discretize(pre_signal_bins, 0:max(pre_signal_bins)+1);
% %bin count data. H(1) has the total number of pre-signal bins where there
% %were 0 counted photons, and so on
% H = zeros(1,max(pre_signal_bins));
% for i=1:max(pre_signal_bins)+1
%     H(i) = sum(tmp(:)==i);
% end

%pre-signal, all bins are independent so we can estimate the probability
%any given bin detects a photon as the average the detection probability of
%all pre-signal bins
P_detect_no_n = 1-mean(pre_signal_bins);
% avg_n_received = sum(pre_signal_bins)/N_sims;
% P_detect_no_n = exp(-avg_n_received);


plot_flag = true;
x = 0:0.001:1;
if plot_flag
    figure
    hold on
    plot(x,x.^(n_td)-x.^(n_td+1), 'b-')
    plot([0,1],[1-P_detect_no_n,1-P_detect_no_n])
    xlabel("$P_r(\lambda=0)$", 'Interpreter','latex')
    ylabel("$P_d$",'Interpreter','latex')
    figure
    plot(x,x.^(n_td)-x.^(n_td+1)-(1-P_detect_no_n));
end


%function we want to find a zero of. Effectively, we want to estimate the
%the noise photon receival probability that gives us the observed
%noise photon detection probability. Note P_n is hte probability a bin does
%NOT receive a photon
fun_n = @(P_n) P_n^(n_td) - P_n^(n_td+1) - (1-P_detect_no_n);

%When detection probability in the pre-signal region is very low, 
%it means that either 
% 1) noise is extremely low. Very few photons received -> very few photons
% detected. (P_receive ~ 0 -> P_detect ~ 0)
% OR
% 2) noise is extremely high, exceeding dead time & suppressing the
% detection probability. (P_receive ~ 1 -> P_detect ~ 0)
% For our case, noise is very low, so limit our search for a zero to [0.3,1]
% fzero gets mad if both endpoints are the same sign (thinks that no zero
% exists)
P_receive_no_n = fzero(fun_n, [0.3,0.999]);

end
