function [temp_result] = temp_function(lowest_constituent,n_td,binned_detections)




pre_signal_bins = binned_detections(1:lowest_constituent-1);
tmp = discretize(pre_signal_bins, 0:max(pre_signal_bins)+1);

%bin count data. H(1) has the total number of pre-signal bins where there
%were 0 counted photons, and so on
H = zeros(1,max(pre_signal_bins));
for i=1:max(pre_signal_bins)+1
    H(i) = sum(tmp(:)==i);
end

%normalize H into a probability approximation of poisson
H = H/sum(H);
P_detect = sum(H(2:end));

x = 0:0.001:1;

temp_result = (x.^n_td - x.^(n_td+1));

% figure
% hold on
% plot([0,1],[P_detect,P_detect], 'k:')
% plot(x,(x.^n_td - x.^(n_td+1)), 'b-')
% xlabel("$P_n(\lambda=0)$", 'Interpreter', 'latex')
% ylabel("$P_{n,detect}$", "Interpreter", "latex")

%function we want to find a zero of. Effectively, we want to estimate the
%the noise photon receival probability that gives us the observed
%noise photon detection probability
fun_n = @(P_n) P_n^(n_td+1) - P_n^(n_td) + P_detect;

% tst = arrayfun(@(x) fun_n(x), 0:.001:1);
%
% figure
% plot(0:.001:1,tst)

% P_receive_no_n = fzero(fun_n, [0.0000000001,.99999999999]);

end
