function [P_receive_no_ns_saturated, P_receive_no_ns_unsaturated] = calculate_receival_probability(P_detect_ns, n_td)


plot_flag = false;
x = 0:0.001:1;
if plot_flag
    figure
    hold on
    plot(x,x.^(n_td)-x.^(n_td+1), 'b-')
    plot([0,1],[P_detect_ns,P_detect_ns])
    xlabel("$P_r(\lambda=0)$", 'Interpreter','latex')
    ylabel("$P_d$",'Interpreter','latex')
    figure
    plot(x,x.^(n_td)-x.^(n_td+1)-(P_detect_ns));
end


%function we want to find a zero of. Effectively, we want to estimate the
%the noise photon receival probability that gives us the observed
%noise photon detection probability. Note P_n is hte probability a bin does
%NOT receive a photon
fun_n = @(P_ns) P_ns^(n_td) - P_ns^(n_td+1) - (P_detect_ns);

%When detection probability in the pre-signal region is very low, 
%it means that either 
% 1) noise is extremely low. Very few photons received -> very few photons
% detected. (P_receive ~ 0 -> P_detect ~ 0)
% OR
% 2) noise is extremely high, exceeding dead time & suppressing the
% detection probability. (P_receive ~ 1 -> P_detect ~ 0)

if P_detect_ns > 0

P_receive_no_ns_saturated = fzero(fun_n, [0.0001,0.5]);
P_receive_no_ns_unsaturated = fzero(fun_n, [0.5,0.99999]);

else
P_receive_no_ns_saturated = 0;
P_receive_no_ns_unsaturated = 1;
end

end
