%average separation between events = 1/num_arrivals
num_arrivals = round(10.^(1:0.5:7));
t_d = 100e-6;
f_d = 1/t_d;
linear_QE = 0.4;


%repeat to reduce variability
sum_par = zeros(size(num_arrivals));
sum_npar = zeros(size(num_arrivals));
n_runs = 100;
for i=1:n_runs

    f_d = 1/t_d;
    f_sig = num_arrivals;

    arrival_times = arrayfun(@(x) sort(rand([1,x])), num_arrivals, 'UniformOutput', false); %n photons arrive randomly between the interval (0,1)

    par_result = cellfun(@(x) sum(PMT_QE(x, t_d, linear_QE, 1)), arrival_times);
    sum_par = sum_par + par_result;

    npar_result = cellfun(@(x) sum(PMT_QE(x, t_d, linear_QE, 0)), arrival_times);
    sum_npar = sum_npar + npar_result;
end

avg_par_result = sum_par/n_runs;
avg_npar_result = sum_npar/n_runs;

figure
hold on
plot(num_arrivals,avg_par_result,'b-')
plot(num_arrivals,avg_npar_result,'r-')
plot([f_d, f_d], [0,max(avg_npar_result)], 'k--')
ax = gca;
ax.XScale = 'log';
xlabel("Input Frequency")
ylabel("Output Frequency")


