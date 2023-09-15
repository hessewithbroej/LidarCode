
stdev = 1;
cutoff_num_stds = 5;
x_start = 75;
x_stop = 120;
dx = 0.1;

x = x_start:dx:x_stop;

gauss = 1/(2*pi*stdev) * exp(-0.5*( ((x_start:dx:x_stop)-(stdev*cutoff_num_stds))/stdev ).^2);
gauss_shift = 1/(2*pi*stdev) * ((x_start:dx:x_stop)).*exp(-0.5*( ((x_start:dx:x_stop)-(stdev*cutoff_num_stds))/stdev ).^2);

alpha = 100;
skew_param = normcdf(alpha*(x/stdev-(stdev*cutoff_num_stds)/stdev));
gauss_skew = gauss .* skew_param;

gauss = gauss/sum(gauss);
gauss_skew = gauss_skew/sum(gauss_skew);



figure
plot(gauss)
hold on
plot(gauss_skew)


%% 
start_x = 1e+3;
stop_x = 120e+3;
dx = 1000;

curve = generate_gaussian(85000,3,[start_x,stop_x,dx],1,10);



x = start_x:dx:stop_x;

figure
plot(x,curve)
