function [curve] = generate_gaussian(mu,stdev,xs,normalize,alpha)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


x_start = xs(1);
x_stop = xs(2);
dx = xs(3);

%generate gaussian
curve = 1/sqrt(2*pi*stdev) * exp(-0.5*( ((x_start:dx:x_stop)-mu)/stdev ).^2);

%skew gaussian by parameter alpha
if alpha ~= 0
    curve = curve .* normcdf(alpha*(x_start:dx:x_stop),  alpha*mu, stdev);
end

if normalize
    curve = curve/sum(curve);
end



end