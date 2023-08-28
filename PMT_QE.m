function [N_recorded] = PMT_QE(N_received, dt, t_d, linear_QE, paralyzable)
%Attenuates the number of received photons according to the specifications
%of the PMT/Discriminator system
%   Currently uses a very basic model that assumes photons are
%   equally-spaced within a dt window

f_received = N_received/dt;
f_dead = 1/t_d;

if paralyzable

    %within linear region
    if f_received <= f_dead

        N_recorded = N_received * linear_QE;

     %saturated paralyzable PMT - will output some suppresed number
    else
        
        %for now, assume that full extinction/paralyzation occurs at 2x
        %(1/deadtime), linear interpolation
        N_recorded = max([0, (f_dead*dt - (f_received-f_dead)*dt) * linear_QE]);
        
    end


else

    %within linear region
    if f_received <= f_dead

        N_recorded = N_received * linear_QE;

    %saturated nonparalyzable PMT - will output one pulse per deadtime
    else

        N_recorded = (f_dead * dt) * linear_QE;

    end

end



end