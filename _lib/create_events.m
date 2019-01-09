%% Create vector of event times with gaussian variation around mean interval
function input_vec = create_events(stim_start,total_stim_dur,EPSC_frq,vec_factor,refrac,Rave)
    interval = 1000/EPSC_frq; % ms
    stim_end = stim_start + total_stim_dur;
    x_thresh  = 1 - Rave/EPSC_frq; % prob of spiking
    
    input_vec = [];
    input_vec(1) = stim_start + invl(interval/2,vec_factor);
    if input_vec(1)<stim_start
        input_vec(1) = stim_start;
    end
    it = 1;
    stim_interval = stim_start;
    while (stim_interval<stim_end-interval)
        stim_interval = stim_interval + interval;
        if rand(1)>x_thresh
            it = it+1;
            input_vec(it) = stim_interval + invl(interval/2,vec_factor);
            while (input_vec(it)-input_vec(it-1))<=refrac
                stim_interval = stim_interval + interval; % add 1 interval
                input_vec(it) = stim_interval + invl(interval/2,vec_factor);
            end
        end        
    end
    function x = invl(mean,vec_factor)
        if (mean <= 0.)
            mean = .01; %(ms)
        end
        std  = mean/vec_factor;
        x = mean+std*randn(1);
        if x>=interval
            x=mod(x,interval);
        elseif x<0
            x=mod(x,interval)+interval;
        end	
    end
end
