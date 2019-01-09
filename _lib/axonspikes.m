function [nspikes, tspikes, vaxon] = axonspikes(vsoma,dt,vth,vreset,tref,tau_axon,gc_c_axon)
    % axonspikes - compute spikes produced in axon when axon has no influence on soma
    iiref   = tref/dt;

    % initialize
    ref     = 0; % flag refactory
    nspikes = 0;
    iispike = -iiref-1;
    vaxon = vreset*ones(size(vsoma));
    tspikes = [];

    for ii = 2:length(vsoma)
        if (ref == 1)
            vaxon(ii) = vreset;
            if (ii-iispike>iiref)
                ref = 0;
            end
        else
            dvaxon = gc_c_axon*(vsoma(ii)-vaxon(ii-1)) + (vreset - vaxon(ii-1))/tau_axon;
            vaxon(ii) = vaxon(ii-1) + dt*dvaxon;

            if(vaxon(ii)>vth)
                ref = 1;
                iispike = ii;
                nspikes = nspikes + 1;
                tspikes(nspikes) = ii*dt;
            end
        end
    end
end

