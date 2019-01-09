% Michiel Remme
% October 2018
% demo run of MSO model plotting inputs, currents, voltage as in Figure 2A from the article:
% Function and energy consumption constrain neuronal biophysics as in a canonical computation: coincidence detection
% Michiel W. H. Remme
% John Rinzel
% Susanne Schreiber
% PLoS Computational Biology 2018

clearvars

% compile c-code for matlab
mex iterate_MSO.c
addpath _lib

% load parameters
run set_parameters

%% SIMULATION BATCH

ITD_range   = [0, 0.5]; % (ms)
EPSG_range  = [0.0196, 0.5]; % (uS)

for k_ITD = 1:length(ITD_range)
    ITD      = ITD_range(k_ITD);    
    run create_input_vectors
    gvec_ITD = zeros(size(gvec));

    k_EPSG = 1;
    EPSG_iters = length(EPSG_range);
    
    for k_EPSG = 1:EPSG_iters       
        A_EPSG = EPSG_range(k_EPSG);
        fprintf('\n')
        fprintf('ITD = %g ms\n',ITD);
        fprintf('EPSG amplitude = %g uS\n',A_EPSG);
        gvec_ITD(1:Ninputs/2,:) = A_EPSG*gvec(1:Ninputs/2,:);
        gvec_ITD(Ninputs/2+1:Ninputs,1+ceil(ITD/dt):end) = A_EPSG*gvec(Ninputs/2+1:Ninputs,1:end-ceil(ITD/dt));  % INTRODUCE ITD BETWEEN IPSI AND CONTRA INPUTS
                
        % RUN SIMULATION
        cellparams  = [L rho tm*1e3 gamma_m tw_fac v_init Ena Ek Es RNinf*1e-6 sdAreaRatio];                
        simparams   = [tend dt dx];
        [vsoma,vdend,ina_syn_result,ik_syn_result,ina_leak_result,ik_leak_result,ik_klt_result] = iterate_MSO(simparams, cellparams, locvec, gvec_ITD);
        areaCorrection  = (area_s+2*len*diam*pi)/(2*nseg + sdAreaRatio)*1e6;
        ina_leak_result = ina_leak_result*gnapas*areaCorrection;
        ik_leak_result  = ik_leak_result*gkpas*areaCorrection;
        ik_klt_result   = ik_klt_result*gpas*areaCorrection;                    
        [nspikes, tspikes, vaxon] = axonspikes(vsoma,dt,vth,v_init,tref,tau_axon,gc_c_axon); % computes number of output spikes and converts to rates

        fprintf('%d spikes in %d ms\n',nspikes,tend);
        
        % PLOT RESULTS
        figure(1)
        clf

        t1 = tend-100;
        t2 = tend;

        gtmp = gvec_ITD(12:-1:1,:)'./A_EPSG;
        input_conductances = [gtmp + repmat(1.5*[0:11],tend/dt,1)]';

        subplot(311)
        plot(t,input_conductances(1:6,:),'g',t,input_conductances(7:12,:),'b')
        axis([t1 t2 -1 19])
        xlabel('Time (ms)')
        ylabel('Synaptic conductances')

        subplot(312)
        plot(t,vsoma,'k',t,vaxon,'r')
        xlim([t1 t2])
        xlabel('Time (ms)')
        ylabel('Membrane potential (mV)')

        subplot(313)
        plot(t,ina_syn_result,'c',t,ik_syn_result,'y',t,ina_leak_result,'b',t,ik_leak_result,'r',t,ik_klt_result,'m')
        axis([t1 t2 -6 3])
        set(gca,'ytick',-6:3:6)
        xlabel('Time (ms)')
        ylabel('Membrae currents (nA)')

        fprintf('\n')
        if k_EPSG<EPSG_iters || k_ITD<length(ITD_range)
            disp('Press a key for the next simulation')
            pause
        end
    end
end
disp('Finished')
