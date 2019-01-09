%% create_input
gvec        = zeros(Ninputs,length(t));
for input_ind = 1:Ninputs
    events_vec  = create_events(stim_start,total_stim_dur,input_freq,vec_factor,refrac,Rave);
    EPSG_ln = length(EPSG_);
    EPSG_inds = events_vec/dt;
    for k = 1:length(EPSG_inds)
        EPSG_ind_vec = round(EPSG_inds(k):(EPSG_inds(k)+EPSG_ln-1));
        if EPSG_ind_vec(end)>length(gvec)
            EPSG_ind_vec = EPSG_ind_vec(EPSG_ind_vec<=length(gvec));
        end
        gvec(input_ind,EPSG_ind_vec) = gvec(input_ind,EPSG_ind_vec)+EPSG_(1:length(EPSG_ind_vec));
    end
end
