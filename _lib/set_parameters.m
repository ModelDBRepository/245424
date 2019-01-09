%%  PARAMETERS
tend    = 1000; % 5000
dt      = 0.01; % ms
t = dt:dt:tend;

% CELL PARAMETERS
gpas    = 0.000858; % (S/cm2)
gklt    = 0.0136; % (S/cm2) peak conductance
len     = 150*1e-4; % (cm) length of branch
seglength = 5*1e-4; % compartment length (cm)
diam    = 2.5*1e-4; % (cm)
area_s  = (20*1e-4)^2*pi; % soma area cm2

tw_fac  = 1; % klt activation time constant - scaling factor

% passive parameters
cm      = 1e-6; % (F/cm2)
Ra      = 200; % (ohm*cm)
tm      = cm/gpas; % (sec)
ra      = 4*Ra./(pi*diam.^2); % (ohm/cm)
lambda  = sqrt(diam./(4*gpas*Ra)); % (cm)
RNinf   = sqrt(Ra/gpas)*(2./(pi*(diam).^(3/2))); % (ohm) : R-in of semi-infinite cable (Koch pg 33) : identical to lambda*ra
Rs      = 1/(area_s*gpas);
L       = len./lambda;
dx      = 0.01;
nseg    = floor(L/dx);
if nseg<3
    dx = L/3;
    nseg = 3;
end
dx      = L/nseg; % also cumbersome, but to make sure that length is not changed through round-off errors
rho     = Rs./RNinf; % ratio of soma resistance over Rinf of semi-infinite cable


v_init  = -60; % resting voltage
Ena     = 53;
Ek      = -106;

% KLT current parameters
whalf   = -57.34;
wk      = -11.7;
tw1     = 6; tw2 = 24; tw3 = -60; twk1 = 7; twk2 = -51; 
tadjw   = 0.22*tw_fac;
tws     = 1.59;
zhalf   = -67; zk = 6.16; h = 0.27;
w_inf   = 1/(1+exp((v_init-whalf)/wk));
z_inf   = (1-h)/(1+exp((v_init-zhalf)/zk))+h;
gamma_m = gklt/gpas;
tw      = (100/(tw1*exp((v_init-tw3)/twk1)+tw2*exp((v_init-tw3)/twk2))+tws) * tadjw * 1e-3;

El      = gamma_m*z_inf*w_inf^4*(v_init-Ek) + v_init;
gnapas  = gpas / ( 1+(Ena-El)/(El-Ek) );
gkpas   = gpas - gnapas;

sdAreaRatio = area_s/(len*diam*pi/nseg); % how much larger is soma area compared to single compartment


% AXONAL SPIKE GENERATION PARAMETERS
vth_above_v_init = 10; % how many mV above rest?
vth             = v_init+vth_above_v_init; % mV 
tref            = 1; % ms 1
tau_axon        = 0.2; % (ms)
gc_c_axon       = 1/0.05; % 1/tau (1/ms) reciprocal of tau coupling : no idea what is appropriate for this (check Josh results) 50 (=20 us)

% STIMULUS PARAMETERS
Ninputs         = 12;   % 12: total of both sides
input_freq      = 500 ; % (Hz)
refrac          = 1;    % (ms)
Rave            = 240;  % (Hz) maximum firing rate of inputs
vec_strength    = 0.988;

stim_start      = dt;   % (ms)
total_stim_dur  = tend; % (ms)
vec_factor      = sqrt(-pi^2/log(vec_strength)/2);

A_EPSG          = 0.0015; % uS
Es              = 0; % mV
taur            = 0.2;
taud            = 0.2001;
tp              = (taur*taud)/(taud - taur) * log(taud/taur); % time to peak
EPSG_factor     = 1/(-exp(-tp/taur) + exp(-tp/taud));
tendEPSG        = 8;
tEPSG           = 0:dt:tendEPSG;
EPSG_           = EPSG_factor*(exp(-tEPSG/taud)-exp(-tEPSG/taur)); % normalized to 1 --> is scaled below

% CREATE VECTOR OF INPUT LOCATIONS: UNIFORMLY SPREAD OVER DISTAL 2/3 OF DENDRITES
xvec        = [(-L+dx/2):dx:0 0 dx/2:dx:L];
nseg        = floor(L/dx); % soma at nseg
locvec = [];
for k = 1:Ninputs/2
    locvec(k)           = ceil(k*nseg*2/3/(Ninputs/2+1)) - 1;
    locvec(Ninputs-k+1) = 2*nseg - locvec(k);
end

