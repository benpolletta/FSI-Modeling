T0 = 2000;  %250 millisecond long trial
dt = .005;
T = floor(T0/dt);
t = (1:T)*dt;

no_e_inputs = 1270; %197 = number of AMPA input synapses per cell in Hjorth et al, times 10 cells
e_rate = 2; % presynaptic firing rate (Hz) in Hjorth et al

no_i_inputs = 930; % = 93 times 10
i_rate = 2;

tau_i1 = 1; tau_ir = 0.5; tau_id = 5; tau_i = 10; tau_r = 1;
tau_e1 = 1; tau_er = 0.5; tau_ed = 2;
no_cells = 10;
CE_e = repmat(eye(no_cells), 1, no_e_inputs/10);
CE_i = repmat(eye(no_cells), 1, no_i_inputs/10);    %Define connectivity from inputs to cells.

% EPSP for spikes at time t = 0.
epsp = tau_i*(exp(-max(t - tau_e1,0)/tau_ed) - exp(-max(t - tau_e1,0)/tau_er))/(tau_ed - tau_er);
epsp = epsp(epsp > eps);
epsp = [zeros(1,length(epsp)) epsp];

% IPSP for spikes at time t = 0.
ipsp = tau_i*(exp(-max(t - tau_i1,0)/tau_id) - exp(-max(t - tau_i1,0)/tau_ir))/(tau_id - tau_ir);
ipsp = ipsp(ipsp > eps);
ipsp = [zeros(1,length(ipsp)) ipsp];
    
%% Generating inputs for all conductances.
e_spikes = rand(no_e_inputs,length(t));
e_spikes = e_spikes < e_rate*dt/1000;

e_spike_arrivals = CE_e*e_spikes; % Calculating presynaptic spikes for each cell.

epsps = nan(size(e_spike_arrivals)); % Calculating EPSP experienced by each cell.
for c = 1:no_cells
    epsps(c,:) = conv(e_spike_arrivals(c,:),epsp,'same');
end

i_spikes = rand(no_i_inputs,length(t));
i_spikes = i_spikes < i_rate*dt/1000;

i_spike_arrivals = CE_i*i_spikes; % Calculating presynaptic spikes for each cell.

ipsps = nan(size(i_spike_arrivals)); % Calculating EPSP experienced by each cell.

for c = 1:no_cells
    ipsps(c,:) = conv(i_spike_arrivals(c,:),ipsp,'same');
end

%% Performing simulations.

firing_rate = zeros(11,10);

figure

for k = 1:11 % Looping over gap junction strength.
    
    [Vs,Vd,s,m,h,n,t] = ing_w_dendritic_gap_jxn(10, epsps - ipsps, T0, [], zeros(10), (.1*(k-1))*(rand(10) > .4));
    
    subplot(3,4,k)
    plot(t', Vs')
    
    Vs_pos = Vs > 0;
    Vs_sign_change = diff(Vs_pos, [], 2);
    spike_indicator = Vs_sign_change == 1;
    
    firing_rate(k, :) = sum(spike_indicator, 2)/(T0/1000);
    
end

firing_avg = sum(firing_rate,2)/10;