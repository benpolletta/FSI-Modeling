clear;
T0 = 1000;  %1 second long trial
dt = .005;
T = floor(T0/dt);
t = (1:T)*dt;
no_cells = 10;

no_e_inputs = 127*no_cells; %127 = number of AMPA input synapses per cell in Hjorth et al, times 10 cells
e_rate = 2; % presynaptic firing rate (Hz) in Hjorth et al

no_i_inputs = 93*no_cells; % = 93 times 10
i_rate = 2;

tau_i1 = 1; tau_ir = 0.5; tau_id = 5; tau_i = 10; tau_r = 1;
tau_e1 = 1; tau_er = 0.5; tau_ed = 2;
delta_t = 5; %number of milliseconds between two spikes to consider them "synchronous"
CE_e = repmat(eye(no_cells), 1, no_e_inputs/no_cells);
CE_i = repmat(eye(no_cells), 1, no_i_inputs/no_cells);    %Define connectivity from inputs to cells.

max_k = 11;

firing_rate = zeros(max_k,10);
spike_pairs = zeros(max_k, 10);
firing_avg = zeros(max_k,1);
pair_avg = zeros(max_k, 10);
spike_indicator = zeros(10,(T0/dt)-1);

[r,c] = subplot_size(max_k);

for j = 1:10 % Computing ten trials (for each level of gap junction conductance).

    e_spikes = rand(no_e_inputs,length(t));
    e_spikes = e_spikes < e_rate*dt/1000;
    
	% EPSP for spikes at time t = 0.
	epsp = tau_i*(exp(-max(t - tau_e1,0)/tau_ed) - exp(-max(t - tau_e1,0)/tau_er))/(tau_ed - tau_er);
	epsp = epsp(epsp > eps);    %?
	epsp = [zeros(1,length(epsp)) epsp]; %?
    
    e_spike_arrivals = CE_e*e_spikes; % Calculating presynaptic spikes for each cell.

    epsps = nan(size(e_spike_arrivals)); % Calculating EPSP experienced by each cell.
    for c = 1:no_cells
      epsps(c,:) = conv(e_spike_arrivals(c,:),epsp,'same');
    end
    
    i_spikes = rand(no_i_inputs,length(t));
    i_spikes = i_spikes < i_rate*dt/1000;
    
    % IPSP for spikes at time t = 0.
    ipsp = tau_i*(exp(-max(t - tau_i1,0)/tau_id) - exp(-max(t - tau_i1,0)/tau_ir))/(tau_id - tau_ir);
    ipsp = ipsp(ipsp > eps);
    ipsp = [zeros(1,length(ipsp)) ipsp];

    i_spike_arrivals = CE_i*i_spikes; % Calculating presynaptic spikes for each cell.

    ipsps = nan(size(i_spike_arrivals)); % Calculating IPSP experienced by each cell.

    for c = 1:no_cells
        ipsps(c,:) = conv(i_spike_arrivals(c,:),ipsp,'same');
    end

    CG = rand(10) > .33;
    
    figure
    
    for k = 1:max_k	

		[Vs,Vd,s,m,h,n,t] = ing_w_dendritic_synapses(no_cells, epsps-ipsps, T0, [], zeros(no_cells), (k - 1)*CG);
		firing_rate(k,j) = 0;
        
        subplot(r, c, k)
        
        plot(t,Vs'), hold on
        
		for a = 1:no_cells
            
			Vs_pos = Vs > 0;
			Vs_sign_change = diff(Vs_pos(a,:), [], 2);
			spike_indicator(a,:) = Vs_sign_change == 1;
			firing_rate(k,j) = firing_rate(k,j) + sum(spike_indicator(a,:));
            synch_indicator = sum(spike_indicator);
			
            for d = 1:((T0/dt)-(delta_t/dt)-1)
                
                synch = sum(synch_indicator(d:d+(delta_t/dt)));
                
                if synch > 0
                
                    synch = synch-1;
                
                end
                
                spike_pairs(k,j) = spike_pairs(k,j) + synch;
            
                synch_time_indicator = t <= d + (delta_t/dt) & t >= d;
                
                plot(t(synch_time_indicator)', Vs(:, synch_time_indicator)')
                
            end
            
        end
        
        display(sprintf('k = %f, firing rate = %f', k, firing_rate(k,j)/10))
        
        firing_rate(k,j) = firing_rate(k,j)/10;
    
    end
    
    firing_avg = sum(firing_rate,2)/10;
	pair_avg = sum(spike_pairs,2)/10;

end