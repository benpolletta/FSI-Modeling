T0 = 250;  %250 millisecond long trial
dt = .005;
T = floor(T0/dt);
t = (1:T)*dt;
no_inputs = 2200; %220 = number of input synapses per cell in Hjorth et al, times 10 cells
e_rate = 2; %presynaptic firing rate (Hz) in Hjorth et al
tau_i1 = 1; tau_ir = 0.5; tau_id = 5; tau_i = 10; tau_r = 1;
tau_e1 = 1; tau_er = 0.5; tau_ed = 2;
no_cells = 10;
delta_t = 5; %number of milliseconds between two spikes to consider them "synchronous"
CE = repmat(eye(no_cells), 1, 220);    %Define connectivity from inputs to cells.

firing_rate = zeros(11,10);
spike_pairs = zeros(11, 10);
firing_avg = zeros(11,1);
pair_avg = zeros(11, 10);
spike_indicator = zeros(10,(T0/dt)-1);
for k = 1:11
	for j = 1:10
		spikes = rand(no_inputs,length(t));
		spikes = spikes < e_rate*dt/1000;

		% EPSP for spikes at time t = 0.
		epsp = tau_i*(exp(-max(t - tau_e1,0)/tau_ed) - exp(-max(t - tau_e1,0)/tau_er))/(tau_ed - tau_er);
		epsp = epsp(epsp > eps);    %?
		epsp = [zeros(1,length(epsp)) epsp]; %?

		spike_arrivals = CE*spikes; % Calculating presynaptic spikes for each cell.

		epsps = nan(size(spike_arrivals)); % Calculating EPSP experienced by each cell. %?
		for c = 1:no_cells
			epsps(c,:) = conv(spike_arrivals(c,:),epsp,'same');  %?
		end
		[Vs,Vd,s,m,h,n,t] = ing_w_dendritic_gap_jxn(10, epsps, 250, [], zeros(10), (.1*(k-1))*(rand(10) > .33));
		firing_rate(k,j) = 0;
		for a = 1:10
			Vs_pos = Vs > 0;
			Vs_sign_change = diff(Vs_pos(a,:), [], 2);
			spike_indicator(a,:) = Vs_sign_change == 1;
			firing_rate(k,j) = firing_rate(k,j) + sum(spike_indicator(a,:));
			for b = 1:10
				for d = 1:((T0/dt)-(delta_t/dt))
					for e = 1:(delta_t/dt)
						if spike_indicator(a,d) && spike_indicator(b,d+e)
                            spike_pairs(k,j) = spike_pairs(k,j) + 1;
                        end
					end
				end
			end
        end
        firing_rate(k,j) = firing_rate(k,j)/10
    end
    firing_avg = sum(firing_rate,2)/10
	pair_avg = sum(spike_pairs,2)/10
end