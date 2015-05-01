clear;
T0 = 5000;  %3000 millisecond long trial
dt = .005; %this is in milliseconds
T = floor(T0/dt);
t = (1:T)*dt;
no_cells = 2;
p_gj = 1; %probability of gap junction between any pair of cells. should be 1 if you're using 2 cells
p_inhib = 1; %probability of inhibitory connection between any pair of cells

no_e_inputs = 127*no_cells; %127 = number of AMPA input synapses per cell in Hjorth et al, times 10 cells
e_size = 0.0053;  %changes magnitude of input. 0.0053 gets you between 5 and 2 Hz firing rate.
%you stop getting firing rate decreases as conductance increases for 10 cells around 0.015. 
%at around 0.5 you can get increased spike pairs with increased conductance in 10 cell networks, but firing rate is deeply weird

no_i_inputs = 93*no_cells; % = 93 times 10
i_size = 0.0053;

inhib_strength = (5*i_size)/sqrt(no_cells); %change 5 to 70 for large

tau_i1 = 1; tau_ir = 0.5; tau_id = 5; tau_i = 10; tau_r = 1;
tau_e1 = 1; tau_er = 0.5; tau_ed = 2;
delta_t = 5; %number of milliseconds between two spikes to consider them "synchronous"

max_j = 5; %number of trials to run
max_k = 10; %number of correlation values to use

gap_conductance = (70*e_size/sqrt(no_cells));
%rhythmicity = zeros(10,10,20000,10);

firing_rate = zeros(max_k,max_j);
spike_pairs = zeros(max_k, max_j);
firing_avg = zeros(max_k,1);
pair_avg = zeros(max_k, max_j);
spike_indicator = zeros(no_cells,(T0/dt)-1);

% EPSP for spikes at time t = 0.
epsp = tau_i*(exp(-max(t - tau_e1,0)/tau_ed) - exp(-max(t - tau_e1,0)/tau_er))/(tau_ed - tau_er);
epsp = epsp(epsp > eps);    %?
epsp = [zeros(1,length(epsp)) epsp]; %?

% IPSP for spikes at time t = 0.
ipsp = tau_i*(exp(-max(t - tau_i1,0)/tau_id) - exp(-max(t - tau_i1,0)/tau_ir))/(tau_id - tau_ir);
ipsp = ipsp(ipsp > eps);
ipsp = [zeros(1,length(ipsp)) ipsp];

%rectangle for convolutions
rect = 3*ones(1,5/dt);
rect = rect(rect > eps);
rect = [zeros(1,length(rect)) rect];

Vs_traces = zeros(max_j, max_k,no_cells,T);
Vd_traces = zeros(max_j, max_k,no_cells,T);

for j = 1:max_j
    
    for k = 1:max_k
   
		fraction_shared = (1/(max_k-1))*(k-1);
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		e_rate = 2; % presynaptic firing rate (Hz) in Hjorth et al
		e_spikes = rand(floor((1-fraction_shared)*(no_e_inputs/no_cells))*no_cells,length(t));
		e_spikes = e_spikes < e_rate*dt/1000;
	
		%Define connectivity from inputs to cells.
		CE_e = repmat(eye(no_cells), 1, floor((1-fraction_shared)*(no_e_inputs/no_cells))) ;
		e_spike_arrivals = CE_e*e_spikes; % Calculating presynaptic spikes for each cell.
		epsps = nan(size(e_spike_arrivals)); % Calculating EPSP experienced by each cell.
		if ~isempty(epsps)
			for c = 1:no_cells
				epsps(c,:) = e_size*conv(e_spike_arrivals(c,:),epsp,'same');
			end
		else
			epsps(c,:) = zeros(no_cells,T);
		end
		
		i_rate = 2;
		i_spikes = rand(floor((1-fraction_shared)*(no_i_inputs/no_cells))*no_cells,length(t));
		i_spikes = i_spikes < i_rate*dt/1000;
	
		%Define connectivity from inputs to cells.
		CE_i = repmat(eye(no_cells), 1, floor((1-fraction_shared)*(no_i_inputs/no_cells)));   
		i_spike_arrivals = CE_i*i_spikes; % Calculating presynaptic spikes for each cell.
		ipsps = nan(size(i_spike_arrivals)); % Calculating IPSP experienced by each cell.
		if ~isempty(ipsps)
			for c = 1:no_cells
				ipsps(c,:) = i_size*conv(i_spike_arrivals(c,:),ipsp,'same');
			end
		else
			ipsps(c,:) = zeros(no_cells,T);
		end
	
		%shared input
		CE_e_shared = ones(no_cells, floor((fraction_shared)*(no_e_inputs/no_cells)));
		shared_e_spikes = rand(floor((fraction_shared)*(no_e_inputs/no_cells)),length(t));
		shared_e_spikes  = shared_e_spikes < e_rate*dt/1000;
		shared_e_arrivals = CE_e_shared*shared_e_spikes;
		shared_e_input = nan(size(shared_e_arrivals)); % Calculating EPSP experienced by each cell. %?
		if ~isempty(shared_e_input)
			for c = 1:no_cells
				shared_e_input(c,:)= e_size*conv(shared_e_arrivals(c,:),epsp,'same');  %?
			end
		else
			shared_e_input = zeros(no_cells,T);
		end
    
		CE_i_shared = ones(no_cells, floor((fraction_shared)*(no_i_inputs/no_cells)));
		shared_i_spikes = rand(floor((fraction_shared)*(no_i_inputs/no_cells)),length(t));
		shared_i_spikes  = shared_i_spikes < i_rate*dt/1000;
		shared_i_arrivals = CE_i_shared*shared_i_spikes;
		shared_i_input = nan(size(shared_i_arrivals)); % Calculating IPSP experienced by each cell. %?
		if ~isempty(shared_e_input)
			for c = 1:no_cells
				shared_i_input(c,:)= i_size*conv(shared_i_arrivals(c,:),ipsp,'same');  %?
			end
		else
			shared_i_input = zeros(no_cells,T);
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		%for i = 1:no_cells, epsps(i, :) = fraction_shared*shared_e_input + (1-fraction_shared)*epsps(i,:); end
        %for i = 1:no_cells, ipsps(i, :) = fraction_shared*shared_i_input + (1-fraction_shared)*ipsps(i,:); end
		
		for i = 1:no_cells, epsps(i, :) = [shared_e_input(1,:) + epsps(i,:)]; end
        for i = 1:no_cells, ipsps(i, :) = [shared_i_input(1,:) + ipsps(i,:)]; end
    
		CG = gap_conductance*(rand(no_cells) < p_gj);
		CG(logical(eye(size(CG)))) = 0;
		CS = inhib_strength*(rand(no_cells) < p_inhib);
		CS(logical(eye(size(CS)))) = 0;
		[Vs,Vd,s,m,h,n,t] = ing_w_dendritic_gap_jxn(no_cells, epsps-ipsps, T0, [], CS, CG);
		Vs_traces(j,k,:,:) = Vs;
        Vd_traces(j,k,:,:) = Vd;
		firing_rate(k,j) = 0;

		for a = 1:no_cells
			Vs_pos = Vs > 0;
			Vs_sign_change = diff(Vs_pos(a,:), [], 2);
			spike_indicator(a,:) = Vs_sign_change == 1;
			firing_rate(k,j) = firing_rate(k,j) + sum(spike_indicator(a,:));
        end
        firing_rate(k,j) = (firing_rate(k,j)/no_cells)*(1000/T0)
		
		wide_spikes = zeros(no_cells,length(t)-1);
		for c = 1:no_cells
			wide_spikes(c,:) = conv(spike_indicator(c,:),rect,'same');
		end
			
		spike_pairs(k,j) = 0;
		for b = 1:no_cells
			for c = 1:no_cells
				if b ~= c %&& CS(b,c) > 0
					wide_sum = wide_spikes(b,:) + wide_spikes(c,:);
					foo = wide_sum > 5;
					wide_synch = diff(foo, [], 2);
					synch_indicator(b,c,:) = wide_synch == 1;
					spike_pairs(k,j) = spike_pairs(k,j) + sum(synch_indicator(b,c,:));
				end
			end
		end
    end
    firing_avg = sum(firing_rate,2)/max_j
end
 %spike_pairs = spike_pairs./(2*firing_rate);
 %pair_avg = sum(spike_pairs,2)/max_j

spike_pairs = spike_pairs./(firing_rate*(T0/1000)); %element-wise
pair_avg = sum(spike_pairs,2)/max_j

intervals = (0:(1/(max_k-1)):1);
figure
plot(intervals,firing_avg)
str = ['Average firing rate, ',num2str(no_cells), ' cells, correlated input, ' num2str(max_j), ' trials']
title(str)
xlabel('Fraction of input shared')
ylabel('Firing rate (Hz)')

figure
plot(intervals,pair_avg)
str = ['Spike pairs, ',num2str(no_cells), ' cells, correlated input, ' num2str(max_j), ' trials']
title(str)
xlabel('Fraction of input shared')
ylabel('Proportion of paired spikes')
 
%for i = 1:5
%    for j = 1:10
%        subplot(10,1,i)
%        plot(squeeze(rhythmicity(i,j,:,:)))
%    end
%end