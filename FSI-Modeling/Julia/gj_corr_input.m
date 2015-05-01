function [mean_pop_firing_rate, mean_pop_norm_spike_pairs, Vs_traces, Vd_traces, s_traces] = gj_corr_input(T0,no_cells,p_gj,max_j,smallvars,p_inhib,inhib_strength,gj_strength)
tic

dt = .005; %this is in milliseconds
T = floor(T0/dt);
t = (1:T)*dt;

no_e_inputs = 127*no_cells; %127 = number of AMPA input synapses per cell in Hjorth et al, times 10 cells
e_rate = 2; % presynaptic firing rate (Hz) in Hjorth et al
e_size = 0.0053; %changes magnitude of input. 0.0053 gets you between 5 and 2 Hz firing rate.
%you stop getting firing rate decreases as conductance increases for 10 cells around 0.015. 
%at around 0.5 you can get increased spike pairs with increased conductance in 10 cell networks, but firing rate is deeply weird

if nargin < 8
    gj_strength = (70*e_size)/sqrt(no_cells); %magnitude of steps of strength of gap junction
end

no_i_inputs = 93*no_cells; % = 93 times 10
i_rate = 2;
i_size = 0.0053;

if nargin < 7
    inhib_strength = 3*(5*i_size)/sqrt(no_cells);
end
if nargin < 6
    p_inhib = 0;
end

tau_i1 = 1; tau_ir = 0.5; tau_id = 5; tau_i = 10; tau_r = 1;
tau_e1 = 1; tau_er = 0.5; tau_ed = 2;

CE_e = repmat(eye(no_cells), 1, no_e_inputs/no_cells);
CE_i = repmat(eye(no_cells), 1, no_i_inputs/no_cells);    %Define connectivity from inputs to cells.

max_k = 11; %number of conductance values to use

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

firing_rate = zeros(no_cells, max_k, max_j);
spike_pairs = zeros(max_k, max_j);

if smallvars > 0
    Vs_traces = zeros(max_j, max_k ,no_cells,T);
    Vd_traces = zeros(max_j, max_k, no_cells,T);
    s_traces = zeros(max_j, max_k, no_cells,T);
end

params = struct('dt', dt, 't', t, 'T', T, 'no_e_inputs', no_e_inputs, 'e_rate', e_rate, 'e_size', e_size, 'gj_strength', gj_strength, 'no_i_inputs', no_i_inputs, 'i_rate', i_rate, 'i_size', i_size, 'CE_e', CE_e, 'CE_i', CE_i, 'max_k', max_k, 'epsp', epsp, 'ipsp', ipsp, 'T0', T0, 'no_cells', no_cells, 'p_gj', p_gj, 'rect', rect, 'p_inhib', p_inhib, 'inhib_strength', inhib_strength, 'smallvars', smallvars);

for j = 1:max_j
	j
    if smallvars > 0
        [firing_rate(:,:,j), spike_pairs(:,j), Vs_traces(j,:,:,:), Vd_traces(j,:,:,:), s_traces(j,:,:,:)] = trial(params);
    else
        [firing_rate(:,:,j), spike_pairs(:,j), ~,~,~] = trial(params);
    end
end

pop_firing_rate = reshape(sum(firing_rate), max_k, max_j)*(1000/(T0*no_cells));
spike_pairs
mean_pop_firing_rate = mean(pop_firing_rate, 2)
spike_pairs = spike_pairs./(2*pop_firing_rate*(T0/1000)); %element-wise
mean_pop_norm_spike_pairs = sum(spike_pairs,2)/max_j

toc

str = ['corrsmalldata', num2str(T0), '_', num2str(no_cells),'_',num2str(p_gj),'_',num2str(max_j),'_',num2str(p_inhib),'_',num2str(inhib_strength),'_',num2str(gj_strength),'.mat'];
save(str,'','-v7')
if smallvars > 0
    str = ['corrbigdata', num2str(T0), '_', num2str(no_cells),'_',num2str(p_gj),'_',num2str(max_j),'_',num2str(p_inhib),'_',num2str(inhib_strength),'_',num2str(gj_strength),'.mat'];
    save(str,'Vs_traces','Vd_traces','s_traces','-v7.3')
end

% intervals = (0:(1/(max_k-1)):1);
% figure
% plot(intervals,mean_pop_firing_rate)
% str = ['Average firing rate, ',num2str(no_cells), ' cells, correlated input, ' num2str(max_j), ' trials'];
% title(str)
% xlabel('Fraction of input shared')
% ylabel('Firing rate (Hz)')
% savefig('10_corr_fr.fig')
% 
% 
% figure
% plot(intervals,mean_pop_norm_spike_pairs)
% str = ['Spike pairs, ',num2str(no_cells), ' cells, correlated input, ' num2str(max_j), ' trials'];
% title(str)
% xlabel('Fraction of input shared')
% ylabel('Proportion of paired spikes')
% savefig('10_corr_pairs.fig')
end

function [firing_rate_jslice, norm_spike_pairs_jslice, Vs_traces_jslice, Vd_traces_jslice, s_traces_jslice] = trial(params)
	firing_rate_jslice = zeros(params.no_cells, params.max_k, 1);
	norm_spike_pairs_jslice = zeros(params.max_k, 1);
	Vs_traces_jslice = zeros(1, params.max_k, params.no_cells, params.T);
	Vd_traces_jslice = zeros(1, params.max_k, params.no_cells, params.T);
    s_traces_jslice = zeros(1, params.max_k, params.no_cells, params.T);
    
    if params.smallvars > 0
        parfor k = 1:params.max_k
        	k
            [firing_rate_jslice(:,k,1), norm_spike_pairs_jslice(k,1), Vs_traces_jslice(1,k,:,:), Vd_traces_jslice(1,k,:,:), s_traces_jslice(1,k,:,:)] = gj_value(params, k);
        end
    else
        parfor k = 1:params.max_k
            k
            [firing_rate_jslice(:,k,1), norm_spike_pairs_jslice(k,1), ~, ~, ~] = gj_value(params, k);
        end
    end
end

function [firing_rate_kslice, spike_pairs_kslice, Vs_traces_kslice, Vd_traces_kslice, s_traces_kslice] = gj_value(params, k)

	firing_rate_kslice = zeros(params.no_cells, 1, 1);
	spike_pairs_kslice = zeros(1,1);
	Vs_traces_kslice = zeros(1, 1, params.no_cells, params.T);
	Vd_traces_kslice = zeros(1, 1, params.no_cells, params.T);
    s_traces_kslice = zeros(1, 1, params.no_cells, params.T);

    spike_indicator = zeros(params.no_cells,(params.T0/params.dt)-1);
	synch_indicator = zeros(params.no_cells, params.no_cells, (params.T0/params.dt)-2);

    fraction_shared = (1/(params.max_k-1))*(k-1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    e_rate = 2; % presynaptic firing rate (Hz) in Hjorth et al
    e_spikes = rand(floor((1-fraction_shared)*(params.no_e_inputs/params.no_cells))*params.no_cells,length(params.t));
    e_spikes = e_spikes < params.e_rate*params.dt/1000;

    %Define connectivity from inputs to cells.
    CE_e = repmat(eye(params.no_cells), 1, floor((1-fraction_shared)*(params.no_e_inputs/params.no_cells))) ;
    e_spike_arrivals = CE_e*e_spikes; % Calculating presynaptic spikes for each cell.
    epsps = nan(size(e_spike_arrivals)); % Calculating EPSP experienced by each cell.
    if ~isempty(epsps)
        for c = 1:params.no_cells
            epsps(c,:) = params.e_size*conv(e_spike_arrivals(c,:),params.epsp,'same');
        end
    else
        for c = 1:params.no_cells
            epsps(c,:) = zeros(params.no_cells,params.T);
        end
    end

    i_rate = 2;
    i_spikes = rand(floor((1-fraction_shared)*(params.no_i_inputs/params.no_cells))*params.no_cells,length(params.t));
    i_spikes = i_spikes < i_rate*params.dt/1000;

    %Define connectivity from inputs to cells.
    CE_i = repmat(eye(params.no_cells), 1, floor((1-fraction_shared)*(params.no_i_inputs/params.no_cells)));   
    i_spike_arrivals = CE_i*i_spikes; % Calculating presynaptic spikes for each cell.
    ipsps = nan(size(i_spike_arrivals)); % Calculating IPSP experienced by each cell.
    if ~isempty(ipsps)
        for c = 1:params.no_cells
            ipsps(c,:) = params.i_size*conv(i_spike_arrivals(c,:),params.ipsp,'same');
        end
    else
        ipsps(c,:) = zeros(params.no_cells,params.T);
    end

    %shared input
    CE_e_shared = ones(params.no_cells, floor((fraction_shared)*(params.no_e_inputs/params.no_cells)));
    shared_e_spikes = rand(floor((fraction_shared)*(params.no_e_inputs/params.no_cells)),length(params.t));
    shared_e_spikes  = shared_e_spikes < params.e_rate*params.dt/1000;
    shared_e_arrivals = CE_e_shared*shared_e_spikes;
    shared_e_input = nan(size(shared_e_arrivals)); % Calculating EPSP experienced by each cell. 
    if ~isempty(shared_e_input)
        for c = 1:params.no_cells
            shared_e_input(c,:)= params.e_size*conv(shared_e_arrivals(c,:),params.epsp,'same'); 
        end
    else
        shared_e_input = zeros(params.no_cells,params.T);
    end

    CE_i_shared = ones(params.no_cells, floor((fraction_shared)*(params.no_i_inputs/params.no_cells)));
    shared_i_spikes = rand(floor((fraction_shared)*(params.no_i_inputs/params.no_cells)),length(params.t));
    shared_i_spikes  = shared_i_spikes < params.i_rate*params.dt/1000;
    shared_i_arrivals = CE_i_shared*shared_i_spikes;
    shared_i_input = nan(size(shared_i_arrivals)); % Calculating IPSP experienced by each cell. 
    if ~isempty(shared_e_input)
        for c = 1:params.no_cells
            shared_i_input(c,:)= params.i_size*conv(shared_i_arrivals(c,:),params.ipsp,'same');  
        end
    else
        shared_i_input = zeros(params.no_cells,params.T);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i = 1:params.no_cells, epsps(i, :) = [shared_e_input(1,:) + epsps(i,:)]; end
    for i = 1:params.no_cells, ipsps(i, :) = [shared_i_input(1,:) + ipsps(i,:)]; end

    CG = params.gj_strength*(rand(params.no_cells) < params.p_gj);
    CG = triu(CG);
    CG = CG + CG.';
    CG(logical(eye(size(CG)))) = 0;

    CS = params.inhib_strength*(rand(params.no_cells) < params.p_inhib);
    CS(logical(eye(size(CS)))) = 0;

    if params.smallvars > 0
        [Vs,Vd_traces_kslice(1,1,:,:),s_traces_kslice(1,1,:,:),~,~,~,~] = ing_w_dendritic_gap_jxn(params.no_cells, epsps-ipsps, params.T0, [], CS, (k-1)*CG);
        Vs_traces_kslice(1,1,:,:) = Vs;
    else
        [Vs,~,~,~,~,~,~] = ing_w_dendritic_gap_jxn(params.no_cells, epsps-ipsps, params.T0, [], CS, (k-1)*CG);
    end
    
	for a = 1:params.no_cells
        Vs_pos = Vs > 0;
        Vs_sign_change = diff(Vs_pos(a,:), [], 2);
        spike_indicator(a,:) = Vs_sign_change == 1;
        firing_rate_kslice(a, 1, 1) = sum(spike_indicator(a, :));
    end
		
    wide_spikes = zeros(params.no_cells,length(params.t)-1);
    for c = 1:params.no_cells
		wide_spikes(c,:) = conv(spike_indicator(c,:),params.rect,'same');
	end
			
	for b = 1:params.no_cells
		for c = 1:params.no_cells
			if CG(b,c) > 0 && b ~= c 
				wide_sum = wide_spikes(b,:) + wide_spikes(c,:);
				foo = wide_sum > 5;
				wide_synch = diff(foo, [], 2);
				synch_indicator(b,c,:) = wide_synch == 1;
				spike_pairs_kslice(1,1) = spike_pairs_kslice(1,1) + sum(synch_indicator(b,c,:));
			end
		end
	end
end