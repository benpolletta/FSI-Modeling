function [V, s_i, s_e, spikes, CI, CE, t] = LIF_interneuron_network(no_cells,no_inputs,T0,e_rate)
%UNTITLED2 Right-hand side of leaky-integrate-and-fire neuron, for
%integration with RK45.
%   T0 is in milliseconds.

dt = .005;

T = floor(T0/dt);

V_rest = -70; V_thresh = -52; V_reset = -59; V_spike = 0;

gE = 0.3; gI = 4;

tau_i1 = 1; tau_ir = 0.5; tau_id = 5; tau_i = 10;
tau_e1 = 1; tau_er = 0.5; tau_ed = 2; tau_e = 20;

t = (1:T)*dt;                     %Define time axis vector (useful for plotting).

% spike_times = poiss_spike_times(T0/1000,1/e_rate,no_inputs); %Compute excitatory inputs.

spikes = rand(no_inputs,length(t));
spikes = spikes < e_rate*dt/1000;

% for i = 1:length(t)
%     
%     cum_spikes(:,i) = sum(spike_times < t(i),2);
%         
% end
% 
% e_spike_times = diff([cum_spikes cum_spikes(:,end)],1,2);

CE = rand(no_cells,no_inputs);    %Define connectivity from inputs to cells.
CE = CE<.1;
CE = CE*gE;

CI = rand(no_cells,no_cells);     %Define connectivity between inhibitory cells.
CI = CI<.1;
CI = CI - diag(diag(CI));
CI = CI*gI;

V = zeros(no_cells,T);                %Make empty variables to hold I-cell results.
s_i = zeros(no_cells,T);                %Make empty variables to hold the synapse results.
s_e = zeros(no_inputs,T);

% i_spike_times = t(end)*ones(no_cells,1);
% e_spike_times = t(end)*ones(no_inputs,1);

% IPSPs for spikes at time t = 0.
ipsp = tau_i*(exp(-max(t - tau_i1,0)/tau_id) - exp(-max(t - tau_i1,0)/tau_ir))/(tau_id - tau_ir);
ipsp = repmat(ipsp,no_cells,1);

% EPSP for spikes at time t = 0.
epsp = tau_i*(exp(-max(t - tau_e1,0)/tau_ed) - exp(-max(t - tau_e1,0)/tau_er))/(tau_ed - tau_er);
epsp = repmat(epsp,no_inputs,1);

V(:,1)= -70.0 + 10*rand(no_cells,1);  	%Set the initial conditions for I-cells.

for i=1:T-1                       %Integrate the equations.
    
    V(:,i+1) = V(:,i) + dt*((V_rest - V(:,i))/tau_i - (CE*s_e(:,i)).*V(:,i) + (CI*s_i(:,i)).*(-70-V(:,i)));
    
    V(V(:,i+1) >= V_thresh,i+1) = V_spike;
    V(V(:,i) == V_spike,i+1) = V_reset;
    
%     s_i_rhs = tau_i*(-exp(-max(t(i) - i_spike_times - tau_i1,0)/tau_id)/tau_id ...
%         + exp(-max(t(i) - i_spike_times - tau_i1,0)/tau_ir)/tau_ir)/(tau_id - tau_ir) ...
%         - (i_spike_times == t(end))*tau_i/(tau_id*tau_ir);
%     
%     s_e_rhs = tau_i*(-exp(-max(t(i) - e_spike_times - tau_e1,0)/tau_ed)/tau_ed ...
%         + exp(-max(t(i) - e_spike_times - tau_e1,0)/tau_er)/tau_er)/(tau_ed - tau_er) ...
%         - (e_spike_times == t(end))*tau_i/(tau_ed*tau_er);
%     
%     s_i(:,i+1) = s_i(:,i) + dt*s_i_rhs;
%     s_e(:,i+1) = s_e(:,i) + dt*s_e_rhs;

%     s_i(i_spike_indices,:) = tau_i*(exp(-max(t - t(i+1)*ones(size(i_spike_indices)) - tau_i1,0)/tau_id) ...
%         - exp(-max(t - t(i+1)*ones(size(i_spike_indices)) - tau_i1,0)/tau_ir))/(tau_id - tau_ir);
%     
%     s_e(e_spike_indices,:) = tau_i*(exp(-max(t - t(i+1)*ones(size(e_spike_indices)) - tau_e1,0)/tau_ed) ...
%         - exp(-max(t - t(i+1)*ones(size(e_spike_indices)) - tau_e1,0)/tau_er))/(tau_ed - tau_er);

    if any(V(:,i+1) >= V_thresh)

        i_spike_indices = find(V(:,i+1) >= V_thresh);

        s_i(i_spike_indices,(i+1):end) = s_i(i_spike_indices,(i+1):end) + ipsp(i_spike_indices,1:T-i);

    end
    
    if any(spikes(:,i+1) == 1)
    
        e_spike_indices = find(spikes(:,i+1) == 1);

        s_e(e_spike_indices,(i+1):end) = s_e(e_spike_indices,(i+1):end) + epsp(e_spike_indices,1:T-i);
    
    end

end

end

