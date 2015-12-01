function psps = multi_Poisson_depressing(no_cells, inputs_per_cell, rate, T, dt, offset)

t = 0:dt:T;

psps = zeros(no_cells, length(t));

for c = 1:no_cells
    
    cell_psps = zeros(length(t), inputs_per_cell);
    
    if isscalar(rate)
    
        spikes = poissrnd(rate*dt/1000, inputs_per_cell, length(t));
    
    else
        
        spikes = poissrnd(repmat(rate*dt/1000, 1, inputs_per_cell))';
        
    end
    
    for i = 1:inputs_per_cell
    
        spike_times = t(find(spikes(i, :) == 1));
        
        if ~isempty(spike_times)
            
            ISIs = diff(spike_times);
            
            amps = depressing_synapse_ISIs_to_amps(ISIs);
            
            input_psps = depressing_synapse_amps_to_PSPs(spike_times, amps, dt, 'i', offset);
            
            cell_psps(1:length(input_psps), i) = input_psps;
            
        end
        
    end
    
    psps(c, :) = sum(cell_psps');
    
end