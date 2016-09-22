function results = spike_phase_plots(data, results, variables)

if ~iscell(variables), variables = {variables}; end

if isempty(results)
    
    for var = 1:length(variables)
        
        results{var} = AnalyzeStudy(data, @phase_metrics, 'variable', variables{var});
        
    end
    
end

no_studies = length(results);

no_freqs = size(results(1).v_spike_phases, 3);

no_periods = size(results(1).v_spike_phases, 2);

vsp_index = 0;

t = data(1).time;

v_mean_spike_phases = nan(no_studies, no_periods, no_freqs);
    
peak_freqs = nan(no_studies, 1);

for s = 1:no_studies
    
    v(:, s) = getfield(data(s), variables{1});
    
    peak_freqs(s) = results(s).peak_freq;
    
    no_spikes = size(results(s).v_spike_phases, 1);
    
    v_spike_phases(vsp_index + (1:no_spikes), :, :) = results(s).v_spike_phases;
    
    v_mean_spike_phases(s, :, :) = circ_mean(results(s).v_spike_phases);
   
    vsp_index = vsp_index + no_spikes;
    
end

colors = {[1 1 0], [0 1 1], [1 0 1]};
    
period_labels = {'st', 'nd', 'rd', 'th'};

freq_labels = {'1.5 Hz', '4.5 Hz', [num2str(mean(peak_freqs), '%.2g'), ' \pm ', num2str(std(peak_freqs), '%.2g')]};

for f = 1:no_freqs
    
    subplot(no_periods + 1, 1, 1)
    
    plot(t/1000, v), xlim([t(1) t(end)]/1000)
    
    for p = 1:no_periods
        
        subplot(no_periods + 1, 2*no_freqs, p*2*no_freqs + 2*f - 1)
        
        h = rose(gca, v_spike_phases(:, p, f));
        
        set(h, 'LineWidth', 2, 'Color', colors{p})
        
        % x = get(h, 'XData'); y = get(h, 'YData');
        %
        % patch(x, y, colors{p})
        
        if f == 1
        
            ylabel(sprintf('%d%s 1/%d of Data', p, period_labels{min(p, 4)}, no_periods), 'FontSize', 12)
        
        end
            
        if p == 1, title({[freq_labels{f}, ' Spike Phases']; '(All Spikes)'}, 'FontSize', 12), end
        
        subplot(no_periods + 1, 2*no_freqs, p*2*no_freqs + 2*f)
        
        h = rose(gca, v_mean_spike_phases(:, p, f));
        
        set(h, 'LineWidth', 2, 'Color', colors{p})
        
        % x = get(h, 'XData'); y = get(h, 'YData');
        %
        % patch(x, y, colors{p})
            
        if p == 1, title('(Mean Spike Phases)', 'FontSize', 12), end
        
        % hold on
        
    end
    
end


