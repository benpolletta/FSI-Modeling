function results = spike_locking_to_input_plot(data, results, name)

if isempty(results)
    
    results = AnalyzeStudy(data, @phase_metrics);
    
end

no_studies = length(results);

no_freqs = size(results(1).v_spike_phases, 3);

no_periods = size(results(1).v_spike_phases, 2);

vsp_index = 0;

% t = data(1).time;

v_mean_spike_mrvs = nan(no_studies, no_periods, no_freqs);
    
peak_freqs = nan(no_studies, 1);

% no_input_freqs = length(unique([data(:).pop1_PPfreq]))
% 
% no_other_conditions = ceil(no_studies/no_input_freqs)

[r, c] = subplot_size(no_studies);

figure

for s = 1:no_studies
    
    % v(:, s) = getfield(data(s), 'pop1_v');
    
    peak_freqs(s) = results(s).peak_freq;
    
    no_spikes = size(results(s).v_spike_phases, 1);
    
    % v_spike_phases(vsp_index + (1:no_spikes), :, :) = results(s).v_spike_phases;
    
    v_mean_spike_mrvs(s, :, :) = nanmean(exp(sqrt(-1)*results(s).v_spike_phases)); % circ_r(results(s).v_spike_phases); % 
    
    % f_index = mod(s - 1, no_input_freqs) + 1
    % 
    % o_index = ceil(s/no_input_freqs)
    
    ax(s) = subplot(r, c, s); % no_input_freqs, no_other_conditions, (f_index - 1)*no_other_conditions + o_index)
    
    rose(gca, results(s).v_spike_phases(:, 1, 1)) % , 60)
    
    hold on
    
    title(sprintf('%.3g Hz Spike Phases', data(s).pop1_PPfreq))
    
    % v_mean_spike_phases(s, :, :) = circ_mean(results(s).v_spike_phases);
   
    vsp_index = vsp_index + no_spikes;
    
end

linkaxes(ax)

for s = 1:no_studies
   
    subplot(r, c, s)
    
    multiplier = max(max(xlim), max(ylim));
    
    compass(multiplier*real(v_mean_spike_mrvs(s, 1, 1)), multiplier*imag(v_mean_spike_mrvs(s, 1, 1)), 'k')
    
end

save_as_pdf(gcf, [name, '_rose'])

figure,

plot([data(:).pop1_PPfreq]', abs(v_mean_spike_mrvs(:, 1, 1)))

axis tight

ylim([0 1])

title('Spike PLV to Input by Freq.', 'FontSize', 16)

xlabel('Freq. (Hz)', 'FontSize', 14)

ylabel('Spike PLV', 'FontSize', 14)

save_as_pdf(gcf, [name, '_MRV'])


