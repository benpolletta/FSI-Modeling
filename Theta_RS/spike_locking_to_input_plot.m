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

no_input_freqs = length(unique([data(:).pop1_PPfreq]));

varied = data(1).simulator_options.vary;

no_varied = size(varied, 1);

[vary_labels, vary_vectors] = deal(cell(no_varied));

for variable = 1:no_varied
    
    vary_labels{variable} = varied{variable, 2};
    
    vary_vectors{variable} = varied{variable, 3};
    
    vary_lengths(variable) = length(varied{variable, 3});
    
end

[vary_lengths, vary_permute] = sort(vary_lengths, 'descend');

vary_labels = vary_labels{vary_permute};

vary_vectors = vary_vectors{vary_permute};

vary_lengths(vary_lengths <= 1) = []; vary_labels{vary_lengths <= 1} = {}; vary_vectors{vary_lengths <= 1} = {};

no_varied = length(vary_lengths) - 2;

no_rows = vary_lengths(1); no_cols = vary_lengths(2); no_figures = prod(vary_lengths(3:end));

if no_figures > 1
    
    vary_lengths_cp = cumprod(vary_lengths(3:end));
    
    figure_params = nan(no_figures, no_varied);
    
    for f = 1:no_figures
        
        figure_params(f, 1) = vary_vectors{3}(mod(f - 1, vary_lengths(3)) + 1);
        
        for v = 2:no_varied
            
            figure_params(f, v) = vary_vectors{v + 2}(ceil(f/vary_lengths_cp(v - 1)));
            
        end
        
    end
    
else
    
    no_figures = 1;
    
end

% no_other_conditions = no_studies/no_input_freqs;
% 
% if no_other_conditions == 1
% 
%     [r, c] = subplot_size(no_studies);
%     
% else
%    
%     r = no_input_freqs;
%     
%     c = no_other_conditions;
%     
% end

for f = 1:no_figures
    
    figure
    
    figure_index = ones(length(data), 1);
    
    for v = 1:no_varied
       
        figure_index = figure_index & ([getfield(data(:), vary_labels{v + 2})] == figure_params(f, v));
        
    end
    
    for r = 1:no_rows % for s = 1:no_studies
        
        row_index = figure_index & ([getfield(data(:), vary_labels{1})] == vary_vectors{1}(r));
        
        for c = 1:no_cols
        
            study_index = row_index & ([getfield(data(:), vary_labels{2})] == vary_vectors{2}(c));
            
            % v(:, s) = getfield(data(s), 'pop1_v');
            
            peak_freqs(r, c, f) = results(study_index).peak_freq;
            
            no_spikes = size(results(study_index).v_spike_phases, 1);
            
            % v_spike_phases(vsp_index + (1:no_spikes), :, :) = results(s).v_spike_phases;
            
            v_mean_spike_mrvs(r, c, f) = nanmean(exp(sqrt(-1)*results(study_index).v_spike_phases)); % circ_r(results(s).v_spike_phases); %
            
            % f_index = mod(s - 1, no_input_freqs) + 1;
            % 
            % o_index = ceil(s/no_input_freqs);
             
            subplot_index = (r - 1)*no_cols + c; % no_other_conditions*(f_index - 1) + o_index;
            
            ax(r, c, f) = subplot(r, c, subplot_index); % no_input_freqs, no_other_conditions, (f_index - 1)*no_other_conditions + o_index)
            
            rose(gca, results(study_index).v_spike_phases(:, 1, 1)) % , 60)
            
            hold on
            
            title(sprintf('%.3g Hz Spike Phases', data(study_index).pop1_PPfreq))
            
            % v_mean_spike_phases(s, :, :) = circ_mean(results(s).v_spike_phases);
            
            vsp_index = vsp_index + no_spikes;
            
        end
        
        linkaxes(reshape(ax(:, :, f), no_rows*no_cols, 1))
        
    end
    
end

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


