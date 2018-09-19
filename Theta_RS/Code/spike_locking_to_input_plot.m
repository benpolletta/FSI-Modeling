function results = spikElocking_to_input_plot(data, results, name, varargin)

v_pop = 'pop1'; i_pop = 'pop1'; label = '';

if ~isempty(varargin)
    
    for v = 1:(length(varargin)/2)
        
        if strcmp(varargin{2*v - 1}, 'i_pop')
            
            i_pop = varargin{2*v};
            
            label = [label, '_', i_pop, 'i'];
        
        elseif strcmp(varargin{2*v - 1}, 'v_pop')
            
            v_pop = varargin{2*v};
            
            label = [label, '_', v_pop, 'v'];
            
        end
        
    end
    
end
            
% input = [i_pop, '_iPeriodicPulses_input'];
% 
% voltage = [v_pop, '_v'];

close('all')

if isempty(results)
    
    results = AnalyzeStudy(data, @phase_metrics, 'v_pop', v_pop, 'i_pop', i_pop);
    
end

no_studies = length(results);

% no_freqs = size(results(1).v_spike_phases, 3);
% 
% no_periods = size(results(1).v_spike_phases, 2);

% vsp_index = 0;

% t = data(1).time;

% no_input_freqs = length(unique([data(:).pop1_PPfreq]));

varied = data(1).varied;

no_varied = length(varied);

[vary_labels, vary_vectors] = deal(cell(no_varied, 1));

for variable = 1:no_varied
    
    vary_labels{variable} = varied{variable};
    
    vary_vectors{variable} = unique([data.(varied{variable})]);
    
    vary_lengths(variable) = length(vary_vectors{variable});
    
end

[vary_lengths, vary_permute] = sort(vary_lengths, 'descend');

vary_labels = vary_labels(vary_permute);

vary_vectors = vary_vectors(vary_permute);

vary_labels(vary_lengths <= 1) = ''; 

vary_vectors(vary_lengths <= 1) = [];

vary_lengths(vary_lengths <= 1) = [];

if vary_lengths(1) <= 10 && vary_lengths(2) <= 10

    no_cols = vary_lengths(1); no_rows = vary_lengths(2);
    
else
    
    [no_rows, no_cols] = subplot_size(vary_lengths(1));
    
    vary_labels(3:(end + 1)) = vary_labels(2:end);
    
    vary_labels(1:2) = vary_labels(1);
    
    vary_vectors(3:(end + 1)) = vary_vectors(2:end);
    
    vary_vectors(1:2) = vary_vectors(1);
    
    vary_lengths(3:(end + 1)) = vary_lengths(2:end);
    
    vary_lengths(1:2) = vary_lengths(1);

end

no_varied = length(vary_lengths) - 2;

no_figures = prod(vary_lengths(3:end));

[peak_freqs, v_mean_spike_mrvs, no_spikes] = deal(nan(no_rows, no_cols, no_figures));

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

figurElabels = cell(no_figures, 1);

nonempty_plots = zeros(no_rows, no_cols, no_figures);

for f = 1:no_figures
    
    figure(f)
    
    figure_index = ones(1, length(data));
    
    figurElabels{f} = 'Spike Locking to Input';
    
    for v = 1:no_varied
       
        figure_index = figure_index & ([data.(vary_labels{v + 2})] == figure_params(f, v));
        
        if v == 1
        
            figurElabels{f} = [vary_labels{v + 2}, ' = ', num2str(figure_params(f, v), '%.3g')];
            
        else
        
            figurElabels{f} = [figurElabels{f}, '; ', vary_labels{v + 2}, ' = ', num2str(figure_params(f, v), '%.3g')];
            
        end
        
    end
    
    for r = 1:no_rows
        
        for c = 1:no_cols
                
            s = (r - 1)*no_cols + c;
            
            if strcmp(vary_labels{1}, vary_labels{2})
                
                if s <= vary_lengths(1)
                    
                    study_index = figure_index & ([data.(vary_labels{1})] == vary_vectors{1}(s));
                    
                    study_label = [vary_labels{1}, ' = ', num2str(vary_vectors{1}(s), '%.3g')];
                    
                else
                    
                    study_index = []; % zeros(size(figure_index));
                    
                    study_label = '';
                    
                end
                
            else
                
                row_index = figure_index & ([data.(vary_labels{2})] == vary_vectors{2}(r));
                
                study_index = row_index & ([data.(vary_labels{1})] == vary_vectors{1}(c));
                
                study_label = [vary_labels{1}, ' = ', num2str(vary_vectors{1}(c), '%.3g'),...
                    ', ', vary_labels{2}, ' = ', num2str(vary_vectors{2}(r), '%.3g')];
                
            end
            
            if ~isempty(study_index) % study_index = logical(study_index);
                
                nonempty_plots(r, c, f) = 1;
                
                peak_freqs(r, c, f) = results(study_index).peak_freq;
                
                % no_spikes = size(results(study_index).v_spike_phases, 1);
                
                % v_spike_phases(vsp_index + (1:no_spikes), :, :) = results(s).v_spike_phases;
                
                mean_spike_mrvs = nanmean(exp(sqrt(-1)*results(study_index).v_spike_phases));
                
                no_spikes(r, c, f) = size(results(study_index).v_spike_phases, 1);
                
                v_mean_spike_mrvs(r, c, f) = mean_spike_mrvs(1); % circ_r(results(s).v_spike_phases); %
                
                % f_index = mod(s - 1, no_input_freqs) + 1;
                %
                % o_index = ceil(s/no_input_freqs);
                
                % subplot_index = (r - 1)*no_cols + c; % no_other_conditions*(f_index - 1) + o_index;
                
                ax(r, c, f) = subplot(no_rows, no_cols, s); % no_input_freqs, no_other_conditions, (f_index - 1)*no_other_conditions + o_index)
                
                rose(gca, results(study_index).v_spike_phases(:, 1, 1)) % , 60)
                
                hold on
                
                title(study_label, 'interpreter', 'none')
                
                % if ~strcmp(vary_labels{1}, vary_labels{2})
                %
                %     ylabel([num2str(data(study_index).(vary_labels{2}), '%.3g'), ' ', vary_labels{2}])
                %
                % end
                
                % v_mean_spike_phases(s, :, :) = circ_mean(results(s).v_spike_phases);
                
                % vsp_index = vsp_index + no_spikes;
                
            end
            
        end
    
    end
    
end

% linkaxes(reshape(ax(:, :, f), no_rows*no_cols, 1))

for f = 1:no_figures
    
    figure(f)
    
    for r = 1:no_rows
        
        for c = 1:no_cols
            
            s = (r - 1)*no_cols + c;
            
            if nonempty_plots(r, c, f)
                
                subplot(no_rows, no_cols, s)
                
                multiplier = max(max(xlim), max(ylim));
                
                compass(multiplier*real(v_mean_spike_mrvs(r, c, f)), multiplier*imag(v_mean_spike_mrvs(r, c, f)), 'k')
                
            end
            
        end
        
    end
    
    mtit(gcf, figurElabels{f}, 'FontSize', 16)

    if no_figures > 1
    
        save_as_pdf(gcf, [name, label, '_rose_', num2str(f)])
    
    else
        
        save_as_pdf(gcf, [name, label, i_label, '_rose'])
        
    end
    
end

% if strcmp(vary_labels{1}, vary_labels{2})
%     
%     v_mean_spike_mrvs = reshape(v_mean_spike_mrvs, no_rows*no_cols, no_figures);
%     
% end

if strcmp(vary_labels{1}, vary_labels{2})
    
    mrv_for_plot = reshape(permute(v_mean_spike_mrvs, [2 1 3:length(size(v_mean_spike_mrvs))]), no_rows*no_cols, no_figures);
    
    mrv_for_plot = mrv_for_plot(1:vary_lengths(1), :);
    
    nspikes_for_plot = reshape(permute(no_spikes, [2 1 3:length(size(no_spikes))]), no_rows*no_cols, no_figures);
    
    nspikes_for_plot = nspikes_for_plot(1:vary_lengths(1), :);
    
else
    
    mrv_for_plot = v_mean_spike_mrvs;
    
    nspikes_for_plot = no_spikes;
    
end

no_figures = size(mrv_for_plot, 3);

for f = 1:no_figures

    figure
    
    subplot(3, 1, 1)
    
    plot(vary_vectors{1}', abs(mrv_for_plot(:, :, f))')
    
    axis tight
    
    box off
    
    ylim([0 1])
    
    title(['Spike PLV to Input by ', vary_labels{1}], 'FontSize', 16, 'interpreter', 'none')
    
    xlabel(vary_labels{1}, 'FontSize', 14, 'interpreter', 'none')
    
    ylabel('Spike PLV', 'FontSize', 14)
    
    legend(figurElabels)
    
    subplot(3, 1, 2)
    
    adjusted_plv = ((abs(mrv_for_plot(:, :, f)).^2).*nspikes_for_plot(:, :, f) - 1)./(nspikes_for_plot(:, :, f) - 1);
    
    plot(vary_vectors{1}, adjusted_plv') % (unique([data(:).(vary_labels{1})])', adjusted_plv')
    
    axis tight
    
    box off
    
    ylim([0 1])
    
    title(['Adjusted Spike PLV to Input by ', vary_labels{1}], 'FontSize', 16, 'interpreter', 'none')
    
    xlabel(vary_labels{1}, 'FontSize', 14, 'interpreter', 'none')
    
    ylabel('Spike PLV', 'FontSize', 14)
    
    subplot(3, 1, 3)
    
    plot(vary_vectors{1}, nspikes_for_plot(:, :, f)')
    
    axis tight
    
    box off
    
    title(['Number of Spikes by ', vary_labels{1}], 'FontSize', 16, 'interpreter', 'none')
    
    xlabel(vary_labels{1}, 'FontSize', 14, 'interpreter', 'none')
    
    ylabel('Number of Spikes', 'FontSize', 14)
    
    if no_figures > 1
        
        save_as_pdf(gcf, [name, label, '_MRV_', num2str(f)])
        
    else
        
        save_as_pdf(gcf, [name, label, '_MRV'])
        
    end
    
end


