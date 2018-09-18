function results = Gutfreund_plots(data, results, variables)

if ~iscell(variables), variables = {variables}; end

if isempty(results)
    
    for var = 1:length(variables)
        
        results{var} = AnalyzeStudy(data, @gutfreund_metrics, 'variable', variables{var});
        
    end
    
end

for var = 1:length(variables)
    
    varied = data(1).simulator_options.vary;
    
    no_varied = size(varied, 1);
    
    [vary_labels, vary_vectors] = deal(cell(no_varied));
    
    for variable = 1:no_varied
        
        vary_labels{variable} = varied{variable, 2};
        
        vary_vectors{variable} = varied{variable, 3};
        
        vary_lengths(variable) = length(varied{variable, 3});
        
    end
    
    rfp_lengths = [vary_lengths([1 2]) prod(vary_lengths(3:end))];
    
    if no_varied > 2
        
        no_figures = rfp_lengths(end);
        
        vary_lengths_cp = cumprod(vary_lengths(3:end));
        
        figure_params = nan(no_figures, no_varied - 2);
        
        for p = 1:no_figures
            
            figure_params(p, 1) = vary_vectors{3}(mod(p - 1, vary_lengths(3)) + 1);
            
            for v = 4:no_varied
                
                figure_params(p, v - 2) = vary_vectors{v}(ceil(p/vary_lengths_cp(v - 3)));
                
            end
            
        end
        
    else
        
        no_figures = 1;
        
    end
    
    for p = 1:no_figures
        
        handle(p) = figure;
        
    end
    
    result = results{var};
    
    result_fields = fieldnames(result);
    
    [r, c] = subplot_size(length(result_fields));
    
    % nz_amp = nan(rfp_lengths(2), rfp_lengths(1));
    
    for f = length(result_fields):-1:1
        
        eval(sprintf('results_for_plot = [result(:).%s];', result_fields{f}))
        
        results_for_plot = reshape(results_for_plot, fliplr(rfp_lengths));
        
        for p = 1:no_figures
            
            rfp_p = reshape(results_for_plot(p, :, :), rfp_lengths([2 1]));
            
            % if f == length(result_fields)
            % 
            %     nz_amp(:, :, p) = rfp_p > .1;
            % 
            % end
            %
            % rfp_p(nz_amp(:, :, p) == 0) = nan;
            
            figure(handle(p))
            
            subplot(r, c, f)
            
            imagesc(vary_vectors{1}, vary_vectors{2}, rfp_p)
            
            colorbar
            
            xlabel(vary_labels{1}, 'Interpreter', 'none')
            
            ylabel(vary_labels{2}, 'Interpreter', 'none')
            
            title(result_fields{f}, 'Interpreter', 'none')
            
        end
        
    end
    
    for p = 1:no_figures
        
        figure(handle(p))
        
        vary_title = '';
        
        vary_name = '';
        
        for variable = 1:2
            
            vary_name = [vary_name, sprintf('_%s_%gto%g', vary_labels{variable}, vary_vectors{variable}([1 end]))];
            
        end
        
        for variable = 1:(no_varied - 2)
            
            vary_title = [vary_title, sprintf('%s = %g ', vary_labels{variable + 2}, figure_params(p, variable))];
            
            vary_name = [vary_name, sprintf('_%s%g', vary_labels{variable + 2}, figure_params(p, variable))];
            
        end
        
        mtit(handle(p), vary_title, 'FontSize', 14, 'yoff', .2)
        
        save_as_pdf(handle(p), ['Figures/gutfreund_metrics', vary_name, '_', variables{var}])
        
    end
    
end

