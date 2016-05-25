function results = Gutfreund_plots(data, results)

varied = data(1).simulator_options.vary;

no_varied = size(varied, 1);

[vary_labels, vary_vectors] = deal(cell(no_varied));

for v = 1:no_varied
    
    vary_labels{v} = varied{v, 2};
    
    vary_vectors{v} = varied{v, 3};
    
    vary_lengths(v) = length(varied{v, 3});
    
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

if isempty(results)
    
    results = AnalyzeStudy(data, @gutfreund_metrics);
    
end

result_fields = fieldnames(results);

[r, c] = subplot_size(length(result_fields));

nz_amp = nan(rfp_lengths(2), rfp_lengths(1));

for f = length(result_fields):-1:1
    
    eval(sprintf('results_for_plot = [results(:).%s];', result_fields{f}))
    
    results_for_plot = reshape(results_for_plot, fliplr(rfp_lengths));
    
    for p = 1:no_figures
        
        rfp_p = reshape(results_for_plot(p, :, :), rfp_lengths([2 1]));
        
        if f == length(result_fields)
            
            nz_amp(:, :, p) = rfp_p > .1;
            
        end
        
        rfp_p(nz_amp(:, :, p) == 0) = nan;
        
        figure(p)

        subplot(r, c, f)
        
        imagesc(vary_vectors{1}, vary_vectors{2}, rfp_p)
        
        colorbar
        
        xlabel(vary_labels{1}, 'Interpreter', 'none')
        
        ylabel(vary_labels{2}, 'Interpreter', 'none')
        
        title(result_fields{f}, 'Interpreter', 'none')
        
    end

end

for p = 1:no_figures
    
    vary_title = '';
    
    vary_name = '';
    
    for v = 1:2
        
        vary_name = [vary_name, sprintf('_%s_%fto%f', vary_labels{v}, vary_vectors{v}([1 end]))];
        
    end
    
    for v = 1:(no_varied - 2)
        
        vary_title = [vary_title, sprintf('%s = %f ', vary_labels{v + 2}, figure_params(p, v))];
        
        vary_name = [vary_name, sprintf('_%s%0.2e', vary_labels{v + 2}, figure_params(p, v))];
        
    end
    
    figure(p)
    
    mtit(p, vary_title, 'FontSize', 14, 'yoff', .2)
    
    save_as_pdf(p, ['gutfreund', vary_name])
    
end

