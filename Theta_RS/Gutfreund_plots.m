function results = Gutfreund_plots(data)

results = AnalyzeStudy(data, @gutfreund_metrics);

varied = data(1).simulator_options.vary;

for v = 1:size(varied, 1)
    
    vary_labels{v} = varied{v, 2};
    
    vary_vectors{v} = varied{v, 3};
    
end

result_fields = fieldnames(results);

[r, c] = subplot_size(length(result_fields));

figure

for f = 1:length(result_fields)

    subplot(r, c, f)
    
    eval(sprintf('results_for_plot = [results(:).%s];', result_fields{f}))
    
    imagesc(vary_vectors{1}, vary_vectors{2}, reshape(results_for_plot, length(vary_vectors{2}), length(vary_vectors{1})))
    
    colorbar
    
    xlabel(vary_labels{1})
    
    ylabel(vary_labels{2})
    
    title(result_fields{f})

end

