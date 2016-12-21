function make_figure_for_presentation(name)

load([name, '.mat'])

time = data(1).time;

for i = 1:20
    
    V(i, :) = data(i).pop1_v'; 

end

for i = 1:20
    
    I(i, :) = data(i).pop1_iPeriodicPulses_Iext';

end

for i = 1:20
    
    subplot(20, 1, i)
    
    [a, b, c] = plotyy(time, V(i, :)', time, -I(i, :)');
    
    box off
    
    if mod(i, 2) == 1
    
        ylabel(a(1), sprintf('%d Hz', i), 'FontSize', 16)
        
    end
    
    set(b, 'LineWidth', 2)
    
    set(c, 'Color', 'r')
    
    axis(a(1),'tight')
    
    set(a,{'ycolor'}, {'b';'r'})
    
    if i~=20
        
        set(a,'XTick',[],'XTickLabel',[])
    
    else
        
        xlabel('Time (ms)', 'FontSize', 16)
        
    end
    
    linkaxes(a,'x')

end