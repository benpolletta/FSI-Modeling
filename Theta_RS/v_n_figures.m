% function v_n_figures(no_rows, no_cols)

no_rows = 5;

no_cols = 15;

figure,

for s = 1:no_rows*no_cols, 

    subplot(no_rows, no_cols, s) 
    
    plot(data(s).pop1_v, data(s).pop1_iKs_n)
    
    box off
    
    ylim([0 .8])
    xlim([-85 0])

    if mod(s, no_cols) == 1
        
        ylabel({['I_{app} = ', num2str(data(s).pop1_I_app, '%.3g')];'iKs Gating Var.'},'FontSize',14) 
    
    end
    
    if ceil(s/no_cols) == no_rows,
    
        xlabel('Voltage','FontSize',14), 
    
    end, 
    
    if ceil(s/no_cols) == 1, 
    
        title(['g_{K_s} = ', num2str(data(s).pop1_gKs, '%.3g')],'FontSize',16)
    
    end
    
end