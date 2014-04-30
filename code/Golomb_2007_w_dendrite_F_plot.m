function [gD,I,F]=Golomb_2007_w_dendrite_F_plot(I_range,no_I_steps,gD_range,no_gD_steps,no_secs,theta_m)

I = linspace(I_range(1),I_range(2),no_I_steps);
no_I_steps = length(I);

gD = linspace(gD_range(1),gD_range(2),no_gD_steps);
no_gD_steps = length(gD);

F=zeros(no_gD_steps,no_I_steps,2);

parfor Ivar = 1:no_I_steps
    
    [Vs,~,~,~] = Golomb_2007_w_dendritic_gap_jxn(no_gD_steps,I(Ivar),no_secs*1000,theta_m,gD',0,0,0);
    
    AP_times = Vs > 0;
    
    freqs = nan(no_gD_steps,1,2);
    
    for gDvar = 1:no_gD_steps
    
        spike_indices = find(diff(AP_times(gDvar,:)') == 1);
        freqs(gDvar,1,1) = 1000./(nanmedian(diff(spike_indices))*.005);
        if ~isempty(nanmax(diff(spike_indices)))
            freqs(gDvar,1,2) = 1000./(nanmax(diff(spike_indices))*.005);
        end
        
    end
    
    F(:,Ivar,:) = freqs;
    
end

save(['Golomb_2007_w_dendrite_F_plot_',num2str(theta_m),'.mat'],'I','gD','F')

figure;
colorplot(F(:,:,1),I,gD,10)
colorbar
title(sprintf('Frequency Plot, theta_m = %g.',theta_m))
ylabel('Slow Potassium Conductance (nS)')
xlabel('Applied Current (mS/cm^2)')

save_as_pdf(gcf,['Golomb_2007_w_dendrite_F_plot_',num2str(theta_m)])

figure;
colorplot(-diff(F,[],3),I,gD,10)
colorbar
title(sprintf('Maximal vs. Minimal ISI Freq., theta_m = %g.',theta_m))
ylabel('Slow Potassium Conductance (nS)')
xlabel('Applied Current (mS/cm^2)')

save_as_pdf(gcf,['Golomb_2007_w_dendrite_burst_plot_',num2str(theta_m)])