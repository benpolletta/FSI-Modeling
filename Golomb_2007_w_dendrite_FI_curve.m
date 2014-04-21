function [Vs_all,I,F]=Golomb_2007_w_dendrite_FI_curve(I_range,no_secs,no_steps,theta_m,gD)

% no_secs=1;

dt = .005;
t = (0:dt:(no_secs*1000-dt))/1000;

% no_steps=5;
I = linspace(I_range(1),I_range(2),no_steps);
no_steps = length(I);

F=zeros(no_steps,4);
Vs_all=zeros(no_steps,length(t));

pow_length_test = pmtm(Vs_all(1,:));
pow_length = length(pow_length_test);

pow_all = zeros(no_steps,pow_length);
freqs = zeros(no_steps,pow_length);

parfor s=1:no_steps
    
    [Vs,~,~,~] = Golomb_2007_w_dendritic_gap_jxn(1,I(s),no_secs*1000,theta_m,gD,0);
    
    AP_times = Vs > 0;
    spike_indices = find(diff(AP_times) == 1);
    freq_1 = 1000/(median(diff(spike_indices))*.005);
    if ~isempty(nanmax(diff(spike_indices)))
        freq_3 = 1000/(nanmax(diff(spike_indices))*.005);
    else
        freq_3 = nan;
    end
    freq_4 = length(spike_indices)/no_secs;
    
    [Vs_pow, f] = pmtm(detrend(Vs),[],[],round(length(t)/no_secs));
    f_restricted = f(10<f & f<150);
    [~,max_index] = max(Vs_pow(10<f & f<150));
    freq_2 = f_restricted(max_index);
    
    Vs_all(s,:) = Vs;
    pow_all(s,:) = Vs_pow;
    freqs(s,:) = f;
    F(s,:) = [freq_1 freq_2 freq_3 freq_4];
    
end

freq = mean(freqs);

save(['Golomb_2007_w_dendrite_FI_curve_thm',num2str(theta_m),'_gd',num2str(gD),'.mat'],'Vs_all','pow_all','freq','I','F')

figure;

plot(I',F,'*')
legend({'From Median ISI','From Spectral Peak','From Maximal ISI','From Spikes/Sec.'})
xlabel('Applied Current')
ylabel('Spike Frequency (Hz)')
save_as_pdf(gcf,['Golomb_2007_w_dendrite_FI_curve_thm',num2str(theta_m),'_gd',num2str(gD)])

label_struct = struct('title',sprintf('FI Curve Voltages, theta_m = %g, g_D = %g',theta_m,gD),'xlabel','Time (ms)','ylabel','I_{app}','yticklabel',I);

plot_mat_1axis(Vs_all,t,label_struct,[],[],'k')

save_as_pdf(gcf,['Golomb_2007_w_dendrite_FI_curve_thm',num2str(theta_m),'_gd',num2str(gD),'_Voltages'])

figure;
imagesc(freq(freq<200),I,zscore(pow_all(:,freq<200)))
title(sprintf('FI Curve Spectra, theta_m = %g, g_D = %g',theta_m,gD))
xlabel('Freq. (Hz)')
ylabel('Applied Current')

save_as_pdf(gcf,['Golomb_2007_FI_curve_thm',num2str(theta_m),'_gd',num2str(gD),'Spectra'])