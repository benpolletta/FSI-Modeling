function [Vs_all,ISIm,I,F]=applied_current_sweep(I_range,no_secs,no_steps,no_reps)

% no_secs=1;

dt = .005;
t = (0:dt:(no_secs*1000-dt))/1000;
f = (1000/dt)*(1:(1000/dt))/(1000/dt);

% no_steps=5;
I_step=(I_range(2)-I_range(1))/no_steps;

I=zeros(1,no_steps+1);
F=zeros(no_reps,no_steps+1);
ISIm=zeros(no_reps,no_steps+1);
Vs_all=zeros(no_steps+1,length(t),2,no_reps);

for r=1:no_reps
    
    parfor i=1:(no_steps+1)
        
        I_range_temp=I_range;
        f_temp=f;
        I_temp=I_range_temp(1)+I_step*(i-1);
        
        [Vs,~,~,~] = ing_w_dendritic_gap_jxn(2,I_temp,no_secs*1000,0,0);
        
        ISIm(r,i)=mean_ISI(Vs(:,(end+1-(1000/dt)):end),dt);
        
        Vs_pow = fft(detrend(Vs(1,(end+1-(1000/dt)):end))).^2;
        [~,max_index] = max(Vs_pow(f_temp<150));
        
        Vs_all(i,:,:,r) = reshape(Vs',1,size(Vs,2),2);
        I(i) = I_temp;
        F(r,i) = f_temp(max_index);
        
    end
    
end

figure()

plot(I,F,'*')
xlabel('Applied Current')
ylabel('Spike Frequency (Hz)')
saveas(gcf,['FI_curve_',num2str(I_range(1)),'to',num2str(I_range(2)),'V.fig'])

figure()

boxplot(ISIm,I)
xlabel('Applied Current')
ylabel('Mean ISI (ms)')
saveas(gcf,['Mean_ISI_',num2str(I_range(1)),'to',num2str(I_range(2)),'V.fig'])

