function [Vs_all,I,F]=gap_jxn_current_sweep(I_range,no_secs,no_steps,no_reps)

% no_secs=1;

dt = .005;
t = (0:dt:(no_secs*1000-dt))/1000;
f = (1000/dt)*(1:(1000/dt))/(1000/dt);

% no_steps=5;
I_step=(I_range(2)-I_range(1))/no_steps;

I=zeros(1,no_steps+1);
F=zeros(no_reps,no_steps+1);
Vs_all=zeros(no_steps+1,length(t));

for r=1:no_reps
    
    parfor i=1:(no_steps+1)
        
        I_range_temp=I_range;
        f_temp=f;
        I_temp=I_range_temp(1)+I_step*(i-1);
        
        [Vs,~,~,~] = ing_w_dendritic_gap_jxn(2,I_temp,no_secs*1000,0,[0 1; 1 0]*.05);
        
        Vs_sum=sum(Vs);
%         Vs_pow = fft(detrend(Vs_sum((end+1-(1000/dt)):end))).^2;
%         [~,max_index] = max(Vs_pow(f_temp<300));
        S(r,i) = sum(abs(diff(Vs)))/sum(sum(abs(Vs)));

%         F(r,i) = f_temp(max_index);
        
        Vs_all(i,:) = Vs_sum;
        I(i) = I_temp;
        
    end
    
end

figure()

boxplot(F,I')
xlabel('Applied Current')
ylabel('Population Frequency (Hz)')
saveas(gcf,['gap_jxn_',num2str(I_range(1)),'to',num2str(I_range(2)),'.fig'])

figure()

Vs_height=(max(max(Vs_all))-min(min(Vs_all)));

for i=1:(no_steps+1)
    
    plot(t,Vs_all(i,:)+Vs_height*(i-1))
    hold on
    
end

set(gca,'YTick',0:Vs_height:(Vs_height*(no_steps+1)),'YTickLabel',I)
saveas(gcf,['V_sum_',num2str(I_range(1)),'to',num2str(I_range(2)),'.fig'])