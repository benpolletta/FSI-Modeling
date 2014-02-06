function [Vs_all,I,F]=FI_curve(I_range,no_secs,no_steps,g_autapse)

% no_secs=1;

dt = .005;
t = (0:dt:(no_secs*1000-dt))/1000;
f = (1000/dt)*(1:(1000/dt))/(1000/dt);

% no_steps=5;
I_step=(I_range(2)-I_range(1))/no_steps;

I=zeros(1,no_steps+1);
F=zeros(1,no_steps+1);
Vs_all=zeros(no_steps+1,length(t));

parfor i=1:(no_steps+1)
   
    I_range_temp=I_range;
    f_temp=f;
    I_temp=I_range_temp(1)+I_step*(i-1);
    
    [Vs,~,~,~] = ing_w_dendritic_gap_jxn(1,I_temp,no_secs*1000,[],g_autapse,0);
    
    Vs_pow = fft(detrend(Vs((end+1-(1000/dt)):end))).^2;
    [~,max_index] = max(Vs_pow(f_temp<150));
    
    Vs_all(i,:) = Vs;
    I(i) = I_temp;
    F(i) = f_temp(max_index);
    
end

figure()

plot(I,F,'*')
xlabel('Applied Current')
ylabel('Spike Frequency (Hz)')
saveas(gcf,['FI_curve_',num2str(I_range(1)),'to',num2str(I_range(2)),'_g_aut_',num2str(g_autapse),'.fig'])

figure()

Vs_height=(max(max(Vs_all))-min(min(Vs_all)));

for i=1:(no_steps+1)
    
    plot(t,Vs_all(i,:)+Vs_height*(i-1))
    hold on
    
end

set(gca,'YTick',0:Vs_height:(Vs_height*(no_steps+1)),'YTickLabel',I)
saveas(gcf,['Voltages_',num2str(I_range(1)),'to',num2str(I_range(2)),'_g_aut_',num2str(g_autapse),'.fig'])