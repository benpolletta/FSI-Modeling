function [Vs_all,SI,SI_times,gc_steps,F]=G_gj_sweep_time_series(gc_range,no_secs,no_steps,no_reps)

% no_secs=1;

dt = .005;
t = (0:dt:(no_secs*1000-dt))/1000;
f = (1000/dt)*(1:(1000/dt))/(1000/dt);

% no_steps=5;
gc_step=(gc_range(2)-gc_range(1))/no_steps;
gc_steps=gc_range(1):gc_step:gc_range(2);

SI=cell(no_reps,no_steps+1);
SI_times=cell(no_reps,no_steps+1);
Vs_all=zeros(no_steps+1,length(t),2,no_reps);
F=zeros(no_reps,no_steps+1);

for r=1:no_reps
    
    for i=1:(no_steps+1)
        
        gc_steps_temp=gc_steps;
        f_temp=f;
        
        [Vs,~,~,~] = ing_w_dendritic_gap_jxn(2,100,no_secs*1000,[],0,[0 1; 1 0]*gc_steps_temp(i));
       
        [SI{r,i},SI_times{r,i}]=sync_time_series(Vs,dt);

        Vs_pow = fft(detrend(Vs(1,(end+1-(1000/dt)):end))).^2;
        [~,max_index] = max(Vs_pow(f_temp<300));        
        F(r,i) = f_temp(max_index);
        
        Vs_all(i,:,:,r) = reshape(Vs',1,size(Vs,2),2);
        
    end
    
end

gc_colors=[repmat(linspace(1,0,no_steps+1),2,1); linspace(0,1,no_steps+1)]';

gc_legend=cell(1,no_steps+1);

for s=1:(no_steps+1)
    
    gc_legend{s}=['g_{gap} = ',num2str(gc_steps(s))];
    
end

h=zeros(no_steps+1,1);

figure()
for r=1:no_reps
    for i=1:(no_steps+1)
            h(i)=plot(SI_times{r,i},SI{r,i},'o','Color',gc_colors(i,:));
            hold on
    end
end
xlabel('Time (ms)')
ylabel('(Inter-Cell ISI)/(Intra-Cell ISI)')
legend(h,gc_legend)
saveas(gcf,['gc_time_series_',num2str(gc_range(1)),'to',num2str(gc_range(2)),'.fig'])

% figure()
% 
% Vs_height=(max(max(Vs_all))-min(min(Vs_all)));
% 
% for i=1:(no_steps+1)
%     
%     plot(t,Vs_all(i,:)+Vs_height*(i-1))
%     hold on
%     
% end
% 
% set(gca,'YTick',0:Vs_height:(Vs_height*(no_steps+1)),'YTickLabel',gc)
% saveas(gcf,['V_sum_',num2str(gc_range(1)),'to',num2str(gc_range(2)),'.fig'])