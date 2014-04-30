function [Vs_all,ISIm,gc,F]=gap_jxn_conductance_sweep(gc_range,no_secs,no_steps,no_reps)

% no_secs=1;

dt = .005;
t = (0:dt:(no_secs*1000-dt))/1000;
f = (1000/dt)*(1:(1000/dt))/(1000/dt);

% no_steps=5;
gc_step=(gc_range(2)-gc_range(1))/no_steps;

gc=zeros(1,no_steps+1);
ISIm=zeros(no_reps,no_steps+1);
Vs_all=zeros(no_steps+1,length(t),2,no_reps);
F=zeros(no_reps,no_steps+1);

for r=1:no_reps
    
    parfor i=1:(no_steps+1)
        
        gc_range_temp=gc_range;
        f_temp=f;
        gc_temp=gc_range_temp(1)+gc_step*(i-1);
        
        [Vs,~,~,~] = ing_w_dendritic_gap_jxn(2,100,no_secs*1000,[],0,[0 1; 1 0]*gc_temp);
       
        ISIm(r,i)=mean_ISI(Vs(:,(end+1-(1000/dt)):end),dt);

        Vs_pow = fft(detrend(Vs(1,(end+1-(1000/dt)):end))).^2;
        [~,max_index] = max(Vs_pow(f_temp<300));        
        F(r,i) = f_temp(max_index);
        
        Vs_all(i,:,:,r) = reshape(Vs',1,size(Vs,2),2);
        gc(r,i) = gc_temp;
        
    end
    
end

figure()
scatter(reshape(gc,no_reps*(no_steps+1),1),reshape(ISIm,no_reps*(no_steps+1),1))
xlabel('Gap Junction Conductance')
ylabel('(Mean Inter-Cell ISI)/(Mean Intra-Cell ISI)')
saveas(gcf,['gc_',num2str(gc_range(1)),'to',num2str(gc_range(2)),'.fig'])

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