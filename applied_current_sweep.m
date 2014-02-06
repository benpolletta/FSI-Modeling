function [Vs_all,ISIm,I,F_cell,F_pop]=applied_current_sweep(I_range,no_secs,no_steps,no_reps)

% no_secs=1;

gsd=0.25;
ggap=.1;
gi=30;

dt = .005;
t = (0:dt:(no_secs*1000-dt))/1000;
f = (1000/dt)*(1:(1000/dt))/(1000/dt);

% no_steps=5;
I_step=(I_range(2)-I_range(1))/no_steps;
I_steps=I_range(1):I_step:I_range(2);

I=zeros(1,no_steps+1);
ISIm=zeros(no_reps,no_steps+1);
Vs_all=zeros(no_steps+1,length(t),2,no_reps);
F_cell=zeros(no_reps,no_steps+1);
F_pop=zeros(no_reps,no_steps+1);

for r=1:no_reps
    
    parfor i=1:(no_steps+1)
        
        I_steps_temp=I_steps;
        f_temp=f;
        
        [Vs,~,~,~] = ing_w_dendritic_gap_jxn(2,I_steps_temp(i),no_secs*1000,gsd,[0 1;1 0]*gi,[0 1; 1 0]*ggap);
       
        ISIm(r,i)=mean_ISI(Vs(:,(end+1-(1000/dt)):end),dt);

        Vs_pow = fft(detrend(Vs(1,(end+1-(1000/dt)):end))).^2;
        [~,max_index] = max(Vs_pow(f_temp<300));        
        F_cell(r,i) = f_temp(max_index);
        
        Vs_pop_pow = fft(detrend(sum(Vs(:,(end+1-(1000/dt)):end)))).^2;
        [~,max_index] = max(Vs_pop_pow(f_temp<300));        
        F_pop(r,i) = f_temp(max_index);
        
        Vs_all(i,:,:,r) = reshape(Vs',1,size(Vs,2),2);
        I(r,i) = I_steps_temp(i);
        
    end
    
end

figure()
scatter(reshape(I,no_reps*(no_steps+1),1),reshape(ISIm,no_reps*(no_steps+1),1))
xlabel('Applied Current')
ylabel('(Mean Inter-Cell ISI)/(Mean Intra-Cell ISI)')
saveas(gcf,['I_',num2str(I_range(1)),'to',num2str(I_range(2)),'_gsd_',num2str(0.25),'_ggap_',num2str(.1),'_gi_',num2str(gi),'_',num2str(no_secs),'s_SI.fig'])

[rows,cols]=subplot_size(no_reps);

figure()

scatter(reshape(I,no_reps*(no_steps+1),1),reshape(F_cell,no_reps*(no_steps+1),1),'o')
xlabel('Applied Current')
ylabel('Cell Frequency')
saveas(gcf,['I_',num2str(I_range(1)),'to',num2str(I_range(2)),'_gsd_',num2str(0.25),'_ggap_',num2str(.1),'_gi_',num2str(gi),'_',num2str(no_secs),'s_Fcell.fig'])

figure()

scatter(reshape(I,no_reps*(no_steps+1),1),reshape(F_pop,no_reps*(no_steps+1),1),'o')
xlabel('Applied Current' )
ylabel('Population Frequency')
saveas(gcf,['I_',num2str(I_range(1)),'to',num2str(I_range(2)),'_gsd_',num2str(0.25),'_ggap_',num2str(.1),'_gi_',num2str(gi),'_',num2str(no_secs),'s_Fpop.fig'])

for i=1:(no_steps+1)
    
    figure()
    
    for j=1:no_reps
        
        Vs=Vs_all(i,:,:,j);
        Vs=reshape(Vs,size(Vs,2),1,2); 
        Vs=reshape(Vs,size(Vs,1),2); 
        
        subplot(rows,cols,j) 
        
        plot(Vs((end+1-10^4):end,:))
        hold on
    
    end
    
    saveas(gcf,['I_',num2str(I_steps(i)),'_gsd_',num2str(0.25),'_ggap_',num2str(.1),'_gi_',num2str(gi),'_',num2str(no_secs),'s_traces.fig'])
   
end

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