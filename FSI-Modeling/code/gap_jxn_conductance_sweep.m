function [Vs_all,ISIm,ggap,F_cell,F_pop]=gap_jxn_conductance_sweep(ggap_range,no_secs,no_steps,no_reps)

% no_secs=1;

I_app=100;
gi=10;
gsd=[];

dt = .005;
t = (0:dt:(no_secs*1000-dt))/1000;
f = (1000/dt)*(1:(1000/dt))/(1000/dt);

% no_steps=5;
ggap_step=(ggap_range(2)-ggap_range(1))/no_steps;
ggap_steps=ggap_range(1):ggap_step:ggap_range(2);

ggap=zeros(no_reps,no_steps+1);
ISIm=zeros(no_reps,no_steps+1);
Vs_all=zeros(no_steps+1,length(t),2,no_reps);
F_cell=zeros(no_reps,no_steps+1);
F_pop=zeros(no_reps,no_steps+1);

for r=1:no_reps
    
    parfor i=1:(no_steps+1)
        
        ggap_steps_temp=ggap_steps;
        f_temp=f;
        
        [Vs,~,~,~] = ing_w_dendritic_gap_jxn(2,I_app,no_secs*1000,gsd,[0 1; 1 0]*gi,[0 1; 1 0]*ggap_steps_temp(i));
       
        ISIm(r,i)=mean_ISI(Vs(:,(end+1-(1000/dt)):end),dt);

        cell_pow = fft(detrend(Vs(1,(end+1-(1000/dt)):end))).^2;
        [~,max_index] = max(cell_pow(f_temp<300));        
        F_cell(r,i) = f_temp(max_index);
        
        pop_pow = fft(detrend(sum(Vs(:,(end+1-(1000/dt)):end)))).^2;
        [~,max_index] = max(pop_pow(f_temp<300));        
        F_pop(r,i) = f_temp(max_index);
        
        Vs_all(i,:,:,r) = reshape(Vs',1,size(Vs,2),2);
        ggap(r,i) = ggap_steps_temp(i);
        
    end
    
end

figure()
scatter(reshape(ggap,no_reps*(no_steps+1),1),reshape(ISIm,no_reps*(no_steps+1),1))
xlabel('Gap Junction Conductance')
ylabel('(Mean Inter-Cell ISI)/(Mean Intra-Cell ISI)')
saveas(gcf,['ggap_',num2str(ggap_range(1)),'to',num2str(ggap_range(2)),'_I_',num2str(I_app),'_gsd_',num2str(gsd),'_gi_',num2str(gi),'_',num2str(no_secs),'s_SI.fig'])

figure()

scatter(reshape(ggap,no_reps*(no_steps+1),1),reshape(F_cell,no_reps*(no_steps+1),1),'o')
xlabel('Applied Current')
ylabel('Cell Frequency')
saveas(gcf,['ggap_',num2str(ggap_range(1)),'to',num2str(ggap_range(2)),'_I_',num2str(I_app),'_gsd_',num2str(gsd),'_gi_',num2str(gi),'_',num2str(no_secs),'s_Fcell.fig'])

figure()

scatter(reshape(ggap,no_reps*(no_steps+1),1),reshape(F_pop,no_reps*(no_steps+1),1),'o')
xlabel('Applied Current' )
ylabel('Population Frequency')
saveas(gcf,['ggap_',num2str(ggap_range(1)),'to',num2str(ggap_range(2)),'_I_',num2str(I_app),'_gsd_',num2str(gsd),'_gi_',num2str(gi),'_',num2str(no_secs),'s_Fpop.fig'])

[rows,cols]=subplot_size(no_reps);

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
    
    saveas(gcf,['ggap_',num2str(ggap_range(1)),'to',num2str(ggap_range(2)),'_I_',num2str(I_app),'_gsd_',num2str(gsd),'_gi_',num2str(gi),'_',num2str(no_secs),'s_traces.fig'])
   
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
% set(gca,'YTick',0:Vs_height:(Vs_height*(no_steps+1)),'YTickLabel',ggap)
% saveas(gcf,['V_sum_',num2str(ggap_range(1)),'to',num2str(ggap_range(2)),'.fig'])