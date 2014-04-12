function Golomb_2007_theta_m_gD_sweep(theta_m_mat,gD_mat,I_app_mat)

T0 = 2000;
dt = 0.005;                       %The time step.
T  = ceil(T0/dt);
t = (1:T)*dt;

[rtm, ctm] = size(theta_m_mat);

[rgd, cgd] = size(gD_mat);

[ria, cia] = size(I_app_mat);

if ~(rtm==rgd && rgd==ria) || ~(ctm==cgd && cgd==cia)
    
    display('All three input matrices must have the same dimensions.')
    
    return
    
else
    
    rows = rtm; columns = ctm;
    
end

VS = nan(T,rows,columns);
yticklabels = cell(rows,columns);

for rvar=1:rows
    
    parfor cvar=1:columns
        
        [Vs_temp,~,~,~,~,~] = Golomb_2007(1,I_app_mat(rvar,cvar),T0,theta_m_mat(rvar,cvar),gD_mat(rvar,cvar));
        
        VS(:,rvar,cvar) = Vs_temp;
        
        yticklabels{rvar,cvar} = sprintf('I_a = %.2g, th_m = %.2g, g_D = %.2g',I_app_mat(rvar,cvar),theta_m_mat(rvar,cvar),gD_mat(rvar,cvar));
        
    end
    
end

date_string = datestr(now,'dd-mm-yy_HH-MM-SS');

try
    
    for cvar = 1:columns
        
        ax_handle = subplot(1,columns,cvar);
        
        plot_mat_1axis(VS(:,:,cvar)',t,struct('title','Golomb 2007 Model FSI','xlabel','Time (ms)','ylabel','Parameters','yticklabel',{yticklabels(:,cvar)}),[],ax_handle,'b');
        
    end
    
    save_as_pdf(gcf,['Golomb_2007_sweep_theta_m_gD_',date_string])

catch error
    
    display(['Could not plot figure: ',error.message])
    save(['Golomb_2007_sweep_theta_m_gD_',date_string,'_error.mat'],'error')
    
end

save(['Golomb_2007_sweep_theta_m',date_string,'.mat'],'VS','theta_m_mat','gD_mat','I_app_mat','t')