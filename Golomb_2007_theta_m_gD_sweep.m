function Golomb_2007_theta_m_gD_sweep(theta_m_mat,gD_mat,I_app_mat)

T0 = 4000;
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

parfor cvar=1:columns
    
    I_app = I_app_mat(:,cvar);
    theta_m = theta_m_mat(:,cvar);
    gD = gD_mat(:,cvar);
    
    [Vs_temp,~,~,~,~,~] = Golomb_2007_RK45_w_dendrite(rows,I_app,T0,theta_m,gD,1,0,0);
    
    VS(:,:,cvar) = Vs_temp';
    
    for rvar = 1:rows
    
        yticklabels{rvar,cvar} = sprintf('I_a = %g, th_m = %g, g_D = %g',I_app(rvar),theta_m(rvar),gD(rvar));
    
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

save(['Golomb_2007_sweep_theta_m_gD_',date_string,'.mat'],'VS','theta_m_mat','gD_mat','I_app_mat','t')