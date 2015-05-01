function Golomb_2007_theta_m_gD_sweep_plotter

T0 = 2000;
dt = 0.005;                       %The time step.
T  = ceil(T0/dt);
t = (1:T)*dt;

[mat_name, mat_dir] = uigetfile('*.mat','Choose a sweep simulation to plot.');
load([mat_dir, mat_name])

[rows, columns] = size(theta_m_mat);

yticklabels = cell(rows,columns);

for rvar=1:rows
    
    for cvar=1:columns
        
        yticklabels{rvar,cvar} = sprintf('I_a = %.2g, th_m = %.2g, g_D = %.2g',I_app_mat(rvar,cvar),theta_m_mat(rvar,cvar),gD_mat(rvar,cvar));
        
    end
    
end

figure()
    
for cvar = 1:columns
    
    ax_handle = subplot(1,columns,cvar);
    
    plot_mat_1axis(VS(:,:,cvar)',t,struct('title','Golomb 2007 Model FSI','xlabel','Time (ms)','ylabel','Parameters','yticklabel',{yticklabels(:,cvar)}),[],ax_handle,'b');
    
end

save_as_pdf(gcf,[mat_dir, mat_name(1:end-4)])