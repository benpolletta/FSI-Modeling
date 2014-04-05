function Golomb_2007_theta_m_gD_sweep(theta_m_vec,gD_vec,I_app_mat)

T0 = 2000;
dt = 0.005;                       %The time step.
T  = ceil(T0/dt);
t = (1:T)*dt;

no_theta_m = length(theta_m_vec);
no_gD = length(gD_vec);

if isempty(I_app_mat)
   
    I_app_mat=4*ones(no_gD,no_theta_m);
    
elseif size(I_app_mat,1)~=no_gD || size(I_app_mat,2)~=no_theta_m
    
    display('I_app_vec must have dimensions length(gD_vec) x length(theta_m_vec).')
    
    return
    
end

VS = nan(T,no_gD,no_theta_m);

for tvar=1:no_theta_m
    
    theta_m=theta_m_vec(tvar);
    
    parfor gDvar=1:no_gD
        
        gD_vec_local=gD_vec;
        
%         index = no_gD*(tvar-1) + gDvar;

        index = no_theta_m*gDvar-tvar+1;
        
        [Vs_temp,~,~,~,~,~] = Golomb_2007(1,I_app_mat(gDvar,tvar),T0,theta_m,gD_vec_local(gDvar));
        
        VS(:,gDvar,tvar) = Vs_temp;
        
    end
    
end

for tvar = 1:no_theta_m

    ax_handle = subplot(no_theta_m,1,no_theta_m-tvar+1);
    
    plot_imf_1axis(VS(:,:,tvar)',t,['\theta_m = ',num2str(theta_m_vec(tvar))],[],ax_handle,'b');
    
end

save(['Golomb_2007_sweep',num2str(10*rand)],'VS','theta_m_vec','gD_vec','I_app_mat')