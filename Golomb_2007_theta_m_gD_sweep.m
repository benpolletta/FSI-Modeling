function Golomb_2007_theta_m_gD_sweep(theta_m_vec,gD_vec,I_app_mat)

% theta_m_vec = -28:4:-24;
no_theta_m = length(theta_m_vec);
% gD_vec = [0.1 0.39 1.8];
no_gD = length(gD_vec);

if isempty(I_app_mat)
   
    I_app_mat=4*ones(no_gD,no_theta_m);
    
elseif size(I_app_mat,1)~=no_gD | size(I_app_mat,2)~=no_theta_m
    
    display('I_app_vec must have dimensions length(gD_vec) x length(theta_m_vec).')
    
    return
    
end

figure()

% index=1;

for tvar=1:no_theta_m
    
    theta_m=theta_m_vec(tvar);
    
    for gDvar=1:no_gD
        
        gD_vec_local=gD_vec;
        
%         index = no_gD*(tvar-1) + gDvar;

        index = no_theta_m*gDvar-tvar+1;
        
        [Vs_temp,~,~,~,~,t] = Golomb_2007(1,I_app_mat(gDvar,tvar),2000,theta_m,gD_vec_local(gDvar));
        
        subplot(no_gD,no_theta_m,index)
        
        plot(t,Vs_temp)
        
%         index=index+1;
        
    end
    
end