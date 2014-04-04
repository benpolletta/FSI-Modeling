function Golomb_2007_theta_m_gD_sweep(theta_m_vec,gD_vec,I_app_vec)

% theta_m_vec = -28:4:-24;
no_theta_m = length(theta_m_vec);
% gD_vec = [0.1 0.39 1.8];
no_gD = length(gD_vec);

if isempty(I_app_vec)
   
    I_app_vec=4*ones(no_gD,1);
    
elseif length(I_app_vec)~=no_gD
    
    display('I_app_vec must be the same length as gD_vec.')
    
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
        
        [Vs_temp,~,~,~,~,t] = Golomb_2007(1,I_app_vec(gDvar),2000,theta_m,gD_vec_local(gDvar));
        
        subplot(no_gD,no_theta_m,index)
        
        plot(t,Vs_temp)
        
%         index=index+1;
        
    end
    
end