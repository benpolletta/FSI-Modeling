function Golomb_2007_gating_variables


theta_m = -28; sigma_m = 11.5;      %Sodium parameters.
theta_h = -58.3; sigma_h = -6.7;    

theta_n = -12.4; sigma_n = 6.8;     %Delayed rectifier parameters.

theta_a = -50; sigma_a = 20;        %Potassium (slow inactivation) parameters.
theta_b = -70; sigma_b = -6;
tau_a = 2; tau_b = 150;

V=-100:.1:100;

vars=[steady(V,theta_m,sigma_m); ones(size(V)); steady(V,theta_h,sigma_h); tau_h(V); steady(V,theta_n,sigma_n); ...
    tau_n(V); steady(V,theta_a,sigma_a); tau_a*ones(size(V)); steady(V,theta_b,sigma_b); tau_b*ones(size(V))];

letters={'x_{\inf}','\tau'};
labels={'m','h','nKv3.1','aKslow','bKslow'};
param_labels={'\alpha','\beta'};

vars_mat = reshape(vars,2,5*length(V));
alpha = vars_mat(1,:)./vars_mat(2,:);
beta = 1./vars_mat(2,:) - alpha;
params_mat = [reshape(alpha,5,length(V))'; reshape(beta,5,length(V))'];
params = reshape(params_mat,length(V),10)';

h = figure;
g = figure;

for i=1:10
    
    figure(h)
    
    subplot(5,2,i)
    plot(V,vars(i,:))
    ylabel([letters{2-mod(i,2)},'(',labels{ceil(i/2)},')'])
    
    figure(g)
    
    subplot(5,2,i)
    plot(V,params(i,:))
    ylabel([param_labels{2-mod(i,2)},'(',labels{ceil(i/2)},')'])
    
end

save_as_pdf(h,'golomb_alpha_beta')

save_as_pdf(g,'golomb_x_inf_tau')

end

%Below, define the auxiliary functions alpha & beta for each gating variable.

function s = steady(V,theta,sigma)
s = 1./(1+exp((theta-V)/sigma));
end

function th = tau_h(V)
theta_th = -60; sigma_th = -12;
th = .5 + 14./(1+exp((theta_th-V)/sigma_th));
end

function tn = tau_n(V)
tn_1 = .087 + (11.4)./(1+exp((V+14.6)/8.6));
tn_2 = .087 + (11.4)./(1+exp((1.3-V)/18.7));
tn = tn_1.*tn_2;
end