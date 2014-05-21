function kotaleski_gating_variables

V= -100:.1:100;

vars=[alpha_m(V); beta_m(V); alpha_h(V); beta_h(V); alpha_n3132(V); beta_n3132(V); alpha_n13(V); beta_n13(V)];

letters={'\alpha','\beta'};
labels={'m','h','nKv3.1','nKv1.3'};
param_labels={'x_{\inf}','\tau'};

vars_mat=reshape(vars,2,4*length(V));
tau=1./sum(vars_mat);
x_inf=vars_mat(1,:).*tau;
params_mat = [reshape(x_inf,4,length(V))'; reshape(tau,4,length(V))'];
params=reshape(params_mat,length(V),8)';

h = figure; g = figure;

for i=1:8
    
    figure(h)
    
    subplot(6,2,i)
    plot(V,vars(i,:))
    ylabel([letters{2-mod(i,2)},'(',labels{ceil(i/2)},')'])
    
    figure(g)
    
    subplot(6,2,i)
    plot(V,params(i,:))
    ylabel([param_labels{2-mod(i,2)},'(',labels{ceil(i/2)},')'])
    
end

clear vars params m_inf

vars = [m_inf(V); tau_m(V); h_inf(V); tau_h(V)];

letters={'x_{\inf}','\tau'};
labels={'mKA','hKA'};
param_labels={'\alpha','\beta'};

vars_mat = reshape(vars,2,length(labels)*length(V));
alpha = vars_mat(1,:)./vars_mat(2,:);
beta = 1./vars_mat(2,:) - alpha;
params_mat = [reshape(alpha,length(labels),length(V))'; reshape(beta,length(labels),length(V))'];
params = reshape(params_mat,length(V),2*length(labels))';

for i=1:4
    
    figure(g)
    
    subplot(6,2,8+i)
    plot(V,vars(i,:))
    ylabel([letters{2-mod(i,2)},'(',labels{ceil(i/2)},')'])
    
    figure(h)
    
    subplot(6,2,8+i)
    plot(V,params(i,:))
    ylabel([param_labels{2-mod(i,2)},'(',labels{ceil(i/2)},')'])
    
end

save_as_pdf(h,'kotaleski_alpha_beta')

save_as_pdf(g,'kotaleski_x_inf_tau')

end

%Below, define the auxiliary functions alpha & beta for each gating variable.

function a_M = alpha_m(V)
num = 3020 - 40*V;
denom = exp((V+75)/13.5) - 1;
a_M = num./denom;
end

function b_M = beta_m(V)
num = 1.2262;
denom = exp(V/42.248);
b_M = num./denom;
end

function a_H = alpha_h(V)
num = .0035;
denom = exp(V/24.186);
a_H = num./denom;
end

function b_H = beta_h(V)
num = -(0.8712 + 0.017*V);
denom = exp((V + 51.25)/(-5.2)) - 1;
b_H = num./denom;
end

function a_N = alpha_n13(V)
num = -(0.616 + 0.014*V);
denom = exp((V+44)/(-2.3)) - 1;
a_N = num./denom;
end

function b_N = beta_n13(V)
num = 0.0043;
denom = exp((V+44)/(34));
b_N = num./denom;
end

function a_N = alpha_n3132(V)
num = 95 - V;
denom = exp((V - 95)/(-11.8)) - 1;
a_N = num./denom;
end

function b_N = beta_n3132(V)
num = 0.025;
denom = exp(V/22.222);
b_N = num./denom;
end

function m_i = m_inf(V)
num = 0.001;
denom = exp((V + 45)/(-13)) + 1;
m_i = num./denom;
end

function t_M = tau_m(V)
t_M = 0.000001*(exp((V + 70)/(-13)) + 1);
end

function h_i = h_inf(V)
num = 0.001;
denom = exp((V + 77)/8) + 1;
h_i = num./denom;
end

function t_h = tau_h(V)
t_h = 14*ones(size(V));
end