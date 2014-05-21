function lewis_rinzel_2003_gating_variables

V=[-100:.1:100];

vars=[alphaM(V); betaM(V); alphaH(V); betaH(V); alphaN31_32(V); betaN31_32(V); alphaN13(V); betaN13(V)];

letters={'\alpha','\beta'};
labels={'m','h','nKv3.1','nKv1.3'};
param_labels={'x_{\inf}','\tau'};

vars_mat=reshape(vars,2,4*length(V));
tau=1./sum(vars_mat);
m_inf=vars_mat(1,:).*tau;
params_mat = [reshape(m_inf,4,length(V))'; reshape(tau,4,length(V))'];
params=reshape(params_mat,length(V),8)';

h = figure; g = figure;

for i=1:8
    
    figure(h)
    
    subplot(4,2,i)
    plot(V,vars(i,:))
    ylabel([letters{2-mod(i,2)},'(',labels{ceil(i/2)},')'])
    
    figure(g)
    
    subplot(4,2,i)
    plot(V,params(i,:))
    ylabel([param_labels{2-mod(i,2)},'(',labels{ceil(i/2)},')'])
    
end

save_as_pdf(h,'lewis_rinzel_alpha_beta')

save_as_pdf(g,'lewis_rinzel_x_inf_tau')

end

%Below, define the auxiliary functions alpha & beta for each gating variable.

function aM = alphaM(V)
aM = (3020 - 40*V)./(exp((75 + V)/(13.5)) - 1);
end

function bM = betaM(V)
bM = 1.2262./exp(V/42.248);
end

function aH = alphaH(V)
aH = 0.0035./exp(V/24.186);
end

function bH = betaH(V)
bH = -(0.8712 + 0.017*V)./(exp((51.25 + V)/(-5.2)) - 1);
end

function aN = alphaN31_32(V)
aN = (95 - V)./(exp((V - 95)/(-11.8)) - 1);
end

function bN = betaN31_32(V)
bN = 0.025./exp(V/22.222);
end

function aN = alphaN13(V)
aN = -(0.616 + 0.014*V)./(exp((44 + V)/(-2.3)) - 1);
end

function bN = betaN13(V)
bN = 0.0043./exp((44 + V)/34);
end