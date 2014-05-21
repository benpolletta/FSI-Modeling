function buszaki_wang_gating_variables

V=-100:.1:100;

vars=[alphaM(V); betaM(V); alphaH(V); betaH(V); alphaN(V); betaN(V)];

letters={'\alpha','\beta'};
labels={'m','h','n'};
param_labels={'x_{\inf}','\tau'};

vars_mat=reshape(vars,2,length(labels)*length(V));
tau=1./sum(vars_mat);
m_inf=vars_mat(1,:).*tau;
params_mat = [reshape(m_inf,length(labels),length(V))'; reshape(tau,length(labels),length(V))'];
params=reshape(params_mat,length(V),2*length(labels))';

h = figure; g = figure;

for i=1:2*length(labels)
    
    figure(h)
    
    subplot(length(labels),2,i)
    plot(V,vars(i,:))
    ylabel([letters{2-mod(i,2)},'(',labels{ceil(i/2)},')'])
    
    figure(g)
    
    subplot(length(labels),2,i)
    plot(V,params(i,:))
    ylabel([param_labels{2-mod(i,2)},'(',labels{ceil(i/2)},')'])
    
end

save_as_pdf(h,'buszaki_wang_alpha_beta')

save_as_pdf(g,'buszaki_wang_x_inf_tau')

end

%Below, define the auxiliary functions alpha & beta for each gating variable.

function aM = alphaM(V)
aM = (2.5-0.1*(V+65)) ./ (exp(2.5-0.1*(V+65)) -1);
end

function bM = betaM(V)
bM = 4*exp(-(V+65)/18);
end

function aH = alphaH(V)
aH = 0.07*exp(-(V+65)/20);
end

function bH = betaH(V)
bH = 1./(exp(3.0-0.1*(V+65))+1);
end

function aN = alphaN(V)
aN = (0.1-0.01*(V+65)) ./ (exp(1-0.1*(V+65)) -1);
end

function bN = betaN(V)
bN = 0.125*exp(-(V+65)/80);
end