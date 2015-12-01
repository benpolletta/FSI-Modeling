% Model: FSI_network
cd /Users/benjaminpittman-polletta/Documents/Science/Research_Projects/Modeling/FSI-Modeling/FSI_network_dnsim;
spec=[];
spec.nodes(1).label = 'soma';
spec.nodes(1).multiplicity = 100;
spec.nodes(1).dynamics = {'v''=current/c'};
spec.nodes(1).mechanisms = {'soma_K','soma_Na','soma_leak'};
spec.nodes(1).parameters = {'c',1,'v_IC',-65};
spec.nodes(2).label = 'dendrite';
spec.nodes(2).multiplicity = 100;
spec.nodes(2).dynamics = {'v''=current/c'};
spec.nodes(2).mechanisms = {'dendrite_K','dendrite_Na','dendrite_input','dendrite_leak','dendrite_iMultiPoissonExp'};
spec.nodes(2).parameters = {'c',1,'v_IC',-65};
spec.connections(1,1).label = 'soma-soma';
spec.connections(1,1).mechanisms = {'soma_soma_iSYN'};
spec.connections(1,1).parameters = [];
spec.connections(1,2).label = 'soma-dendrite';
spec.connections(1,2).mechanisms = {'soma_dendrite_iCOM'};
spec.connections(1,2).parameters = [];
spec.connections(2,1).label = 'dendrite-soma';
spec.connections(2,1).mechanisms = {'dendrite_soma_iCOM'};
spec.connections(2,1).parameters = [];
spec.connections(2,2).label = 'dendrite-dendrite';
spec.connections(2,2).mechanisms = {'dendrite_dendrite_iGAP'};
spec.connections(2,2).parameters = [];
%dnsim(spec); % open model in DNSim GUI

% DNSim simulation and plots:
data = runsim(spec,'timelimits',[0 10000],'dt',.02,'SOLVER','euler','savedata_flag',0,'timesurfer_flag',0); % simulate DNSim models
% plotv(data,spec,'varlabel','v'); % quickly plot select variables
v_soma = data.soma_v;
time = data.time;

subplot(3, 1, 1)
imagesc(time(5001:end), 1:spec.nodes(1).multiplicity, v_soma(5001:end, :)')
colorbar
ylabel('Cell #')
xlabel('Time (ms)')

v_mean = mean(v_soma(5001:end, :)');

subplot(3, 1, 2)
plot(time(5001:end), mean(v_soma(5001:end, :)'))   
axis tight
xlabel('Time (ms)')
ylabel('Mean Voltage (mV)')

[v_mean_hat, f] = pmtm(v_mean, [], [], 1000/.02);

subplot(3, 1, 3)
loglog(f, v_mean_hat)
axis tight
xlabel('Frequency (Hz)')
ylabel('Spectral Power')

save_as_eps(gcf, 'FSI_network_dnsim_output')
%visualizer(data); % ugly interactive tool hijacked to visualize sim_data

% Sweep over parameter values:
% model=buildmodel(spec); % parse DNSim spec structure
% simstudy(model,{'soma'},{'N'},{'[1 2]'},'timelimits',[0 100],'dt',.02,'SOLVER','euler'); % N = # of cells


% Manual simulation and plots:
%{
%-----------------------------------------------------------
% Auxiliary variables:
	soma_soma_iSYN_width = inf;
	soma_soma_iSYN_Nmax  = max((100),(100));
	soma_soma_iSYN_srcpos = linspace(1,soma_soma_iSYN_Nmax,(100))'*ones(1,(100));
	soma_soma_iSYN_dstpos = (linspace(1,soma_soma_iSYN_Nmax,(100))'*ones(1,(100)))';
	soma_soma_iSYN_netcon = rand((100),(100))<0.3;
	soma_soma_iSYN_netcon = soma_soma_iSYN_netcon.*(1-eye((100)));
	dendrite_soma_iCOM_compcon = eye((100),(100));
	dendrite_iMultiPoissonExp_Ge = multi_Poisson((100),(127),(2),(10),(1),2,.5,(2000),(0.01));
	dendrite_iMultiPoissonExp_Gi = multi_Poisson((100),(73),(2),(10),(1),5,.5,(2000),(0.01));
	soma_dendrite_iCOM_compcon = eye((100),(100));
	dendrite_dendrite_iGAP_UB = max((100),(100));
	dendrite_dendrite_iGAP_Xpre = linspace(1,dendrite_dendrite_iGAP_UB,(100))'*ones(1,(100));
	dendrite_dendrite_iGAP_Xpost = (linspace(1,dendrite_dendrite_iGAP_UB,(100))'*ones(1,(100)))';
	dendrite_dendrite_iGAP_mask = rand((100),(100))<0.3;
	dendrite_dendrite_iGAP_mask = dendrite_dendrite_iGAP_mask.*(1-eye((100)));

% Anonymous functions:
	soma_K_an            = @(soma_v) .01*(soma_v+55)./(1-exp(-(soma_v+55)/10));
	soma_K_bn            = @(soma_v) .125*exp(-(soma_v+65)/80);    
	soma_K_ik            = @(soma_v,soma_K_n) (36)*(soma_K_n.^4).*(soma_v-(-77));
	soma_Na_am           = @(soma_v) .1*(soma_v+40)./(1-exp(-(soma_v+40)/10));
	soma_Na_bm           = @(soma_v) 4*exp(-(soma_v+65)/18);       
	soma_Na_ah           = @(soma_v) .07*exp(-(soma_v+65)/20);     
	soma_Na_bh           = @(soma_v) 1./(1+exp(-(soma_v+35)/10));  
	soma_Na_ina          = @(soma_v,soma_Na_m,soma_Na_h) (120)*soma_Na_h.*(soma_Na_m.^3).*(soma_v-(50));
	soma_leak_ileak      = @(soma_v) (0.3)*(soma_v-(-54.4));       
	soma_soma_iSYN_ISYN  = @(V,soma_soma_iSYN_s) ((0).*(soma_soma_iSYN_netcon*soma_soma_iSYN_s).*(V-(0)));
	dendrite_soma_iCOM_dV = @(IN,OUT) ((IN*ones(1,size(IN,1)))'-(OUT*ones(1,size(OUT,1))));
	dendrite_soma_iCOM_ICOM = @(IN,OUT) (0.15).*sum((((IN*ones(1,size(IN,1)))'-(OUT*ones(1,size(OUT,1))))).*dendrite_soma_iCOM_compcon,2);
	dendrite_K_an        = @(dendrite_v) .01*(dendrite_v+55)./(1-exp(-(dendrite_v+55)/10));
	dendrite_K_bn        = @(dendrite_v) .125*exp(-(dendrite_v+65)/80);
	dendrite_K_ik        = @(dendrite_v,dendrite_K_n) (3.6)*(dendrite_K_n.^4).*(dendrite_v-(-77));
	dendrite_Na_am       = @(dendrite_v) .1*(dendrite_v+40)./(1-exp(-(dendrite_v+40)/10));
	dendrite_Na_bm       = @(dendrite_v) 4*exp(-(dendrite_v+65)/18);
	dendrite_Na_ah       = @(dendrite_v) .07*exp(-(dendrite_v+65)/20);
	dendrite_Na_bh       = @(dendrite_v) 1./(1+exp(-(dendrite_v+35)/10));
	dendrite_Na_ina      = @(dendrite_v,dendrite_Na_m,dendrite_Na_h) (12)*dendrite_Na_h.*(dendrite_Na_m.^3).*(dendrite_v-(50));
	dendrite_input_I     = @(t) 10;                                
	dendrite_leak_ileak  = @(dendrite_v) (0.03)*(dendrite_v-(-54.4));
	dendrite_iMultiPoissonExp_Gte = @(t) (0.1).*dendrite_iMultiPoissonExp_Ge(:,max(1,round(t/(0.01))));
	dendrite_iMultiPoissonExp_Gti = @(t) (0.1).*dendrite_iMultiPoissonExp_Gi(:,max(1,round(t/(0.01))));
	dendrite_iMultiPoissonExp_Itrain_e = @(t,dendrite_v) ((0.1).*dendrite_iMultiPoissonExp_Ge(:,max(1,round(t/(0.01))))).*(dendrite_v - (0));
	dendrite_iMultiPoissonExp_Itrain_i = @(t,dendrite_v) ((0.1).*dendrite_iMultiPoissonExp_Gi(:,max(1,round(t/(0.01))))).*(dendrite_v - (-85));
	dendrite_iMultiPoissonExp_Itrain = @(t,dendrite_v) (((0.1).*dendrite_iMultiPoissonExp_Ge(:,max(1,round(dendrite_v/(0.01))))).*(dendrite_v - (0))) + (((0.1).*dendrite_iMultiPoissonExp_Gi(:,max(1,round(dendrite_v/(0.01))))).*(dendrite_v - (-85)));
	soma_dendrite_iCOM_dV = @(IN,OUT) ((IN*ones(1,size(IN,1)))'-(OUT*ones(1,size(OUT,1))));
	soma_dendrite_iCOM_ICOM = @(IN,OUT) (0.15).*sum((((IN*ones(1,size(IN,1)))'-(OUT*ones(1,size(OUT,1))))).*soma_dendrite_iCOM_compcon,2);
	dendrite_dendrite_iGAP_IGAP = @(V1,V2) (0.5).*sum(((V1*ones(1,size(V1,1)))'-(V2*ones(1,size(V2,1)))).*dendrite_dendrite_iGAP_mask,2);

% ODE Handle, ICs, integration, and plotting:
ODEFUN = @(t,X) [((-((36)*(X(101:200).^4).*(X(1:100)-(-77))))+((-((120)*X(301:400).*(X(201:300).^3).*(X(1:100)-(50))))+((-((0.3)*(X(1:100)-(-54.4))))+((-(((0).*(soma_soma_iSYN_netcon*X(401:500)).*(X(1:100)-(0)))))+((((0.15).*sum((((X(501:600)*ones(1,size(X(501:600),1)))'-(X(1:100)*ones(1,size(X(1:100),1))))).*dendrite_soma_iCOM_compcon,2)))+0)))))/(1);(.01*(X(1:100)+55)./(1-exp(-(X(1:100)+55)/10))).*(1-X(101:200))-(.125*exp(-(X(1:100)+65)/80)).*X(101:200);(.1*(X(1:100)+40)./(1-exp(-(X(1:100)+40)/10))).*(1-X(201:300))-(4*exp(-(X(1:100)+65)/18)).*X(201:300);(.07*exp(-(X(1:100)+65)/20)).*(1-X(301:400))-(1./(1+exp(-(X(1:100)+35)/10))).*X(301:400);-X(401:500)./(1) + ((1-X(401:500))/(0.25)).*(1+tanh(X(1:100)/10));((-((3.6)*(X(601:700).^4).*(X(501:600)-(-77))))+((-((12)*X(801:900).*(X(701:800).^3).*(X(501:600)-(50))))+(((10))+((-((0.03)*(X(501:600)-(-54.4))))+((-((((0.1).*dendrite_iMultiPoissonExp_Ge(:,max(1,round(X(501:600)/(0.01))))).*(X(501:600) - (0))) + (((0.1).*dendrite_iMultiPoissonExp_Gi(:,max(1,round(X(501:600)/(0.01))))).*(X(501:600) - (-85)))))+((((0.15).*sum((((X(1:100)*ones(1,size(X(1:100),1)))'-(X(501:600)*ones(1,size(X(501:600),1))))).*soma_dendrite_iCOM_compcon,2)))+((((0.5).*sum(((X(501:600)*ones(1,size(X(501:600),1)))'-(X(501:600)*ones(1,size(X(501:600),1)))).*dendrite_dendrite_iGAP_mask,2)))+0)))))))/(1);(.01*(X(501:600)+55)./(1-exp(-(X(501:600)+55)/10))).*(1-X(601:700))-(.125*exp(-(X(501:600)+65)/80)).*X(601:700);(.1*(X(501:600)+40)./(1-exp(-(X(501:600)+40)/10))).*(1-X(701:800))-(4*exp(-(X(501:600)+65)/18)).*X(701:800);(.07*exp(-(X(501:600)+65)/20)).*(1-X(801:900))-(1./(1+exp(-(X(501:600)+35)/10))).*X(801:900);];
IC = [-65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.317        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.052        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596        0.596];

[t,y]=ode23(ODEFUN,[0 100],IC);   % numerical integration
figure; plot(t,y);           % plot all variables/functions
try legend('soma\_v','soma\_K\_n','soma\_Na\_m','soma\_Na\_h','soma\_soma\_iSYN\_s','dendrite\_v','dendrite\_K\_n','dendrite\_Na\_m','dendrite\_Na\_h'); end
%-
%}
