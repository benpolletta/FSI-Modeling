function STO_phase_experiment

tic; [data, name] = Gutfreund_original(0,0,'PPwidth',[100 250 500 750 1000],...
    'PPstim',-[0 .25 .5 .75 1]*10^(-4),'kernel_type',[1 3 7 15],...
    'I_app',.125,'gNa',0,'gKDR',0,'gl',0); toc;





