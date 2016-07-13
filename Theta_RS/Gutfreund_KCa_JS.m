function data = Gutfreund_KCa_JS(I_const, I_app, varargin)

% Set tau_fast = 7, look at I_app = 2.5, ..., 3.5 to see transition from
% subthreshold oscillations to intermittent spiking to continuous spiking.

vary_label = ''; vary_cell = cell(length(varargin)/2, 3);

if ~isempty(varargin)
    
    for a = 1:(length(varargin)/2)
        
        vary_label = [vary_label, sprintf('_%s_%fto%f', varargin{2*a - 1}, varargin{2*a}(1), varargin{2*a}(end))];
        
        vary_cell(a, :) = {'pop1', varargin{2*a - 1}, varargin{2*a}};
        
    end
    
end

model_eqns = ['dv/dt=I_const+I(t)+@current/Cm; Cm=.25; v(0)=-65;',...
    '{iNaP,iKsKDyn,iKDRGKDyn,iNaG,gleak,CaDyn,iCa,KDyn,iKCa,};',...
    sprintf('gKs=0.084; gNaP=0.025; gl=0.025; I_const=%f;', I_const),...    %  halfKs=-60; halfNaP=-60; gNaP=0.0125; 
    'tau_fast=5; tau_h=tau_fast; tau_m=tau_fast;',...
    'offset=0; Koffset=offset; Noffset=offset;',...     % gKDR=5/3; gNa=12.5/3; gl=0;
    sprintf('I(t)=I_app*(ton<t&t<toff)*(1+rand*.25); ton=500; toff=3500; I_app=%f;', I_app),...
    'monitor functions'];

if ~isempty(varargin)
    
    if strcmp(version('-release'), '2012a')
    
        data = SimulateModel(model_eqns, 'tspan', [0 4000], 'vary', vary_cell);
    
    else
        
        data = SimulateModel(model_eqns, 'tspan', [0 4000], 'vary', vary_cell, 'compile_flag', 1);
    
    end
    
else
    
    if strcmp(version('-release'), '2012a')
    
        data = SimulateModel(model_eqns, 'tspan', [0 4000]);
    
    else
    
        data=SimulateModel(model_eqns, 'tspan', [0 4000], 'compile_flag', 1);
        
    end
     
end
    
try 
    
    PlotData(data)

    save_as_pdf(gcf, [sprintf('Figures/gutfreund_KCa_JS_Iconst%f_Iapp%f', I_const, I_app), vary_label])

catch error
    
    display('PlotData failed:')
    
    display(error)

end