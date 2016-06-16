function data = Gutfreund(I_const, I_app, varargin)

vary_label = ''; vary_cell = cell(length(varargin)/2, 3);

if ~isempty(varargin)
    
    for a = 1:(length(varargin)/2)
        
        vary_label = [vary_label, sprintf('_%s_%fto%f', varargin{2*a - 1}, varargin{2*a}(1), varargin{2*a}(end))];
        
        vary_cell(a, :) = {'pop1', varargin{2*a - 1}, varargin{2*a}};
        
    end
    
end

model_eqns = ['dv/dt=I_const+I(t)+@current/Cm; Cm=.25; v(0)=-65;',...
    sprintf('{iNaP,iKs,iKDRG,iNaG,gleak}; halfKs=-60; halfNaP=-60; gNaP=0.0125; I_const=%f;', I_const),...
    'offset=-17; Koffset=offset; Noffset=offset; gdenom=2.5; gKDR=5/gdenom; gNa=12.5/gdenom; gl=0;',... % gdenom=3; gKDR=5/gdenom; gNa=12.5/gdenom; gl=0;'
    'tau_fast=3; tau_h=tau_fast; tau_m=tau_fast;',...
    sprintf('I(t)=I_app*(ton<t&t<toff)*(1+rand*.25); ton=500; toff=3500; I_app=%f;', I_app),...
    'monitor functions'];

if ~isempty(varargin)
    
    if strcmp(version('-release'), '2012a')
    
        data = SimulateModel(s, 'tspan', [0 4000], 'vary', vary_cell);
    
    else
        
        data = SimulateModel(s, 'tspan', [0 4000], 'vary', vary_cell, 'compile_flag', 1);
    
    end
    
else
    
    if strcmp(version('-release'), '2012a')
    
        data = SimulateModel(s, 'tspan', [0 4000]);
    
    else
    
        data=SimulateModel(s, 'tspan', [0 4000], 'compile_flag', 1);
        
    end
     
end
    
try 
    
    PlotData(data)

    save_as_pdf(gcf, [sprintf('gutfreund_Iconst%f_Iapp%f_', I_const, I_app), vary_label])

catch error
    
    display('PlotData failed:')
    
    display(error)

end