function data = Gutfreund_original(I_const, I_app, varargin)

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
    sprintf('{iNaP,iKs,iKDRG,iNaG,gleak}; gKs=0.084; gNaP=0.025; gl=0.025; I_const=%f;', I_const),...    %  halfKs=-60; halfNaP=-60; gNaP=0.0125;
    'tau_fast=7; tau_h=tau_fast; tau_m=tau_fast;',...
    'offset=0; Koffset=offset; Noffset=offset;',...     % gKDR=5/3; gNa=12.5/3; gl=0;
    sprintf('I(t)=I_app*(ton<t&t<toff)*(1+rand*.25); ton=500; toff=3500; I_app=%f;', I_app),...
    'monitor functions'];

if ~isempty(varargin)
    
    data=SimulateModel(model_eqns, 'tspan', [0 4000], 'vary', vary_cell, 'compile_flag', 1);
    
else
    
    data=SimulateModel(model_eqns, 'tspan', [0 4000], 'compile_flag', 1); %{'pop1','gKs',.04:.002:.06 % ;'pop1','gNaP',.010:.001:.020}); % {'pop1','gd',5:10;'pop1','I_app',10:20});

end
    
try 
    
    PlotData(data)

    save_as_pdf(gcf, [sprintf('gutfreund_Iconst%f_Iapp%f_', I_const, I_app), vary_label])

catch error
    
    display('PlotData failed:')
    
    display(error)

end