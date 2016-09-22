function STO_phase_experiment

I_app = 0.2;
tauCa = 300;
gKCa = 0.005;
aKCa_factor = 0.1;
CAF = 24;
pulse = 0;

[data, name] = Gutfreund_original(0, 0, 'I_app', I_app, 'pulse', pulse);



