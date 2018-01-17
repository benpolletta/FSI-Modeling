
  [data, name] = Gutfreund_KCa_PRC_model(0, 6000, 0, 'I_app', 0.87, 'STPstim', 0:-0.1:-0.5, 'STPwidth', 30, 'STPshift', 0:15:600, 'STPkernelType', 25, 'STPonset', 1000);
  
  results = Gutfreund_KCa_PRC_plot(data, [], name);