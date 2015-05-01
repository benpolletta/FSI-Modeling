function [inhib_gj_fr, inhib_gj_pairs,inhib_shared_fr, inhib_shared_pairs,gj_shared_fr, gj_shared_pairs] = inhib_tests(T0,no_cells,p_gj,max_j,p_inhib)
smallvars = 1;
max_i = 10;
max_k = 11;
n = (15*0.0053)/sqrt(no_cells);

%inhib conductance vs gj conductance vs firing rate, spike pairs. 
%inhib conductance = 0 to 3x highest gj conductance, in same steps as gj
inhib_gj_fr = zeros(max_i, max_k);
inhib_gj_pairs= zeros(max_i, max_k);
for i = (1:max_i)
    inhib_strength = i * n;
    [mean_pop_firing_rate, mean_pop_norm_spike_pairs, ~, ~, ~] = gj_uncorr_input(T0,no_cells,p_gj,max_j,smallvars,p_inhib,inhib_strength);
    inhib_gj_fr(i,:) = mean_pop_firing_rate;
    inhib_gj_pairs(i,:) = mean_pop_norm_spike_pairs;
end

%inhib conductance vs shared input vs firing rate, spike pairs
%gj conductance = whatever it is with correlated input tests
inhib_shared_fr = zeros(max_i, max_k);
inhib_shared_pairs= zeros(max_i, max_k);
for i = (1:max_i)
    inhib_strength = i * n;
    [mean_pop_firing_rate, mean_pop_norm_spike_pairs, ~, ~, ~] = gj_corr_input(T0,no_cells,p_gj,max_j,smallvars,p_inhib,inhib_strength);
    inhib_shared_fr(i,:) = mean_pop_firing_rate;
    inhib_shared_pairs(i,:) = mean_pop_norm_spike_pairs;
end

%might be interesting: gj conductance vs shared input vs firing rate, spike pairs
gj_shared_fr = zeros(max_i, max_k);
gj_shared_pairs= zeros(max_i, max_k);

inhib_strength = n;
for i = (1:max_i)
    gj_strength = i * n;
    [mean_pop_firing_rate, mean_pop_norm_spike_pairs, ~, ~, ~] = gj_corr_input(T0,no_cells,p_gj,max_j,smallvars,p_inhib,inhib_strength,gj_strength);
    gj_shared_fr(i,:) = mean_pop_firing_rate;
    gj_shared_pairs(i,:) = mean_pop_norm_spike_pairs;
end

str = ['inhib_test_data', num2str(T0), '_', num2str(no_cells),'_',num2str(p_gj),'_',num2str(max_j),'_',num2str(p_inhib),'.mat'];
save(str,'','-v7')
end