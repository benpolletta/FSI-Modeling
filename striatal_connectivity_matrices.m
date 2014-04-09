function [CS,CG] = striatal_connectivity_matrices(no_cells,plot_opt)

% Constructs synaptic and gap junctional connectivity matrices for FSI networks,
% according to percentages of connections reported in Russo & Taverna 2012.
% OUTPUT: CS = synaptic connectivity matrix. CG = gap-junction connectivity
% matrix.
% INPUT: no_cells = number of cells.

% Pairs of cell indices.
forward_pairs = nchoosek(1:no_cells,2);
no_pairs = size(forward_pairs,1);

reverse_pairs = fliplr(forward_pairs);

%% Parameters.

% Proportions of cell pairs having gap junctions, chemical synapses, or
% both.
p_both = .33;
p_gj_only = .18;
p_syn_only = .09;

p_gj = p_both + p_gj_only;

% Proportions of chemical synapses which are bidirectional.
p_both_bi = .23;
p_syn_only_bi = .31;

% Stats for gap junction conductances.
mean_gj = [1.7; 1.4];
var_gj = [0.2; 0.1];
mean_gj_diff = 0.48;
var_gj_diff = 0.07;

%% Gap Junctions.
% Choosing pairs with gap junctions first.
gj_randos = rand(no_pairs,1);
gj_indicator = gj_randos < p_gj;

% Choosing which direction has larger conductance.
gj_direction_randos = rand(no_pairs,1);
gj_forward_indicator = gj_indicator & (gj_direction_randos < .5); % Half are forward.
gj_reverse_indicator = gj_indicator & ~gj_forward_indicator; % Half are reverse.

% Collecting pairs of cells involved.
big_gj_pairs = [forward_pairs(gj_forward_indicator,:); reverse_pairs(gj_reverse_indicator,:)];
small_gj_pairs = fliplr(big_gj_pairs);
no_gjs = size(big_gj_pairs,1);

% Choosing gap junction conductances as multivariate normal.
covar_large_small = (var_gj_diff^2 - norm(var_gj))/2;
g_gj = mvnrnd(mean_gj, diag(var_gj) + [0 covar_large_small; covar_large_small 0], no_gjs);

% Making connectivity matrix.
CG = zeros(no_cells);

for g = 1:no_gjs
    CG(big_gj_pairs(g,2),big_gj_pairs(g,1)) = g_gj(g,1);
    CG(small_gj_pairs(g,2),small_gj_pairs(g,1)) = g_gj(g,2);
end

%% Chemical Synapses.
% Choosing pairs of synapses, conditional on whether a pair exhibits a gap
% junction or not.
syn_randos = rand(no_pairs,1);
p_syn_given_no_gj = p_syn_only/(1-p_gj);
syn_only_indicator = ~gj_indicator & (syn_randos < p_syn_given_no_gj); % Those exhibiting only chemical synapses are chosen from cells without GJs.
p_syn_given_gj = p_both/p_gj;
both_indicator = gj_indicator & (syn_randos < p_syn_given_gj); % Those exhibiting both are chosen from cells with GJs.
syn_indicator = syn_only_indicator | both_indicator; % Population of pairs combined to get all synaptic connections.

% Choosing which chemical synapses (in the absence of electrical synapses)
% are bidirectional.
syn_only_randos = rand(no_pairs,1);
syn_only_bi_indicator = syn_only_indicator & (syn_only_randos < p_syn_only_bi);

% Choosing which chemical synapses (in the presence of electrical synapses)
% are bidirectional.
both_randos = rand(no_pairs,1);
both_bi_indicator = both_indicator & (both_randos < p_both_bi);

% Populations of bidirectional synapse pairs combined to get all bidirectional synapses.
syn_bi_indicator = syn_only_bi_indicator | both_bi_indicator;

% Choosing, from unidirectional synapses, which are forward and which
% are reverse. (Not sure if this is necessary, but better safe than sorry.)
syn_direction_randos = rand(no_pairs,1);
syn_forward_indicator = syn_indicator & ~syn_bi_indicator & (syn_direction_randos < .5);
syn_reverse_indicator = syn_indicator & ~syn_bi_indicator & ~syn_forward_indicator;

% Retrieving which pairs of cells are involved.
syn_pairs = [forward_pairs(syn_bi_indicator | syn_forward_indicator,:); reverse_pairs(syn_bi_indicator | syn_reverse_indicator,:)];
no_syns = size(syn_pairs,1);

% Making connectivity matrix.
CS = zeros(no_cells);

g_s = 10*rand(no_syns,1) + 5; % Random synaptic conductance between 5 and 15 nS.
for s = 1:no_syns
    CS(syn_pairs(s,2),syn_pairs(s,1)) = g_s(s);
end

%% Analyze Output.

if plot_opt > 0
    
    figure(),
    subplot(1,2,1), imagesc(CS), title('Synaptic Connectivity Matrix')
    subplot(1,2,2), imagesc(CG), title('Gap Junction Connectivity Matrix')
    
    m_p_gj = sum(sum(CG>0))/(100^2);
    m_p_syn = sum(sum(triu(CS>0 | CS'>0)))/nchoosek(100,2);
    m_p_both = sum(sum(triu(CG>0 & (CS>0 | CS'>0))))/nchoosek(100,2);
    m_p_syn_only = sum(sum(triu(CG==0 & (CS>0 | CS'>0))))/nchoosek(100,2);
    m_p_gj_only = sum(sum(triu(CG>0 & ~(CS>0 | CS'>0))))/nchoosek(100,2);
    m_mean_gj_diff = mean(mean(abs(CG-CG')));
    
    fprintf('\t Expected: \t Measured: \n')
    fprintf('p_gj: \t %f \t %f \n',p_gj,m_p_gj)
    fprintf('p_syn: \t %f \t %f \n',p_syn_only+p_both,m_p_syn)
    fprintf('p_both: \t %f \t %f \n',p_both,m_p_both)
    fprintf('p_syn_only: \t %f \t %f \n',p_syn_only,m_p_syn_only)
    fprintf('p_gj_only: \t %f \t %f \n',p_gj_only,m_p_gj_only)
    fprintf('mean_gj_diff: \t %f \t %f \n',mean_gj_diff,m_mean_gj_diff)
    
end
	