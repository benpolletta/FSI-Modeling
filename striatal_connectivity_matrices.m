function [CS,CG] = striatal_connectivity_matrices(no_cells)

% Constructs synaptic and gap junctional connectivity matrices for FSI networks,
% according to percentages of connections reported in Russo & Taverna 2012.
% OUTPUT: CS = synaptic connectivity matrix. CG = gap-junction connectivity
% matrix.
% INPUT: no_cells = number of cells.

% Proportions of cell pairs having gap junctions, chemical synapses, or
% both.
pct_both = .33;
pct_gj_only = .18;
pct_syn_only = .09;

pct_gj = pct_both + pct_gj_only;

% Proportions of connections which are bidirectional.
pct_gj_bi = .23;
pct_syn_bi = .31;

% Pairs of cell indices.
forward_pairs = nchoosek(1:no_cells,2);
no_pairs = size(forward_pairs,1);

reverse_pairs = fliplr(forward_pairs);

% Choosing pairs with gap junctions first.
gj_randos = rand(no_pairs,1);
gj_indicator = gj_randos < pct_gj;

% Choosing pairs of synapses, conditional on whether a pair exhibits a gap
% junction or not.
syn_randos = rand(no_pairs,1);
syn_only_indicator = ~gj_indicator & (syn_randos < pct_syn_only); % Those exhibiting only chemical synapses are chosen from cells without GJs.
both_indicator = gj_indicator & (syn_randos < pct_both); % Those exhibiting both are chosen from cells with GJs.
syn_indicator = syn_only_indicator | both_indicator; % Population of pairs combined to get all synaptic connections.

% Choosing which gap junctions are bidirectional.
gj_bi_randos = rand(no_pairs,1);
gj_bi_indicator = gj_indicator & (gj_bi_randos < pct_gj_bi);

% Choosing, from unidirectional gap junctions, which are forward and which
% are reverse. (Not sure if this is necessary, but better safe than sorry.)
gj_direction_randos = rand(no_pairs,1);
gj_forward_indicator = gj_indicator & ~gj_bi_indicator & (gj_direction_randos < .5); % Half of unidirectional are forward.
gj_reverse_indicator = gj_indicator & ~gj_bi_indicator & ~gj_forward_indicator; % Half of unidirectional are reverse.

% Choosing which gap junctions are bidirectional.
syn_bi_randos = rand(no_pairs,1);
syn_bi_indicator = syn_indicator & (syn_bi_randos < pct_syn_bi);

% Choosing, from unidirectional synapses, which are forward and which
% are reverse. (Not sure if this is necessary, but better safe than sorry.)
syn_direction_randos = rand(no_pairs,1);
syn_forward_indicator = syn_indicator & ~syn_bi_indicator & (syn_direction_randos < .5);
syn_reverse_indicator = syn_indicator & ~syn_bi_indicator & ~syn_forward_indicator;

% Retrieving which pairs of cells are involved.
gj_pairs = [forward_pairs(gj_bi_indicator | gj_forward_indicator,:); reverse_pairs(gj_bi_indicator | gj_reverse_indicator,:)];
syn_pairs = [forward_pairs(syn_bi_indicator | syn_forward_indicator,:); reverse_pairs(syn_bi_indicator | syn_reverse_indicator,:)];

% Making connectivity matrices.
CG = zeros(no_cells);
CS = zeros(no_cells);

for g = 1:size(gj_pairs,1)
    CG(gj_pairs(g,2),gj_pairs(g,1)) = 1;
end

for s = 1:size(syn_pairs,1)
    CS(syn_pairs(s,2),syn_pairs(s,1)) = 1;
end
    
end
	
	
	
	