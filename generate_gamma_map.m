
%%The input variable alpha is the map variance parameter (in the standard
%%units of M^-1), while chr and map are genetic map variables as
%%described in evenly_spaced_blocks.m (in typical applications, the
%%pseudo-count prior is applied before this perturbation script is run).


function [gamma_map1] = generate_gamma_map(alpha,chr,map)
gamma_map1 = map;
for i = 2:size(map,1)
    if chr(i) == chr(i-1)
        gamma_map1(i) = gamma_map1(i-1) + gamrnd(alpha*(map(i)-map(i-1))/100,1/alpha)*100;
    end
end
end
