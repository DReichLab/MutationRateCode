
%%The input variables are:
%%prior = pseudo-count prior in units of cM/kb (e.g., 0.00009)
%%chr, pos, map = base map variables (see note in evenly_spaced_blocks.m)
%%type = 'Eur' when using shared AA map and 'decode' for DECODE map


function [cM_Eur_rescaled] = scale_Eur_map(prior,chr,pos,map,type)
cM_Eur_rescaled = map;
for i = 2:size(map,1)
    if chr(i) == chr(i-1)
        diff = map(i) - map(i-1);
        diff_ps = pos(i) - pos(i-1);
        if strcmp(type,'Eur')
            cM_Eur_rescaled(i) = cM_Eur_rescaled(i-1) + (29/(29+2.55*prior/0.0001))*(diff+prior*max(0,diff_ps)/1000);
        else
            cM_Eur_rescaled(i) = cM_Eur_rescaled(i-1) + (31.8/(31.8+2.67*prior/0.0001))*(diff+prior*max(0,diff_ps)/1000);
        end
    end
end
end
