
%%This script assumes that the map ("Eur") is already defined in the
%%workspace.  This consists of three variables of dimension Lx1, where L
%%is the number of markers on which the map is defined.  For the ith
%%marker, chr_Eur(i) records the chromosome number (1 to 22), pos_Eur(i)
%%records the physical position (in bases, but encoded as the base number
%%plus chr_Eur(i)*10^9 so as to keep the entire list monotonic), and
%%cM_Eur records the map position (in cM, but encoded as the genetic
%%distance plus chr_Eur(i)*10^3 so as to keep the entire list monotonic).  


%Make evenly spaced blocks...derived from Eur_map.m
block_starts_even_Eur = [];
block_ends_even_Eur = [];
%next_space = [];
i = 2;
while i < size(cM_Eur,1)
    for j = 1:size(cM_Eur,1)-i
        if cM_Eur(i+j) - cM_Eur(i+j-1) > 0.1 %0.01
            i=i+j;
            break
        end
        if pos_Eur(i+j)-pos_Eur(i) > 10^5
            block_starts_even_Eur = [block_starts_even_Eur; i];
            block_ends_even_Eur = [block_ends_even_Eur; i+j];
            i = i+j;
            break
        end
        if j == size(cM_Eur,1)-i
            i=i+1;
            break
        end
    end
end

blocks_left_even_Eur = zeros(size(block_starts_even_Eur));
blocks_right_even_Eur = zeros(size(block_starts_even_Eur));
for i = 1:size(block_starts_even_Eur)
    s = block_starts_even_Eur(i);
    for j = 1:s-1
        if chr_Eur(s-j) < chr_Eur(s)
            blocks_left_even_Eur(i) = s-j+1;
            break
        end
        if cM_Eur(s)-cM_Eur(s-j) > 0.1
            blocks_left_even_Eur(i) = s-j;
            break
        end
    end
    if s-j == 1
        blocks_left_even_Eur(i) = s-j;
    end
end
for i = 1:size(block_ends_even_Eur)
    s = block_ends_even_Eur(i);
    for j = 1:size(cM_Eur,1)-s
        if chr_Eur(s+j) > chr_Eur(s)
            blocks_right_even_Eur(i) = s+j-1;
            break
        end
        if cM_Eur(s+j)-cM_Eur(s) > 0.1
            blocks_right_even_Eur(i) = s+j;
            break
        end
    end
    if s+j == size(cM_Eur,1)
        blocks_right_even_Eur(i) = s+j;
    end
end    

max_phys = 10^6; %enforce to keep the simulations under control

%Pare down the list of cold spots here before making the simulation script

under_max = ((pos_Eur(blocks_right_even_Eur)-pos_Eur(blocks_left_even_Eur))) < max_phys;
block_starts_even_Eur = block_starts_even_Eur(under_max);
block_ends_even_Eur = block_ends_even_Eur(under_max);
blocks_right_even_Eur = blocks_right_even_Eur(under_max);
blocks_left_even_Eur = blocks_left_even_Eur(under_max);


