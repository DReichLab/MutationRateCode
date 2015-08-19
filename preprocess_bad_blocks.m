
%%This script takes in a "bad_blocks" file
%%[data_dir]/[prefix].bad_blocks.tab and processes it into variables used
%%for filtering the data.  The .tab file should have three columns; each
%%segment to be filtered out has one line, with the chromosome number in
%%the first column and the start and end (physical) positions (inclusive)
%%in the second and third columns.  The variable "extend" (typically set
%%to 2) specifies whether to extend the raw segments (and by how many
%%bases in either direction).  The results are saved in the out_dir folder.


function [den_bad,num_overlap,num_overlap_core,mask_starts,mask_ends] = preprocess_bad_blocks(extend,data_dir,prefix,out_dir,chr,phys,map,block_starts,blocks_left,blocks_right)

max_d = 0.1; 
nbins = 60;
dgrid = 0:max_d/nbins:max_d;

%initialization
den_bad = zeros(size(blocks_left,1),size(dgrid,2));
num_overlap = zeros(size(blocks_left,1),1);
num_overlap_core = zeros(size(blocks_left,1),1);
mask_starts = {};
mask_ends = {};

    %load data
    badfile = sprintf('%s/%s.bad_blocks.tab',data_dir,prefix);
    bad_regions = load(badfile,'ascii');
    bad_starts = bad_regions(:,2) + 10^9*bad_regions(:,1);
    bad_ends = bad_regions(:,3) + 10^9*bad_regions(:,1);
    
    %%Extend bad blocks a few bases in each direction
    if extend > 0
        for i = 2:size(bad_starts,1)-1
            bad_starts(i) = max(bad_starts(i)-extend,bad_ends(i-1)+1);
            bad_ends(i) = min(bad_ends(i)+extend,bad_starts(i+1)-1);
        end
    end

    for i = 1:size(blocks_left,1)
        overlap_starts = bad_starts(bad_starts < phys(blocks_right(i)));
        overlap_ends = bad_ends(bad_starts < phys(blocks_right(i)));
        overlap_starts = overlap_starts(overlap_ends > phys(blocks_left(i)));
        overlap_ends = overlap_ends(overlap_ends > phys(blocks_left(i)));
        core_starts = bad_starts(bad_starts < phys(block_starts(i))+100000);
        core_ends = bad_ends(bad_starts < phys(block_starts(i))+100000);
        core_starts = core_starts(core_ends > phys(block_starts(i)));
        core_ends = core_ends(core_ends > phys(block_starts(i)));
        if size(overlap_starts,1) > 0
            if overlap_starts(1) < phys(blocks_left(i))
                overlap_starts(1) = phys(blocks_left(i));
            end
            if overlap_ends(end) > phys(blocks_right(i))
                overlap_ends(end) = phys(blocks_right(i));
            end
            num_overlap(i) = sum(overlap_ends-overlap_starts+1);
        end
        if size(core_starts,1) > 0
            if core_starts(1) < phys(block_starts(i))
                core_starts(1) = phys(block_starts(i));
            end
            if core_ends(end) > phys(block_starts(i))+100000
                core_ends(end) = phys(block_starts(i))+100000;
            end
            num_overlap_core(i) = sum(core_ends-core_starts+1);
        end
        mask_starts{i} = overlap_starts;
        mask_ends{i} = overlap_ends;
        
            %denominator adjustment for un-called sections
            overlap_bases = [];
            midp_phys = 50000+phys(block_starts(i));
            midp_cM = interp1q(phys(blocks_left(i):blocks_right(i)),map(blocks_left(i):blocks_right(i)),midp_phys);
            for k=1:size(overlap_starts,1)
                overlap_bases = [overlap_bases overlap_starts(k):overlap_ends(k)];
            end
            overlap_cM = interp1q(phys(blocks_left(i):blocks_right(i)),map(blocks_left(i):blocks_right(i)),overlap_bases');
            d = abs(midp_cM-overlap_cM);
            if size(d,1) > 0
                for k = 1:size(d,1)
                    if d(k) < max_d
                        den_bad(i,floor(d(k)*nbins/max_d)+1) = den_bad(i,floor(d(k)*nbins/max_d)+1) + 1;
                    end
                end
            end
        
        
    end

var_file = sprintf('%s/%s.filter_even_extend2_eur.mat',out_dir,prefix);
save(var_file,'den_bad','num_overlap','num_overlap_core','mask_starts','mask_ends');

end

