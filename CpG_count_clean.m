
%%This function counts (in order of output variables) CpG and non-CpG
%%hets in ascertained regions, total CpG sites and total sites in
%%ascertained regions, and CpG hets and non-CpG hets genome-wide. The
%%structure and input variables are very similar to the data-loading
%%functions.  Here in addition to the hetfile and filterfile, the program
%%needs a sitefile input with two columns (chromosome and physical
%%position) for all CpG sites, as well as a badfile (as in
%%preprocess_bad_blocks.m). 


function [CpG,nonCpG,tot_CpG_sites,tot_sites,tot_CpG,tot_nonCpG,list] = CpG_count_clean(half_window,het_dir,filter_dir,het_suff,site_suff,filter_suff,prefixes,max_missing,mincount,maxcount,chr,phys,block_starts,blocks_left,blocks_right)

%%Currently set up to display all of the variables (except for "list")
%%but not to save them to disk 

CpG = zeros(23,1);
nonCpG = zeros(23,1);
tot_CpG_sites = zeros(23,1);
tot_sites = zeros(23,1);
tot_CpG = 0;
tot_nonCpG = 0;
list = []; %CpG hets, total hets, CpG sites, total sites

for j = 1:size(prefixes,2)
    %Load data
    hetfile  = sprintf('%s/%s.%s',het_dir,prefixes{j},het_suff);
    sitefile  = sprintf('%s/%s.%s',het_dir,prefixes{j},site_suff);
    filterfile = sprintf('%s/%s.%s',filter_dir,prefixes{j},filter_suff);
    genome_data = load(hetfile,'ascii');
    site_data = load(sitefile,'ascii');
    genome_pos = genome_data(:,2) + 10^9*genome_data(:,1);
    site_pos = site_data(:,2) + 10^9*site_data(:,1);
    CpG_list = genome_data(:,3);
    tot_CpG = tot_CpG + sum(CpG_list);
    tot_nonCpG = tot_nonCpG + sum(1-CpG_list);
    load(filterfile);
    good_blocks = ones(1,size(blocks_left,1));
    num_overlap_window = zeros(size(blocks_left,1),1);
    
    if half_window < 50000 %redo calculations as in preprocess_bad_blocks.m
        %make new count of filtered-out bases
        badfile = sprintf('%s/%s.bad_blocks.tab',het_dir,prefixes{j});
        bad_regions = load(badfile,'ascii');
        bad_starts = bad_regions(:,2) + 10^9*bad_regions(:,1);
        bad_ends = bad_regions(:,3) + 10^9*bad_regions(:,1);
    
        %extend bad blocks 2 bases in each direction
        for i = 2:size(bad_starts,1)-1
            bad_starts(i) = max(bad_starts(i)-2,bad_ends(i-1)+1);
            bad_ends(i) = min(bad_ends(i)+2,bad_starts(i+1)-1);
        end
    end

for i = 1:size(blocks_left,1)
        len = phys(blocks_right(i))-phys(blocks_left(i));
            if num_overlap(i)/len > max_missing
                good_blocks(i) = 0;
            end
            if num_overlap_core(i)/100000 > max_missing
                good_blocks(i) = 0;
            end

    if good_blocks(i) == 1;    
        if half_window < 50000
            window_starts = bad_starts(bad_starts < phys(block_starts(i))+50000+half_window);
            window_ends = bad_ends(bad_starts < phys(block_starts(i))+50000+half_window);
            window_starts = window_starts(window_ends > phys(block_starts(i))+50000-half_window);
            window_ends = window_ends(window_ends > phys(block_starts(i))+50000-half_window);
       
            if size(window_starts,1) > 0
                if window_starts(1) < phys(block_starts(i))+50000-half_window
                    window_starts(1) = phys(block_starts(i))+50000-half_window;
                end
                if window_ends(end) > phys(block_starts(i))+50000+half_window
                    window_ends(end) = phys(block_starts(i))+50000+half_window;
                end
                num_overlap_window(i) = sum(window_ends-window_starts+1);
            end
        end
        
        datarow = genome_pos((genome_pos > phys(blocks_left(i)))&(genome_pos < phys(blocks_right(i))));
        site_phys = site_pos((site_pos > phys(blocks_left(i)))&(site_pos < phys(blocks_right(i))));
        CpGrow = CpG_list((genome_pos > phys(blocks_left(i)))&(genome_pos < phys(blocks_right(i))));

        datarow = (datarow - phys(blocks_left(i)))/len;
        seg_phys = datarow*len+phys(blocks_left(i));        
        midp_phys = 50000+phys(block_starts(i));
          
        CpG_core = CpGrow((abs(seg_phys-midp_phys) < half_window)); %default would be half_window = 50K
        site_core = site_phys((abs(site_phys-midp_phys) < half_window));
        cold_pos = sum((seg_phys-midp_phys < 50000 & seg_phys-midp_phys >= 0));
        if isempty(cold_pos) == 1
            cold_pos = 0;
        end
        cold_neg = sum((seg_phys-midp_phys > -50000 & seg_phys-midp_phys < 0));
        if isempty(cold_neg) == 1
            cold_neg = 0;
        end
        coldcount = cold_pos + cold_neg;
        if cold_pos > 0 && cold_neg > 0 && coldcount >= mincount*(100000-num_overlap_core(i))/100000 && coldcount <= maxcount*(100000-num_overlap_core(i))/100000
            jack_ind = 1:23;
            jack_ind(chr(block_starts(i))) = [];
            CpG(jack_ind) = CpG(jack_ind) + sum(CpG_core);
            nonCpG(jack_ind) = nonCpG(jack_ind) + sum(1-CpG_core);
            tot_CpG_sites(jack_ind) = tot_CpG_sites(jack_ind) + size(site_core,1);
            if half_window < 50000
                tot_sites(jack_ind) = tot_sites(jack_ind) + 2*half_window - num_overlap_window(i);
            else
                tot_sites(jack_ind) = tot_sites(jack_ind) + 100000 - num_overlap_core(i);
            end
            list = [list; sum(CpG_core) sum(1-CpG_core)+sum(CpG_core) ...
                    size(site_core,1) 2*half_window-num_overlap_window(i)]; 
        end
    end
end
end
disp(CpG)
disp(nonCpG)
disp(tot_CpG_sites)
disp(tot_sites)
disp(tot_CpG)
disp(tot_nonCpG)
end




