
%%The input for this function is almost the same as that for
%%load_real_data.m, but with one extra variable CpG = 1/0 depending on
%%which category of het sites is to be used.  Also, the het file should
%%now have an extra third column with the same 1/0 coding for each site.



function [dnum_data,dden_data,het_rate_data] = load_real_data_CpG(CpG,save,het_dir,filter_dir,prefixes,het_suff,filter_suff,max_missing,mincount,maxcount,chr,phys,map,block_starts,blocks_left,blocks_right,data_reps)

max_d = 0.1; %max value of d for plotting
nbins = 60; %number of bins for plot
dgrid = 0:max_d/nbins:max_d;

%Initialization
dnum_data = zeros(data_reps+1,size(dgrid,2));
dden_data = zeros(data_reps+1,size(dgrid,2));
tot_segsites_data = 0;
tot_sites_data = 0;

%%Main loop - loop over genomes first
for j = 1:size(prefixes,2)
    %load data
    hetfile  = sprintf('%s/%s.%s',het_dir,prefixes{j},het_suff);
    filterfile = sprintf('%s/%s.%s',filter_dir,prefixes{j},filter_suff);
    genome_data = load(hetfile,'ascii');
    genome_pos = genome_data(:,2) + 10^9*genome_data(:,1);
    CpG_info = genome_data(:,3);
    genome_pos = genome_pos(CpG_info==CpG);
    load(filterfile);
    good_blocks = ones(1,size(blocks_left,1));

    %check missing data threshold    
for i = 1:size(blocks_left,1)
        len = phys(blocks_right(i))-phys(blocks_left(i));
            if num_overlap(i)/len > max_missing
                good_blocks(i) = 0;
            end
            if num_overlap_core(i)/100000 > max_missing
                good_blocks(i) = 0;
            end
        
    if good_blocks(i) == 1;    
        %put het data in proper format
        datarow = genome_pos((genome_pos > phys(blocks_left(i)))&(genome_pos < phys(blocks_right(i))));
        datarow = (datarow - phys(blocks_left(i)))/len;
        seg_phys = datarow*len+phys(blocks_left(i)); %physical positions of het sites
        midp_phys = 50000+phys(block_starts(i));
        midp_cM = interp1q(phys(blocks_left(i):blocks_right(i)),map(blocks_left(i):blocks_right(i)),midp_phys);
        %add up hets in each half of the region
        cold_pos = sum((seg_phys-midp_phys < 50000 & seg_phys-midp_phys >= 0));
        if isempty(cold_pos) == 1
            cold_pos = 0;
        end
        cold_neg = sum((seg_phys-midp_phys > -50000 & seg_phys-midp_phys < 0));
        if isempty(cold_neg) == 1
            cold_neg = 0;
        end
        coldcount = cold_pos + cold_neg;
        %hack to override the pos/neg check for CpGs given low density
        if CpG == 1
            cold_pos = 1;
            cold_neg = 1;
        end
        tot_segsites_data = tot_segsites_data + size(datarow,1);
        tot_sites_data = tot_sites_data + len - num_overlap(i);
        %check whether this region gets used or not in this genome
        if cold_pos > 0 && cold_neg > 0 && coldcount >= mincount*(100000-num_overlap_core(i))/100000 && coldcount <= maxcount*(100000-num_overlap_core(i))/100000
            %which replictes get the data - jackknife over chromosomes
            jack_ind = 1:data_reps+1;
            jack_ind(chr(block_starts(i))) = [];
            jack_den = zeros(1,data_reps+1);           
            jack_den(jack_ind) = jack_den(jack_ind) + 1;
            %numerator
            seg_cM = interp1q(phys(blocks_left(i):blocks_right(i)),map(blocks_left(i):blocks_right(i)),seg_phys);
            d = abs(midp_cM-seg_cM);
            if size(d,1) > 0
                for k = 1:size(d,1)
                    if d(k) < max_d
                        dnum_data(jack_ind,floor(d(k)*nbins/max_d)+1) = dnum_data(jack_ind,floor(d(k)*nbins/max_d)+1) + 1;
                    end
                end
            end
        
    
            %denominator 
            bin_low_phys = zeros(nbins+1,1);
            bin_low_phys(1) = midp_phys;
            bin_high_phys = zeros(nbins+1,1);
            bin_high_phys(1) = midp_phys;
            for b = 1:nbins
                if midp_cM-map(blocks_left(i)) > max_d*b/nbins
                    bin_low_phys(b+1) = interp1q(map(blocks_left(i):blocks_right(i)),phys(blocks_left(i):blocks_right(i)),midp_cM-max_d*b/nbins);
                    dden_data(:,b) = dden_data(:,b) + jack_den'*(bin_low_phys(b) - bin_low_phys(b+1));
                else
                    dden_data(:,b) = dden_data(:,b) + jack_den'*(bin_low_phys(b) - phys(blocks_left(i)));
                    break
                end
            end
            for b = 1:nbins
                if map(blocks_right(i))-midp_cM > max_d*b/nbins
                    bin_high_phys(b+1) = interp1q(map(blocks_left(i):blocks_right(i)),phys(blocks_left(i):blocks_right(i)),midp_cM+max_d*b/nbins);
                    dden_data(:,b) = dden_data(:,b) + jack_den'*(bin_high_phys(b+1) - bin_high_phys(b));
                else
                    dden_data(:,b) = dden_data(:,b) + jack_den'*(phys(blocks_right(i)) - bin_high_phys(b));
                    break
                end
            end
    
            %denominator adjustment for un-called sections
            for k = 1:size(dgrid,2)
                dden_data(jack_ind,k) = dden_data(jack_ind,k) - den_bad(i,k);
            end

        end
    end
end

end

het_rate_data = tot_segsites_data/tot_sites_data;

%save variables
if save == 1
var_file = sprintf('%s/data_bins_10to20.mat',out_dir);
save(var_file,'dnum_data','dden_data','het_rate_data');
end

end

