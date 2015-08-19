
%%Description of input variables:
%%print = whether to write .psmcfa files and/or save variables.  Options
%%are 0 (don't save anything), 0.5 (save variables but don't print
%%.psmcfa), 1 (save variables and write .psmcfa for all jackknife
%%replicates), and 23 (save variables and only print .psmcfa for
%%full-data rep).
%%inter_var = whether to set different randomized mutation rates for
%%different regions.  If 0, not implemented; if greater than 0, then each
%%region's mutation rate is picked uniformly at random on the interval
%%[inter_var*10^-8 2*mu-inter_var*10^-8]. 
%%intra_var = whether to implement within-region rate variability.  If 0,
%%not implemented.  If -1, uses auxiliary file (next variable; e.g.,
%%diversity in African populations).  If positive, generates randomized
%%variability parametrized by intra_var value.
%%varfile = file with regional variability information (if desired)
%%mu = value of mutation rate for simulations (e.g., 2.5*10^-8)
%%N = population size for simulations (e.g., 10000)
%%sim_history = string of population size changes (e.g., '-eN 0.025 0.1
%%-eN 0.05 1')
%%mincount, maxcount = bounds for number of het sites per 100 kb (e.g.,
%%5,10)
%%ms_path = path for msHOT binary
%%out_dir = directory in which to save output
%%out_name = partial file name for output
%%chr, phys, map = base map variables (see evenly_spaced_blocks.m)
%%block_starts, blocks_left, blocks_right = coordinates of regions, as from
%%evenly_spaced_blocks.m
%%sim_reps = number of simulated genomes to create
%%sim_map = map to use for simulating the data (e.g., with perturbation
%%and pseudo-count applied)



function [dnum_sim,dden_sim,het_rate_sim] = simulate_data_basic(print,inter_var,intra_var,varfile,mu,N,sim_history,mincount,maxcount,ms_path,out_dir,out_name,chr,phys,map,block_starts,blocks_left,blocks_right,sim_reps,sim_map)

max_d = 0.1; %max value of d for plotting
nbins = 60; %number of bins for plot
dgrid = 0:max_d/nbins:max_d;
psmc_binsize = 100;

%compute total physical length of map
total_phys = 0;
for i = 2:size(phys,1)
    if chr(i) == chr(i-1)
        total_phys = total_phys + phys(i) - phys(i-1);
    end
end

%%For variable mutation rate simulations
if intra_var == -1    
    load('varfile');
    sum_SNPs = sum(sum(SNP_totals,2));
    sum_sites = sum(sum(site_totals,2));
end

%%For PSMC inference, only use one copy of each chromosome from among all simulated samples, chosen at random
gen_select = randi(sim_reps,1,22);

if print >= 1
fid_psmc = zeros(24-print,1);
for i = print:23
    psmc_name = sprintf('%s/%s_rep%d.psmcfa',out_dir,out_name,i);
    fid_psmc(i) = fopen(psmc_name, 'w');
end
end

%compute total genetic length of map
total_cM_sim = 0;
for i = 2:size(sim_map,1)
    if chr(i) == chr(i-1)
        total_cM_sim = total_cM_sim + sim_map(i) - sim_map(i-1);
    end
end

%initialize important variables
r = total_cM_sim/total_phys/100;
dnum_sim = zeros(23,size(dgrid,2));
dden_sim = zeros(23,size(dgrid,2));
tot_segsites_sim = 0;
tot_sites_sim = 0;

%%Simulated data loop
for i = 1:size(blocks_left,1) %loop over regions
    len = phys(blocks_right(i))-phys(blocks_left(i)); %physical length of region
    num_spots = blocks_right(i)-blocks_left(i); %number of "hotspots," i.e., changes in recombination rate within region 
    
    %%Variable mutation rate
    if inter_var == 0
        mu_run = mu; 
    else
        mu_run = 2*(mu-inter_var*10^-8)*rand(1)+inter_var*10^-8;
    end
    if intra_var > 0
        %intra-block variability
        mu_run = 2/(1+intra_var)*mu_run;      
        mu_boundary1 = rand;
        mu_boundary2 = rand;
        bound_low = min(mu_boundary1, mu_boundary2);
        bound_high = max(mu_boundary1, mu_boundary2);
        inside = randi([0,1]);
    elseif intra_var == -1
        ratios = [SNP_totals(i,1)/max([site_totals(i,1) 1]) SNP_totals(i,2)/max([site_totals(i,2) 1]) SNP_totals(i,3)/max([site_totals(i,3) 1])];
        for th = 1:3
            if site_totals(i,th) < 10000
                ratios(th) = sum_SNPs/sum_sites;
            end
        end
        mu_run = mu_run*max(ratios)/(sum_SNPs/sum_sites);
        bound_low = 1/3;
        bound_high = 2/3;
    end
    
    %defining msHOT command to run
    theta = 4*N*mu_run*len;
    rho = 4*N*r*len;
    seeds = randi(1000000,3,1);
    command = sprintf('%s 2 %d -seeds %d %d %d -t %f -r %f %d %s -v %d',ms_path,sim_reps,seeds(1),seeds(2),seeds(3),theta,rho,len,sim_history,num_spots);
    for s = 1:num_spots
        spot_left = phys(blocks_left(i)+s-1)-phys(blocks_left(i))+1;
        spot_right = phys(blocks_left(i)+s)-phys(blocks_left(i));
        spot_weight = (sim_map(blocks_left(i)+s)-sim_map(blocks_left(i)+s-1))/(phys(blocks_left(i)+s)-phys(blocks_left(i)+s-1))/r/100;
        spot = sprintf('%d %d %f',spot_left,spot_right,spot_weight);
        command = sprintf('%s %s',command,spot);
    end
    command = strcat(command,' | egrep ''positions:|segsites: 0''');
    
    %call msHOT
    [~, result] = system(command);
    ends = [0 strfind(result,char(10))]; %keep track of the line breaks in the msHOT output, which mark the boundaries between the output for each genome
    
    %bin setup
    midp_phys = 50000+phys(block_starts(i));
    jack_ind = 1:23;
    jack_ind(chr(block_starts(i))) = []; %for the jackknife
    midp_cM = interp1q(phys(blocks_left(i):blocks_right(i)),map(blocks_left(i):blocks_right(i)),midp_phys);
    jack_den = zeros(1,23);
    
    %aggregate the results
    for j = 1:size(ends,2)-1 %loop over genomes
        datarow = sscanf(result(12+ends(j):ends(j+1)-1),'%f'); 

     if intra_var > 0    
        %%Down-sample for variable mutation rate: random variability
        sample = zeros(size(datarow));
        if isempty(datarow) == 0;
        if inside == 1
            for k = 1:size(datarow,1)
                if datarow(k) > bound_low && datarow(k) < bound_high
                    keep = rand;
                    if keep > intra_var
                        sample(k) = 1;
                    end
                end
            end
        else
            for k = 1:size(datarow,1)
                if datarow(k) < bound_low || datarow(k) > bound_high
                    keep = rand;
                    if keep > intra_var
                        sample(k) = 1;
                    end
                end
            end
        end
        datarow(sample==1) = [];
        end
     end     
     
     if intra_var == -1    
        %%Down-sample for variable mutation rate: auxiliary diversity file
        sample = zeros(size(datarow));
        if isempty(datarow) == 0;
            for k = 1:size(datarow,1)
                keep = rand;
                if datarow(k) < bound_low
                    if keep > ratios(1)/max(ratios)
                        sample(k) = 1;
                    end
                elseif datarow(k) < bound_high
                    if keep > ratios(2)/max(ratios)
                        sample(k) = 1;
                    end
                else
                    if keep > ratios(3)/max(ratios)
                        sample(k) = 1;
                    end
                end
            end
            datarow(sample==1) = [];
        end
     end
     
        %%Write .psmcfa output for selected genomes
        if print >= 1
        for rep = print:23
        if gen_select(chr(block_starts(i))) == j && chr(block_starts(i)) ~= rep 
            hets = zeros(1,floor(len/psmc_binsize));
            base_pos = len*datarow;
            bin_num = floor(base_pos/psmc_binsize)+1;
            hets(bin_num) = 1;
            fid = fid_psmc(rep);
            fprintf(fid,'>ms');
            fprintf(fid,'%d',i);
            fprintf(fid,'\n');
            for k = 1:size(hets,2)
                fprintf(fid,'%d',hets(k));
            end
            fprintf(fid,'\n');
        end
        end
        end   
   
        seg_phys = datarow*len+phys(blocks_left(i)); %physical positions of het sites
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
        tot_segsites_sim = tot_segsites_sim + size(datarow,1);
        tot_sites_sim = tot_sites_sim + len;
        %check whether this region gets used or not in this genome
        if cold_pos > 0 && cold_neg > 0 && coldcount >= mincount && coldcount <= maxcount
            jack_den(jack_ind) = jack_den(jack_ind) + 1;
            
            %numerator
            seg_cM = interp1q(phys(blocks_left(i):blocks_right(i)),map(blocks_left(i):blocks_right(i)),seg_phys);
            d = abs(midp_cM-seg_cM);
            if size(d,1) > 0
                for k = 1:size(d,1)
                    if d(k) < max_d
                        dnum_sim(jack_ind,floor(d(k)*nbins/max_d)+1) = dnum_sim(jack_ind,floor(d(k)*nbins/max_d)+1) + 1;
                    end
                end
            end
        end
    end
    
    %denominator - includes all genomes at once
    bin_low_phys = zeros(nbins+1,1);
    bin_low_phys(1) = midp_phys;
    bin_high_phys = zeros(nbins+1,1);
    bin_high_phys(1) = midp_phys;
    for b = 1:nbins
        if midp_cM-map(blocks_left(i)) > max_d*b/nbins
            bin_low_phys(b+1) = interp1q(map(blocks_left(i):blocks_right(i)),phys(blocks_left(i):blocks_right(i)),midp_cM-max_d*b/nbins);
            dden_sim(:,b) = dden_sim(:,b) + jack_den'*(bin_low_phys(b) - bin_low_phys(b+1));
        else
            dden_sim(:,b) = dden_sim(:,b) + jack_den'*(bin_low_phys(b) - phys(blocks_left(i)));
            break
        end
    end
    for b = 1:nbins
        if map(blocks_right(i))-midp_cM > max_d*b/nbins
            bin_high_phys(b+1) = interp1q(map(blocks_left(i):blocks_right(i)),phys(blocks_left(i):blocks_right(i)),midp_cM+max_d*b/nbins);
            dden_sim(:,b) = dden_sim(:,b) + jack_den'*(bin_high_phys(b+1) - bin_high_phys(b));
        else
            dden_sim(:,b) = dden_sim(:,b) + jack_den'*(phys(blocks_right(i)) - bin_high_phys(b));
            break
        end
    end
   
end

het_rate_sim = tot_segsites_sim/tot_sites_sim;

%save variables
if print > 0    
var_file = sprintf('%s/sim_data_bins.mat',out_dir);
save(var_file,'dnum_sim','dden_sim','het_rate_sim');
end
if print >= 1
for i = print:23
    fclose(fid_psmc(i));
end
end

end

