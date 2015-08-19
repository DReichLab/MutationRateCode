
%%This function generates calibration data to match real test data.
%%Overall it's very similar to calibrate_simulated_data.m and/or
%%load_real_data.m, with only one slightly differerent input: now
%%cal_genomes_per_sample specifies the number of genomes' worth of
%%calibration data per curve per real genome (e.g., with 8 real genomes,
%%cal_genomes_per_sample = 5 means 40 cal genomes per curve).



function [mu_cal,cal,het_rate_cal] = calibrate_real_data_new(seed_base,filter_dir,prefixes,filter_suff,chromosome,max_missing,cal_values,mincount,maxcount,ms_path,dnum_data,dden_data,het_rate_data,theta_base,history,chr,phys,map,block_starts,blocks_left,blocks_right,cal_genomes_per_sample,cal_chr,cal_phys,cal_map,cal_blocks_left,cal_blocks_right)

%setup
max_d = 0.1; %max value of d for plotting
nbins = 60; %number of bins for plot
dgrid = 0:max_d/nbins:max_d;
%compute total physical length of map
total_phys = 0;
for n = 2:size(phys,1)
    if chr(n) == chr(n-1)
        total_phys = total_phys + phys(n) - phys(n-1);
    end
end
%compute total genetic length of map
total_cM_cal = 0;
for n = 2:size(cal_map,1)
    if cal_chr(n) == cal_chr(n-1)
        total_cM_cal = total_cM_cal + cal_map(n) - cal_map(n-1);
    end
end
r = total_cM_cal/total_phys/100;
num_curves = size(cal_values,2);
cal = zeros(1,size(dgrid,2),num_curves);

%%Adjust the pop size history to account for multiple sites per PSMC bin
psmc_intervals = 26; 
old_hist_times = zeros(psmc_intervals,1);
new_hist_times = zeros(psmc_intervals,1);
old_hist_sizes = zeros(psmc_intervals,1);
new_hist_sizes = zeros(psmc_intervals,1);
histvec = sscanf(history,'%s %f %f %f');
for m = 1:psmc_intervals
    old_hist_times(m) = histvec(6*(m-1)+4);
    old_hist_sizes(m) = histvec(6*m);
end
new_hist_times(1) = 2*(1-(1-(theta_base)*old_hist_times(1)/2*200)^0.01)/(2*theta_base);
new_hist_sizes(1) = old_hist_sizes(1)*new_hist_times(1)/old_hist_times(1);
for m = 2:psmc_intervals
    new_hist_times(m) = new_hist_times(m-1)+2*((1-(1-(theta_base)*(old_hist_times(m)+old_hist_times(m-1))/2*200)^0.01)/(2*theta_base)-new_hist_times(m-1));
    new_hist_sizes(m) = old_hist_sizes(m)*(new_hist_times(m)-new_hist_times(m-1))/(old_hist_times(m)-old_hist_times(m-1));  
end
history = '';
for m =  1:psmc_intervals
    history = sprintf('%s-eN %f %f ',history,new_hist_times(m),new_hist_sizes(m));
end

%run the simulation subfunction
for z = 1:num_curves
    cal(:,:,z) = simulate(cal_values(z)*10^-8,z);
end


%function to do the simulations
function curve = simulate(mu_input,z)    
    
%setup
dnum_cal = zeros(1,size(dgrid,2));
dden_cal = zeros(1,size(dgrid,2));
tot_segsites_cal = 0;
tot_sites_cal = 0;

%simulation loop - first loop over genomes for filtering
for p = 1:size(prefixes,2)
    filterfile = sprintf('%s/%s.%s',filter_dir,prefixes{p},filter_suff);
    load(filterfile,'den_bad','num_overlap','num_overlap_core', ...
         'mask_starts','mask_ends');

for i = 1:size(blocks_left,1) %loop over regions
if chromosome ~= chr(block_starts(i)) %for calibrating jackknife replicates
    len = phys(blocks_right(i))-phys(blocks_left(i)); %physical length of region    

    if num_overlap(i)/len <= max_missing && num_overlap_core(i)/100000 <= max_missing
    num_spots = cal_blocks_right(i)-cal_blocks_left(i); 

    %setup
    midp_phys = 50000+phys(block_starts(i));
    midp_cM = interp1q(phys(blocks_left(i):blocks_right(i)),map(blocks_left(i):blocks_right(i)),midp_phys);
    jack_den = 0; %no jackknife here - run separate instances of calibration
    jack_ind = 1;
    theta = theta_base*len;    
    rho = theta_base*len*r/mu_input;
    
    %%Simulating one sequence at a time for randomization purposes
    for j = 1:cal_genomes_per_sample
             
        %fixed seed value for each region
        if seed_base > 0
            rng(500*25000*seed_base+100*25000*z+100*i+cal_genomes_per_sample*p+j);
        end         
        seeds = randi(1000000,3,1); 
    
        command = sprintf('%s 2 1 -seeds %d %d %d -t %f -r %f %d %s -v %d',ms_path,seeds(1),seeds(2),seeds(3),theta,rho,len,history,num_spots);
        for s = 1:num_spots
            spot_left = cal_phys(cal_blocks_left(i)+s-1)-cal_phys(cal_blocks_left(i))+1;
            spot_right = cal_phys(cal_blocks_left(i)+s)-cal_phys(cal_blocks_left(i));
            spot_weight = (cal_map(cal_blocks_left(i)+s)-cal_map(cal_blocks_left(i)+s-1))/(cal_phys(cal_blocks_left(i)+s)-cal_phys(cal_blocks_left(i)+s-1))/r/100;
            spot = sprintf('%d %d %f',spot_left,spot_right,spot_weight);
            command = sprintf('%s %s',command,spot);
        end
        command = strcat(command,' | egrep ''positions:|segsites: 0''');

        %call msHOT    
        [~, result] = system(command);
        ends = strfind(result,char(10));
        
        datarow = sscanf(result(12:ends-1),'%f');    
        seg_phys = datarow*len+phys(blocks_left(i));

        %%Mask the bad blocks
        mask_hets = zeros(size(datarow));
        if isempty(datarow) == 0;
        for k = 1:size(datarow,1)
            overlap = mask_ends{i}(mask_starts{i} < seg_phys(k));
            overlap = overlap(overlap > seg_phys(k));
            if isempty(overlap) == 0
                mask_hets(k) = 1;
            end
        end
        datarow(mask_hets==1) = [];
        end
        seg_phys = datarow*len+phys(blocks_left(i));

        cold_pos = sum((seg_phys-midp_phys < 50000 & seg_phys-midp_phys >= 0));
        if isempty(cold_pos) == 1
            cold_pos = 0;
        end
        cold_neg = sum((seg_phys-midp_phys > -50000 & seg_phys-midp_phys < 0));
        if isempty(cold_neg) == 1
            cold_neg = 0;
        end
        coldcount = cold_pos + cold_neg;
        tot_segsites_cal = tot_segsites_cal + size(datarow,1);
        tot_sites_cal = tot_sites_cal + len - num_overlap(i);
        %check whether this region gets used or not in this genome
        if coldcount >= mincount*(100000-num_overlap_core(i))/100000 && ...
                coldcount <= maxcount*(100000-num_overlap_core(i))/100000

            jack_den(jack_ind) = jack_den(jack_ind) + 1;
            
            %numerator
            seg_cM = interp1q(phys(blocks_left(i):blocks_right(i)),map(blocks_left(i):blocks_right(i)),seg_phys);
            d = abs(midp_cM-seg_cM);
            if size(d,1) > 0
                for k = 1:size(d,1)
                    if d(k) < max_d
                        dnum_cal(jack_ind,floor(d(k)*nbins/max_d)+1) = dnum_cal(jack_ind,floor(d(k)*nbins/max_d)+1) + 1;
                    end
                end
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
            dden_cal(:,b) = dden_cal(:,b) + jack_den'*(bin_low_phys(b) - bin_low_phys(b+1));
        else
            dden_cal(:,b) = dden_cal(:,b) + jack_den'*(bin_low_phys(b) - phys(blocks_left(i)));
            break
        end
    end
    for b = 1:nbins
        if map(blocks_right(i))-midp_cM > max_d*b/nbins
            bin_high_phys(b+1) = interp1q(map(blocks_left(i):blocks_right(i)),phys(blocks_left(i):blocks_right(i)),midp_cM+max_d*b/nbins);
            dden_cal(:,b) = dden_cal(:,b) + jack_den'*(bin_high_phys(b+1) - bin_high_phys(b));
        else
            dden_cal(:,b) = dden_cal(:,b) + jack_den'*(phys(blocks_right(i)) - bin_high_phys(b));
            break
        end
    end

    %Denominator adjustment for un-called sections
      for k = 1:size(dgrid,2)
          dden_cal(:,k) = dden_cal(:,k) - jack_den'*den_bad(i,k);
      end

    end
end
end
end

%%Adjust asymptote and translate to match the first bin
het_rate_cal = tot_segsites_cal/tot_sites_cal;
curve = dnum_cal./dden_cal+(dnum_cal./dden_cal-dnum_cal(1)/dden_cal(1))*(het_rate_sim/het_rate_cal-1);

end

%%Interpolation estimate of mutation rate
x0 = 2.3;
binlen = size(dnum_sim,2)-1;
var = dnum_sim(:,1:end-1)./dden_sim(:,1:end-1).^2;
mu_cal = 0;
for jj = 1:1;
    mu_cal(jj) = lsqnonlin(@int_fun,x0);
end

    function F = int_fun(x)
        F = zeros(binlen,1);
        for ii = 1:binlen
            interp_vec = zeros(1,num_curves);
            for kk = 1:num_curves
                interp_vec(kk) = cal(jj,ii,kk);
            end
            F(ii) = dnum_sim(jj,ii)/dden_sim(jj,ii) - interp1(cal_values,interp_vec,x,'spline');
        end
        F = F./var(jj,:)'.^0.5;
    end

end
