%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Step 4: 'S4_k_cluster_threshold.m' : compute activation threshold through a k-mean clustering approach.
	input: -'Info_IX.mat': all matlab structures used to compute an activation threshold.
	output: -'k_XX': resulting clustering groups
			'k_mean': mean fluorescence signal of each clustering group
			'k_min': min fluorescence signal of each clustering group
		-'Ik': final image where each pixel values corresponds to that of the cluster it has been assigned.
	dependencies: 'k_cluster_seg.m'
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

proj = 0; %parameter set to 1 to analyze projection Info structures

info = whos;
del = zeros(size(info,1),1);
for ii = 1:size(info,1)
    if ~contains(info(ii).name,'Info_')
        del(ii) = 1;
    else
        pos = strfind(info(ii).name, '_I');
        info(ii).num = str2double(info(ii).name(pos+2:end));
    end
end
info(del==1) = [];

if proj == 0
    vals = [];
    for ii = 1:size(info,1)
        info_temp = eval(sprintf('%1$s', info(ii).name));
        for jj = 1:size(info_temp,2)
            temp = mean(info_temp(jj).r_NAc_af(:,:,info_temp(jj).select_planes),3,'omitnan');
            temp(isnan(temp)) = 0;
            vals = [vals,temp];
        end
    end
    max_val = prctile(vals(:),99);
    min_val = prctile(vals(:),1);

    [Ik,k_mean] = k_mean_seg(vals,'v_range',[min_val max_val], 'nk', 5);
end


%% NORM
if proj == 0
    vals_norm = [];
    for ii = 1:size(info,1)
        info_temp = eval(sprintf('%1$s', info(ii).name));
        for jj = 1:size(info_temp,2)
            temp= info_temp(jj).r_NAc_af_norm;
            temp(isnan(temp)) = 0;
            vals_norm = [vals_norm,temp];
        end
    end
    width = size(info_temp(jj).r_NAc_af_norm,2);

else
    vals_norm = [];
    for ii = 1:size(info,1)
        info_temp = eval(sprintf('%1$s', info(ii).name));
        temp_core = mean(info_temp(2).r_core(:,:,info_temp(jj).select_planes),3,'omitnan');
        temp_sum_core = sum(sum(temp_core,'omitnan'),'omitnan');
        temp_num_core = sum(sum(double(~isnan(temp_core))));
        temp_core(isnan(temp_core)) = 0;
        
        
        temp_shell =  mean(info_temp(2).r_shell(:,:,info_temp(jj).select_planes),3,'omitnan');
        temp_sum_shell = sum(sum(temp_shell,'omitnan'),'omitnan');
        temp_num_shell = sum(sum(double(~isnan(temp_shell))));
        temp_shell(isnan(temp_shell)) = 0;
      
        
        temp_final = temp_core+temp_shell;
        
        temp_final =  temp_final./((temp_sum_shell+temp_sum_core)/(temp_num_shell+temp_num_core));
        vals_norm = [vals_norm, temp_final];
    end
    width = size(info_temp(2).r_core,2);

end

max_val_norm = prctile(vals_norm(:),99);
min_val_norm = prctile(vals_norm(:),1);



[Ik_norm,k_norm_mean] = k_mean_seg(vals_norm,'v_range',[min_val_norm max_val_norm], 'nk', 5);
k_norm_min = [0,(k_norm_mean(2:5)+k_norm_mean(1:4))./2];



num_im = size(info,1);
height= size(Ik_norm,1)*3;
if proj ==0
    stp = 0;
    count = 1;
    while stp == 0
        try
            Ik_norm(:, (count*width)+1:(count+1)*width) = [];
            count = count+1;
        catch
            stp = 1;
        end
    end
end
Ik_norm_amp = [Ik_norm,zeros(size(Ik_norm,1), (3-rem(num_im,3))*width)];
total_width = size(Ik_norm_amp,2)/3;

Ik_norm_amp = [Ik_norm_amp(:,1:total_width); ...
Ik_norm_amp(:,total_width+1:2*total_width); ...
Ik_norm_amp(:,2*total_width+1:end)];


figure
imagesc(reshape(Ik_norm_amp,height,[]))
