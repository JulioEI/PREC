%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Step 3 (optional): S3_plot_change_in_Info_files.m :code used to visualize change in images after substracting the autofluorescence signal
	input: 	-'Info_IX.mat' matlab structure containing experimental processed data (after Step 2)
		-'Info_autofluo.mat' (see Step 1). 
	params: -'orig_color_th': upper threshold colormap for original data
		-'norm_color_th': upper threshold colormap for normalized data

%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
 %WARNING: MANUALLY CHANGE VARIABLE NAME: Info_IX IN ALL THE
%SCRIPT (use control+F to replace)


orig_color_th = 5000;
norm_color_th = 3;
basal = mean(Info_I1(1).r_core_og(:,:,Info_I1(1).select_planes),3,'omitnan');
basal(isnan(basal)) = 0;
temp = mean(Info_I1(1).r_shell_og(:,:,Info_I1(1).select_planes),3,'omitnan');
temp(isnan(temp)) = 0;
basal = basal + temp;
optostim = mean(Info_I1(2).r_core_og(:,:,Info_I1(2).select_planes),3,'omitnan');
optostim(isnan(optostim)) = 0;
temp = mean(Info_I1(2).r_shell_og(:,:,Info_I1(2).select_planes),3,'omitnan');
temp(isnan(temp)) = 0;
optostim = optostim + temp;

basal_af = mean(Info_I1(1).r_core_af(:,:,Info_I1(1).select_planes),3,'omitnan');
basal_af(isnan(basal_af)) = 0;
temp = mean(Info_I1(1).r_shell_af(:,:,Info_I1(1).select_planes),3,'omitnan');
temp(isnan(temp)) = 0;
basal_af = basal_af + temp;

optostim_af = mean(Info_I1(2).r_core_af(:,:,Info_I1(2).select_planes),3,'omitnan');
optostim_af(isnan(optostim_af)) = 0;
temp = mean(Info_I1(2).r_shell_af(:,:,Info_I1(2).select_planes),3,'omitnan');
temp(isnan(temp)) = 0;
optostim_af = optostim_af + temp;

basal_norm = mean(Info_I1(1).r_core(:,:,Info_I1(1).select_planes),3,'omitnan');
basal_norm(isnan(basal_norm)) = 0;
temp = mean(Info_I1(1).r_shell(:,:,Info_I1(1).select_planes),3,'omitnan');
temp(isnan(temp)) = 0;
basal_norm = basal_norm + temp;
optostim_norm = mean(Info_I1(2).r_core(:,:,Info_I1(2).select_planes),3,'omitnan');
optostim_norm(isnan(optostim_norm)) = 0;
temp = mean(Info_I1(2).r_shell(:,:,Info_I1(2).select_planes),3,'omitnan');
temp(isnan(temp)) = 0;
optostim_norm = optostim_norm + temp;


figure

subplot(2,4,1)
imagesc(basal)
title('Original Image')
caxis([0,orig_color_th])
ylabel('Basal', 'Fontsize', 14)
subplot(2,4,2)
imagesc(basal_norm)
title('Normalized Image')
caxis([0,norm_color_th])
subplot(2,4,3)
imagesc(Info_autofluo(1).r_proj_og)
title('Autofluorescence image')
caxis([0,orig_color_th])
subplot(2,4,4)
imagesc(basal_af)
title('Norm Image after autofluo substraction')
caxis([0,norm_color_th])

subplot(2,4,5)
imagesc(optostim)
caxis([0,orig_color_th])
ylabel('Optostim', 'Fontsize', 14)
subplot(2,4,6)
imagesc(optostim_norm)
caxis([0,norm_color_th])
subplot(2,4,7)
imagesc(Info_autofluo(2).r_proj_og)
caxis([0,orig_color_th])
subplot(2,4,8)
imagesc(optostim_af)
caxis([0,norm_color_th])
