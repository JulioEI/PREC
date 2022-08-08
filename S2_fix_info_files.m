%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Step 2: S2_fix_info_files.m: substract reference autofluorescence signal from experimental Info_IX structures.

	input: 	-'Info_IX.mat': matlab structure containing experimental data (see Step 0).
		-'Info_autofluo.mat' (see Step 1). 
	params: intensity difference between control images (used to compute autofluorescence signal) and
		experimental ones (computed ad hoc by matching background signals outside NAc, default: 0.4)
	output: 'Info_IX.mat' updated matlab structure with the following new fields:
			-'r_NAc_af': red signal on each grid of the NAc after accounting for autofluorescence signal
			-'r_core_af': red signal on each grid on the core region of the NAc after accounting for autofluoresnce signal	
			-'r_shell_af': red signal on each grid on the shell region of the NAc after accounting for autofluoresnce signal
		'IX_Basal.tif' image containing processed grid basal matrix
		'IX_TBS.tif' image containing processed grid tbs matrix	

%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%substract autofluo ref to Info
Info = Info_I1;
ratio = 0.4;
%create prec pre-normalization
for ii = 1:2
    Info(ii).r_core_og = Info(ii).r_core;
    Info(ii).r_shell_og = Info(ii).r_shell;
end
%compute prec pre-normalization
for jj = 1:2
    for ii = 1:10
        Info(jj).r_core_og(:,:,ii) = Info(jj).r_core(:,:,ii).*Info(jj).r_bg(ii,1);
        Info(jj).r_shell_og(:,:,ii) = Info(jj).r_shell(:,:,ii).*Info(jj).r_bg(ii,1);
    end
end
%find ratio such that BG_autofluo = 0.4BG & substrat BG_autofluo &
%recompute prec norm
for jj = 1:2
    Info(jj).afluo_ratio = ratio.* Info(jj).r_bg./Info_autofluo(jj).r_bg;
    %add core and shell before removing auto and then split again (avoid
    %edge problem)
    NAc = Info(jj).r_core_og;
    NAc(isnan(NAc)) = 0;
    temp = Info(jj).r_shell_og;
    temp(isnan(temp)) = 0;
    NAc = NAc + temp;
    %define core mask and shell mask to divide
    core_mask = double(~isnan(sum(Info(jj).r_core_og,3)));
    shell_mask = double(~isnan(sum(Info(jj).r_shell_og,3)));
    NAc_mask = core_mask+shell_mask;
    NAc_mask(NAc_mask == 0) = NaN;
    NAc_mask(~isnan(NAc_mask)) = 1;

    shell_mask(shell_mask~=1) = NaN;
    core_mask(shell_mask==1) = NaN;
    core_mask(core_mask~=1) = NaN;
    %compute new images
    Info(jj).r_core_af = nan(size(Info(jj).r_core_og));
    Info(jj).r_shell_af = nan(size(Info(jj).r_shell_og));
    Info(jj).r_NAc_af = nan(size(Info(jj).r_shell_og));
    for plane = 1:size(Info(jj).r_core,3)
        NAc(:,:,plane) = (NAc(:,:,plane) - ...
            Info(jj).afluo_ratio(plane,1)*Info_autofluo(jj).r_proj_og)./...
            abs(Info(jj).r_bg(plane,1));%-Info(jj).afluo_ratio(plane,1)*Info_autofluo(jj).r_bg);
        
        Info(jj).r_NAc_af(:,:,plane) = NAc(:,:,plane).*NAc_mask;
        Info(jj).r_core_af(:,:,plane) = NAc(:,:,plane).*core_mask;
        Info(jj).r_shell_af(:,:,plane) = NAc(:,:,plane).*shell_mask;
    end
    Info(jj).r_core_af(Info(jj).r_core_af<0) = NaN;
    Info(jj).r_shell_af(Info(jj).r_shell_af<0) = NaN;
    Info(jj).r_NAc_af(Info(jj).r_NAc_af<0) = NaN;
    
    if strcmpi(Info(jj).prot, 'basal')
        temp = mean(Info(jj).r_NAc_af(:,:,Info(jj).select_planes),3);
        basal_bg = nanmean(temp(:));
    end
    Info(jj).r_NAc_af_norm = nanmean(Info(jj).r_NAc_af(:,:,Info(jj).select_planes),3)/basal_bg;
    name_image = char(strcat(Info(jj).file(1:3) , '_' , Info(jj).prot, '.tif'));
    imwrite(Info(jj).r_NAc_af_norm/6,name_image,'tif')
    
end
Info_I1 = Info;
save('Info_I1.mat','Info_I1')