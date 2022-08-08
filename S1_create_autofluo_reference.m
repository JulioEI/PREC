%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Step 1: S1_create_autofluo_reference.m : code used to create a reference autoflorescence signal 
for each grid from control z-stack images with no viral infection.

	input:  -'Info_IX.mat': matlab structure containing control data (Info_I1, Info_I2, ..., Info_I9). Need
			to be loaded on the workspace.
	params: none
	output: -'Info_autofluo.mat': matlab structure containing selfdescriptive fields regarding autofluorescence signal. Info duplicated
			in two rows to ease later analysis (Basal & TBS). Key fields:
			-'r_bg': red background fluorescence outside NAc to account for microscope conditions
			-r_proj_og': mean red autofluorescence signal on each 50umx50um grid
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Create Info total BG
Info_autofluo = Info_I1;
Info_autofluo(1).file = 'autofluorescencia total';
Info_autofluo(2).file = 'autofluorescencia total';

%compute mean bg signal for each plane
Info_autofluo(1).r_bg = nan(10,9);
Info_autofluo(1).r_bg = [Info_I1(1).r_bg(:,1),Info_I2(1).r_bg(:,1),...
    Info_I3(1).r_bg(:,1),Info_I4(1).r_bg(:,1),Info_I5(1).r_bg(:,1),...
    Info_I6(1).r_bg(:,1),Info_I7(1).r_bg(:,1),Info_I8(1).r_bg(:,1),...
    Info_I9(1).r_bg(:,1)];

Info_autofluo(1).r_bg = mean(mean(Info_autofluo(1).r_bg,2,'omitnan'),'omitnan');
Info_autofluo(2).r_bg = Info_autofluo(1).r_bg;

%create prec pre-normalization
for ii = 1:2
    Info_I1(ii).r_core_og = Info_I1(ii).r_core;
    Info_I1(ii).r_shell_og = Info_I1(ii).r_shell;
    Info_I2(ii).r_core_og = Info_I2(ii).r_core;
    Info_I2(ii).r_shell_og = Info_I2(ii).r_shell;
    Info_I3(ii).r_core_og = Info_I3(ii).r_core;
    Info_I3(ii).r_shell_og = Info_I3(ii).r_shell;
    Info_I4(ii).r_core_og = Info_I4(ii).r_core;
    Info_I4(ii).r_shell_og = Info_I4(ii).r_shell;
    Info_I5(ii).r_core_og = Info_I5(ii).r_core;
    Info_I5(ii).r_shell_og = Info_I5(ii).r_shell;
    Info_I6(ii).r_core_og = Info_I6(ii).r_core;
    Info_I6(ii).r_shell_og = Info_I6(ii).r_shell;
    Info_I7(ii).r_core_og = Info_I7(ii).r_core;
    Info_I7(ii).r_shell_og = Info_I7(ii).r_shell;
    Info_I8(ii).r_core_og = Info_I8(ii).r_core;
    Info_I8(ii).r_shell_og = Info_I8(ii).r_shell;
    Info_I9(ii).r_core_og = Info_I9(ii).r_core;
    Info_I9(ii).r_shell_og = Info_I9(ii).r_shell;
end
%compute prec pre-normalization
for jj = 1:2
    for ii = 1:10
        Info_I1(jj).r_core_og(:,:,ii) = Info_I1(jj).r_core(:,:,ii).*Info_I1(jj).r_bg(ii,1);
        Info_I1(jj).r_shell_og(:,:,ii) = Info_I1(jj).r_shell(:,:,ii).*Info_I1(jj).r_bg(ii,1);
        Info_I2(jj).r_core_og(:,:,ii) = Info_I2(jj).r_core(:,:,ii).*Info_I2(jj).r_bg(ii,1);
        Info_I2(jj).r_shell_og(:,:,ii) = Info_I2(jj).r_shell(:,:,ii).*Info_I2(jj).r_bg(ii,1);
        Info_I3(jj).r_core_og(:,:,ii) = Info_I3(jj).r_core(:,:,ii).*Info_I3(jj).r_bg(ii,1);
        Info_I3(jj).r_shell_og(:,:,ii) = Info_I3(jj).r_shell(:,:,ii).*Info_I3(jj).r_bg(ii,1);
        Info_I4(jj).r_core_og(:,:,ii) = Info_I4(jj).r_core(:,:,ii).*Info_I4(jj).r_bg(ii,1);
        Info_I4(jj).r_shell_og(:,:,ii) = Info_I4(jj).r_shell(:,:,ii).*Info_I4(jj).r_bg(ii,1);
        Info_I5(jj).r_core_og(:,:,ii) = Info_I5(jj).r_core(:,:,ii).*Info_I5(jj).r_bg(ii,1);
        Info_I5(jj).r_shell_og(:,:,ii) = Info_I5(jj).r_shell(:,:,ii).*Info_I5(jj).r_bg(ii,1);
        Info_I6(jj).r_core_og(:,:,ii) = Info_I6(jj).r_core(:,:,ii).*Info_I6(jj).r_bg(ii,1);
        Info_I6(jj).r_shell_og(:,:,ii) = Info_I6(jj).r_shell(:,:,ii).*Info_I6(jj).r_bg(ii,1);
        Info_I7(jj).r_core_og(:,:,ii) = Info_I7(jj).r_core(:,:,ii).*Info_I7(jj).r_bg(ii,1);
        Info_I7(jj).r_shell_og(:,:,ii) = Info_I7(jj).r_shell(:,:,ii).*Info_I7(jj).r_bg(ii,1);
        Info_I8(jj).r_core_og(:,:,ii) = Info_I8(jj).r_core(:,:,ii).*Info_I8(jj).r_bg(ii,1);
        Info_I8(jj).r_shell_og(:,:,ii) = Info_I8(jj).r_shell(:,:,ii).*Info_I8(jj).r_bg(ii,1);
        Info_I9(jj).r_core_og(:,:,ii) = Info_I9(jj).r_core(:,:,ii).*Info_I9(jj).r_bg(ii,1);
        Info_I9(jj).r_shell_og(:,:,ii) = Info_I9(jj).r_shell(:,:,ii).*Info_I9(jj).r_bg(ii,1);
    end
end

for jj = 1:2
    Info_I1(jj).r_core_og_proj = mean(Info_I1(jj).r_core_og(:,:,Info_I1(jj).select_planes),3,'omitnan');
    Info_I1(jj).r_shell_og_proj = mean(Info_I1(jj).r_shell_og(:,:,Info_I1(jj).select_planes),3,'omitnan');
    Info_I2(jj).r_core_og_proj = mean(Info_I2(jj).r_core_og(:,:,Info_I2(jj).select_planes),3,'omitnan');
    Info_I2(jj).r_shell_og_proj = mean(Info_I2(jj).r_shell_og(:,:,Info_I2(jj).select_planes),3,'omitnan');
    Info_I3(jj).r_core_og_proj = mean(Info_I3(jj).r_core_og(:,:,Info_I3(jj).select_planes),3,'omitnan');
    Info_I3(jj).r_shell_og_proj = mean(Info_I3(jj).r_shell_og(:,:,Info_I3(jj).select_planes),3,'omitnan');
    Info_I4(jj).r_core_og_proj = mean(Info_I4(jj).r_core_og(:,:,Info_I4(jj).select_planes),3,'omitnan');
    Info_I4(jj).r_shell_og_proj = mean(Info_I4(jj).r_shell_og(:,:,Info_I4(jj).select_planes),3,'omitnan');
    Info_I5(jj).r_core_og_proj = mean(Info_I5(jj).r_core_og(:,:,Info_I5(jj).select_planes),3,'omitnan');
    Info_I5(jj).r_shell_og_proj = mean(Info_I5(jj).r_shell_og(:,:,Info_I5(jj).select_planes),3,'omitnan');
    Info_I6(jj).r_core_og_proj = mean(Info_I6(jj).r_core_og(:,:,Info_I6(jj).select_planes),3,'omitnan');
    Info_I6(jj).r_shell_og_proj = mean(Info_I6(jj).r_shell_og(:,:,Info_I6(jj).select_planes),3,'omitnan');
    Info_I7(jj).r_core_og_proj = mean(Info_I7(jj).r_core_og(:,:,Info_I7(jj).select_planes),3,'omitnan');
    Info_I7(jj).r_shell_og_proj = mean(Info_I7(jj).r_shell_og(:,:,Info_I7(jj).select_planes),3,'omitnan');
    Info_I8(jj).r_core_og_proj = mean(Info_I8(jj).r_core_og(:,:,Info_I8(jj).select_planes),3,'omitnan');
    Info_I8(jj).r_shell_og_proj = mean(Info_I8(jj).r_shell_og(:,:,Info_I8(jj).select_planes),3,'omitnan');
    Info_I9(jj).r_core_og_proj = mean(Info_I9(jj).r_core_og(:,:,Info_I9(jj).select_planes),3,'omitnan');
    Info_I9(jj).r_shell_og_proj = mean(Info_I9(jj).r_shell_og(:,:,Info_I9(jj).select_planes),3,'omitnan');
end
for jj = 1:2
    Info_I1(jj).r_core_proj = mean(Info_I1(jj).r_core(:,:,Info_I1(jj).select_planes),3,'omitnan');
    Info_I1(jj).r_shell_proj = mean(Info_I1(jj).r_shell(:,:,Info_I1(jj).select_planes),3,'omitnan');
    Info_I2(jj).r_core_proj = mean(Info_I2(jj).r_core(:,:,Info_I2(jj).select_planes),3,'omitnan');
    Info_I2(jj).r_shell_proj = mean(Info_I2(jj).r_shell(:,:,Info_I2(jj).select_planes),3,'omitnan');
    Info_I3(jj).r_core_proj = mean(Info_I3(jj).r_core(:,:,Info_I3(jj).select_planes),3,'omitnan');
    Info_I3(jj).r_shell_proj = mean(Info_I3(jj).r_shell(:,:,Info_I3(jj).select_planes),3,'omitnan');
    Info_I4(jj).r_core_proj = mean(Info_I4(jj).r_core(:,:,Info_I4(jj).select_planes),3,'omitnan');
    Info_I4(jj).r_shell_proj = mean(Info_I4(jj).r_shell(:,:,Info_I4(jj).select_planes),3,'omitnan');
    Info_I5(jj).r_core_proj = mean(Info_I5(jj).r_core(:,:,Info_I5(jj).select_planes),3,'omitnan');
    Info_I5(jj).r_shell_proj = mean(Info_I5(jj).r_shell(:,:,Info_I5(jj).select_planes),3,'omitnan');
    Info_I6(jj).r_core_proj = mean(Info_I6(jj).r_core(:,:,Info_I6(jj).select_planes),3,'omitnan');
    Info_I6(jj).r_shell_proj = mean(Info_I6(jj).r_shell(:,:,Info_I6(jj).select_planes),3,'omitnan');
    Info_I7(jj).r_core_proj = mean(Info_I7(jj).r_core(:,:,Info_I7(jj).select_planes),3,'omitnan');
    Info_I7(jj).r_shell_proj = mean(Info_I7(jj).r_shell(:,:,Info_I7(jj).select_planes),3,'omitnan');
    Info_I8(jj).r_core_proj = mean(Info_I8(jj).r_core(:,:,Info_I8(jj).select_planes),3,'omitnan');
    Info_I8(jj).r_shell_proj = mean(Info_I8(jj).r_shell(:,:,Info_I8(jj).select_planes),3,'omitnan');
    Info_I9(jj).r_core_proj = mean(Info_I9(jj).r_core(:,:,Info_I9(jj).select_planes),3,'omitnan');
    Info_I9(jj).r_shell_proj = mean(Info_I9(jj).r_shell(:,:,Info_I9(jj).select_planes),3,'omitnan');
end

for jj = 1:2
    %norm
    Info_autofluo(jj).r_proj = nan(size(Info_I1(jj).r_core_proj,1),size(Info_I1(jj).r_core_proj,2),size(Info_I1(jj).r_core_proj,3),9);
    temp = nan(size(Info_I1(1).r_core_proj,1),size(Info_I1(1).r_core_proj,2),9);
    temp(:,:,1) = Info_I1(1).r_core_proj;
    temp(:,:,2) = Info_I2(1).r_core_proj;
    temp(:,:,3) = Info_I3(1).r_core_proj;
    temp(:,:,4) = Info_I4(1).r_core_proj;
    temp(:,:,5) = Info_I5(1).r_core_proj;
    temp(:,:,6) = Info_I6(1).r_core_proj;
    temp(:,:,7) = Info_I7(1).r_core_proj;
    temp(:,:,8) = Info_I8(1).r_core_proj;
    temp(:,:,9) = Info_I9(1).r_core_proj;
    temp_core = temp;
    temp_core(isnan(temp_core)) = 0;
    temp = nan(size(Info_I1(1).r_shell_proj,1),size(Info_I1(1).r_shell_proj,2),9);
    temp(:,:,1) = Info_I1(1).r_shell_proj;
    temp(:,:,2) = Info_I2(1).r_shell_proj;
    temp(:,:,3) = Info_I3(1).r_shell_proj;
    temp(:,:,4) = Info_I4(1).r_shell_proj;
    temp(:,:,5) = Info_I5(1).r_shell_proj;
    temp(:,:,6) = Info_I6(1).r_shell_proj;
    temp(:,:,7) = Info_I7(1).r_shell_proj;
    temp(:,:,8) = Info_I8(1).r_shell_proj;
    temp(:,:,9) = Info_I9(1).r_shell_proj;
    temp_shell = temp;
    temp_shell(isnan(temp_shell)) = 0;
    Info_autofluo(jj).r_proj = mean(temp_shell + temp_core,3,'omitnan');

    %og
    Info_autofluo(jj).r_proj_og = nan(size(Info_I1(jj).r_core_og_proj,1),size(Info_I1(jj).r_core_og_proj,2),size(Info_I1(jj).r_core_og_proj,3),9);
    temp = nan(size(Info_I1(1).r_core_og_proj,1),size(Info_I1(1).r_core_og_proj,2),9);
    temp(:,:,1) = Info_I1(1).r_core_og_proj;
    temp(:,:,2) = Info_I2(1).r_core_og_proj;
    temp(:,:,3) = Info_I3(1).r_core_og_proj;
    temp(:,:,4) = Info_I4(1).r_core_og_proj;
    temp(:,:,5) = Info_I5(1).r_core_og_proj;
    temp(:,:,6) = Info_I6(1).r_core_og_proj;
    temp(:,:,7) = Info_I7(1).r_core_og_proj;
    temp(:,:,8) = Info_I8(1).r_core_og_proj;
    temp(:,:,9) = Info_I9(1).r_core_og_proj;
    temp_core = temp;
    temp_core(isnan(temp_core)) = 0;
    temp = nan(size(Info_I1(1).r_shell_og_proj,1),size(Info_I1(1).r_shell_og_proj,2),9);
    temp(:,:,1) = Info_I1(1).r_shell_og_proj;
    temp(:,:,2) = Info_I2(1).r_shell_og_proj;
    temp(:,:,3) = Info_I3(1).r_shell_og_proj;
    temp(:,:,4) = Info_I4(1).r_shell_og_proj;
    temp(:,:,5) = Info_I5(1).r_shell_og_proj;
    temp(:,:,6) = Info_I6(1).r_shell_og_proj;
    temp(:,:,7) = Info_I7(1).r_shell_og_proj;
    temp(:,:,8) = Info_I8(1).r_shell_og_proj;
    temp(:,:,9) = Info_I9(1).r_shell_og_proj;
    temp_shell = temp;
    temp_shell(isnan(temp_shell)) = 0;
    Info_autofluo(jj).r_proj_og = mean(temp_shell + temp_core,3,'omitnan');
end
%save everything
save('Info_I1.mat', 'Info_I1');
save('Info_I2.mat', 'Info_I2');
save('Info_I3.mat', 'Info_I3');
save('Info_I4.mat', 'Info_I4');
save('Info_I5.mat', 'Info_I5');
save('Info_I6.mat', 'Info_I6');
save('Info_I7.mat', 'Info_I7');
save('Info_I9.mat', 'Info_I9');
save('Info_autofluo.mat','Info_autofluo')