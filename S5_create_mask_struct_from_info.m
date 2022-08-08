%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Step 5:	'S5_create_mask_struct_from_info.m' : compute structure with containing mean data for a given experimental condition (e.g. amyg, mPFC). 
	input: 	-'Info_IX.mat': all matlab structures corresponding to a given experimental condition (see Step 2).
	params: -'threshold': threshold to determine active/inactive regions (default: 1.71 computed through a k-cluster approach, see Methods).
	output: -'struct_output.mat': structure output containing selfdescriptive fields with the results. Key fields include:
			-'proj_images': inner structure containing processed data for each individual pair of images of the condition. 
			-'IT_basal_norm': average matrix containing mean fluorescence signal for each grid on the experimental condition before optostimulation
					normalized to mean fluorescence basal condition.
			-'IT_tbs_norm': average matrix containing mean fluorescence signal for each grid on the experimental condition after optostimulation
					normalized to mean fluorescence basal condition.
			-'mask_norm': average binary matrix containing active regions inside NAc (computed according to threshold param).
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

threshold = 1.70543272569049;  
threshold_norm = 1.70543272569049; 

%threshold CaMPARI signals = 1.70543272569049
%threshold glutamatergic afferents = 1.17317683
%threshold dopaminergic afferents = 0.6722

z_scale = [0 6]; %z-axis limits for 3D plots 
c_scale = [0.8 2]; %colormap limits for plots
z_scale_norm = [0 6];
c_scale_norm = [0.8 2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
flag = 0;
%%
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

num_images = size(info,1);
fprintf('- %i image(s) have been detected: ',num_images);
fprintf('\n- The images are: ');
for ii = 1:size(info,1)
    fprintf('\n\t%i/%i: ''%s''', ii, num_images,info(ii).name);
end
%% Check planes
for ii = 1:size(info,1)
   info_check = eval(sprintf('%1$s', info(ii).name));
   if ~isfield(info_check,'select_planes')
       fprintf('Image %i does not contain the selected planes in Info. Finalizing code.', images(ii).num);
       flag = 1;
       return
   end
end
if flag==0
    group = input('\nWhat condition/area do these images belong to?: ','s');
    %%Create temp files
    proj_images = struct([]);
    IT_basal = [];
    IT_basal_core = [];
    IT_basal_shell = [];
    IT_tbs = [];
    IT_tbs_core = [];
    IT_tbs_shell = [];

    IT_basal_norm = [];
    IT_basal_core_norm = [];
    IT_basal_shell_norm = [];
    IT_tbs_norm = [];
    IT_tbs_core_norm = [];
    IT_tbs_shell_norm = [];

    %get core and shell masks
    core_mask = [];
    shell_mask = [];
    for ii = 1:size(info,1)
        info_temp = eval(sprintf('%1$s', info(ii).name));
        if isempty(core_mask)
            core_mask(:,:,end) = double(~isnan(info_temp(1).r_core_af(:,:,1)));
            core_mask(:,:,end+1) = double(~isnan(info_temp(2).r_core_af(:,:,1)));
            shell_mask(:,:,end) = double(~isnan(info_temp(1).r_shell_af(:,:,1)));
            shell_mask(:,:,end+1) = double(~isnan(info_temp(2).r_shell_af(:,:,1)));
        else
            core_mask(:,:,end+1) = double(~isnan(info_temp(1).r_core_af(:,:,1)));
            core_mask(:,:,end+1) = double(~isnan(info_temp(2).r_core_af(:,:,1)));
            shell_mask(:,:,end+1) = double(~isnan(info_temp(1).r_shell_af(:,:,1)));
            shell_mask(:,:,end+1) = double(~isnan(info_temp(2).r_shell_af(:,:,1)));
        end
    end
    core_mask = sum(core_mask,3);
    core_mask = double(core_mask~=0);
    shell_mask = sum(shell_mask,3);
    shell_mask = double(shell_mask>=size(info,1));
    core_mask(shell_mask==1) = 0;
    for ii = 1:size(info,1)
        fprintf('Creating temporal images according to planes in Info: %i/%i\n',ii, size(info,1))
        info_temp = eval(sprintf('%1$s', info(ii).name));

        %BASAL
        proj_images(2*ii-1).file = info_temp(1).file;
        proj_images(2*ii-1).select_planes = info_temp(1).select_planes;
        proj_images(2*ii-1).prot = info_temp(1).prot;
        proj_images(2*ii-1).group = group;
        proj_images(2*ii-1).num = info(ii).num;

        proj_images(2*ii-1).proj = mean(info_temp(1).r_core_af(:,:,info_temp(1).select_planes),3,'omitnan');
        proj_images(2*ii-1).proj(isnan(proj_images(2*ii-1).proj)) = 0;
        temp = mean(info_temp(1).r_shell_af(:,:,info_temp(1).select_planes),3,'omitnan');
        temp(isnan(temp)) = 0;

        proj_images(2*ii-1).proj = proj_images(2*ii-1).proj + temp;
        proj_images(2*ii-1).proj_core = proj_images(2*ii-1).proj.*core_mask;
        proj_images(2*ii-1).proj_shell = proj_images(2*ii-1).proj.*shell_mask;

        bg = proj_images(2*ii-1).proj;
        bg(bg==0) = NaN;
        bg = nanmean(bg(:));
        proj_images(2*ii-1).bg = bg;
        proj_images(2*ii-1).proj_norm = proj_images(2*ii-1).proj./bg;
        proj_images(2*ii-1).proj_core_norm = proj_images(2*ii-1).proj_core./bg;
        proj_images(2*ii-1).proj_shell_norm = proj_images(2*ii-1).proj_shell./bg;

        %TBS
        proj_images(2*ii).file = info_temp(2).file;
        proj_images(2*ii).select_planes = info_temp(2).select_planes;
        proj_images(2*ii).prot = info_temp(2).prot;
        proj_images(2*ii).group = group;
        proj_images(2*ii).num = info(ii).num;
        proj_images(2*ii).proj = mean(info_temp(2).r_core_af(:,:,info_temp(2).select_planes),3,'omitnan');
        proj_images(2*ii).proj(isnan(proj_images(2*ii).proj)) = 0;

        temp = mean(info_temp(2).r_shell_af(:,:,info_temp(2).select_planes),3,'omitnan');
        temp(isnan(temp)) = 0;
        proj_images(2*ii).proj = proj_images(2*ii).proj + temp;
        proj_images(2*ii).proj_core = proj_images(2*ii).proj.*core_mask;
        proj_images(2*ii).proj_shell = proj_images(2*ii).proj.*shell_mask;

        proj_images(2*ii).bg = bg;
        proj_images(2*ii).proj_norm = proj_images(2*ii).proj./bg;
        proj_images(2*ii).proj_core_norm = proj_images(2*ii).proj_core./bg;
        proj_images(2*ii).proj_shell_norm = proj_images(2*ii).proj_shell./bg;
        %To calculate the total later
            %first image
        if contains(lower(proj_images(2*ii-1).prot),'basal')
            if isempty(IT_basal)
                IT_basal(:,:,end) = proj_images(2*ii-1).proj;
                IT_basal_core(:,:,end) = proj_images(2*ii-1).proj_core;
                IT_basal_shell(:,:,end) = proj_images(2*ii-1).proj_shell;


                IT_basal_norm(:,:,end) = proj_images(2*ii-1).proj_norm;
                IT_basal_core_norm(:,:,end) = proj_images(2*ii-1).proj_core_norm;
                IT_basal_shell_norm(:,:,end) = proj_images(2*ii-1).proj_shell_norm;

            else
                IT_basal(:,:,end+1) = proj_images(2*ii-1).proj;
                IT_basal_core(:,:,end+1) = proj_images(2*ii-1).proj_core;
                IT_basal_shell(:,:,end+1) = proj_images(2*ii-1).proj_shell;

                IT_basal_norm(:,:,end+1) = proj_images(2*ii-1).proj_norm;
                IT_basal_core_norm(:,:,end+1) = proj_images(2*ii-1).proj_core_norm;
                IT_basal_shell_norm(:,:,end+1) = proj_images(2*ii-1).proj_shell_norm;
            end
        else
            if isempty(IT_tbs)
                IT_tbs(:,:,end) = proj_images(2*ii-1).proj;
                IT_tbs_core(:,:,end) = proj_images(2*ii-1).proj_core;
                IT_tbs_shell(:,:,end) = proj_images(2*ii-1).proj_shell;

                IT_tbs_norm(:,:,end) = proj_images(2*ii-1).proj_norm;
                IT_tbs_core_norm(:,:,end) = proj_images(2*ii-1).proj_core_norm;
                IT_tbs_shell_norm(:,:,end) = proj_images(2*ii-1).proj_shell_norm;
            else
                IT_tbs(:,:,end+1) = proj_images(2*ii-1).proj;
                IT_tbs_core(:,:,end+1) = proj_images(2*ii-1).proj_core;
                IT_tbs_shell(:,:,end+1) = proj_images(2*ii-1).proj_shell;

                IT_tbs_norm(:,:,end+1) = proj_images(2*ii-1).proj_norm;
                IT_tbs_core_norm(:,:,end+1) = proj_images(2*ii-1).proj_core_norm;
                IT_tbs_shell_norm(:,:,end+1) = proj_images(2*ii-1).proj_shell_norm;
            end
        end
            %second image
        if contains(lower(proj_images(2*ii).prot),'basal')
            if isempty(IT_basal)
                IT_basal(:,:,end) = proj_images(2*ii).proj;
                IT_basal_core(:,:,end) = proj_images(2*ii).proj_core;
                IT_basal_shell(:,:,end) = proj_images(2*ii).proj_shell;


                IT_basal_norm(:,:,end) = proj_images(2*ii).proj_norm;
                IT_basal_core_norm(:,:,end) = proj_images(2*ii).proj_core_norm;
                IT_basal_shell_norm(:,:,end) = proj_images(2*ii).proj_shell_norm;

            else
                IT_basal(:,:,end+1) = proj_images(2*ii).proj;
                IT_basal_core(:,:,end+1) = proj_images(2*ii).proj_core;
                IT_basal_shell(:,:,end+1) = proj_images(2*ii).proj_shell;

                IT_basal_norm(:,:,end+1) = proj_images(2*ii).proj_norm;
                IT_basal_core_norm(:,:,end+1) = proj_images(2*ii).proj_core_norm;
                IT_basal_shell_norm(:,:,end+1) = proj_images(2*ii).proj_shell_norm;
            end
        else
            if isempty(IT_tbs)
                IT_tbs(:,:,end) = proj_images(2*ii).proj;
                IT_tbs_core(:,:,end) = proj_images(2*ii).proj_core;
                IT_tbs_shell(:,:,end) = proj_images(2*ii).proj_shell;

                IT_tbs_norm(:,:,end) = proj_images(2*ii).proj_norm;
                IT_tbs_core_norm(:,:,end) = proj_images(2*ii).proj_core_norm;
                IT_tbs_shell_norm(:,:,end) = proj_images(2*ii).proj_shell_norm;
            else
                IT_tbs(:,:,end+1) = proj_images(2*ii).proj;
                IT_tbs_core(:,:,end+1) = proj_images(2*ii).proj_core;
                IT_tbs_shell(:,:,end+1) = proj_images(2*ii).proj_shell;

                IT_tbs_norm(:,:,end+1) = proj_images(2*ii).proj_norm;
                IT_tbs_core_norm(:,:,end+1) = proj_images(2*ii).proj_core_norm;
                IT_tbs_shell_norm(:,:,end+1) = proj_images(2*ii).proj_shell_norm;
            end
        end
    end
    %Plot all temp
    figure('units','normalized','outerposition',[0 0 1 1])
    tax = [];
    for ii = 1:size(proj_images,2)
        if rem(ii,2) == 1
            ax = subplot(2,size(proj_images,2)/2,(ii+1)/2);
        else
            ax = subplot(2,size(proj_images,2)/2,(ii/2)+ size(proj_images,2)/2);
        end
        im = proj_images(ii).proj;
        im(isnan(im)) = 0;
        surf(im,'FaceAlpha',1);
        set(gca, 'TickDir', 'out')
        zlim(z_scale)
        caxis(c_scale)

        colormap('parula');
        set(gca, 'XGrid', 'off')
        set(gca, 'YGrid', 'off')
        set(gca, 'ZGrid', 'off')
        set(gca, 'XTick', []);
        set(gca, 'YTick', []);
        set(gca, 'TitleFontSizeMultiplier',1.5)
        set(gca,'fontsize',14)
        title("I" + info(floor((ii+1)/2)).num + " -- "  + proj_images(ii).prot,'Color','k')
        zlabel('Photoconversion Index [a.u.]')
        set(gca, 'Color', 'k')
        set(gca, 'XColor', 'w')
        set(gca, 'YColor', 'w')
        set(gca, 'ZColor', 'w')
        axis off
        tax = [tax,ax];
    end
    Link = linkprop(tax,{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
    setappdata(gcf, 'StoreTheLink', Link);
    set(gcf,'Color', 'w')

        %Plot all temp norm
    figure('units','normalized','outerposition',[0 0 1 1])
    tax = [];
    for ii = 1:size(proj_images,2)
        if rem(ii,2) == 1
            ax = subplot(2,size(proj_images,2)/2,(ii+1)/2);
        else
            ax = subplot(2,size(proj_images,2)/2,(ii/2)+ size(proj_images,2)/2);
        end
        im = proj_images(ii).proj_norm;
        im(isnan(im)) = 0;
        surf(im,'FaceAlpha',1);
        set(gca, 'TickDir', 'out')
        zlim(z_scale_norm)
        caxis(c_scale_norm)

        colormap('parula');
        set(gca, 'XGrid', 'off')
        set(gca, 'YGrid', 'off')
        set(gca, 'ZGrid', 'off')
        set(gca, 'XTick', []);
        set(gca, 'YTick', []);
        set(gca, 'TitleFontSizeMultiplier',1.5)
        set(gca,'fontsize',14)
        title("I" + info(floor((ii+1)/2)).num + " Norm -- "  + proj_images(ii).prot,'Color','k')
        zlabel('Photoconversion Index [a.u.]')
        set(gca, 'Color', 'k')
        set(gca, 'XColor', 'w')
        set(gca, 'YColor', 'w')
        set(gca, 'ZColor', 'w')
        axis off
        tax = [tax,ax];
    end
    Link = linkprop(tax,{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
    setappdata(gcf, 'StoreTheLink', Link);
    set(gcf,'Color', 'w')
    %% Create total files
    IT_basal = squeeze(mean(IT_basal,3,'omitnan'));
    IT_basal_core = squeeze(mean(IT_basal_core,3,'omitnan'));
    IT_basal_shell = squeeze(mean(IT_basal_shell,3,'omitnan'));
    IT_tbs = squeeze(mean(IT_tbs,3,'omitnan'));
    IT_tbs_core = squeeze(mean(IT_tbs_core,3,'omitnan'));
    IT_tbs_shell = squeeze(mean(IT_tbs_shell,3,'omitnan'));

    IT_basal_norm = squeeze(mean(IT_basal_norm,3,'omitnan'));
    IT_basal_core_norm = squeeze(mean(IT_basal_core_norm,3,'omitnan'));
    IT_basal_shell_norm = squeeze(mean(IT_basal_shell_norm,3,'omitnan'));
    IT_tbs_norm = squeeze(mean(IT_tbs_norm,3,'omitnan'));
    IT_tbs_core_norm = squeeze(mean(IT_tbs_core_norm,3,'omitnan'));
    IT_tbs_shell_norm = squeeze(mean(IT_tbs_shell_norm,3,'omitnan'));

    %Plot total images norm
    figure
    ax1 = subplot(2,3,1);
    surf(IT_basal_norm, 'FaceAlpha', 1);
    shading interp
    set(gca, 'TickDir', 'out')
    title('IT Basal NORM')
    zlim(z_scale_norm)
    caxis(c_scale_norm)
    ax2 = subplot(2,3,4);
    surf(IT_tbs_norm, 'FaceAlpha', 1);
    shading interp
    set(gca, 'TickDir', 'out')
    title('IT TBS NORM')
    zlim(z_scale_norm)
    caxis(c_scale_norm)


    ax3 = subplot(2,3,2);
    surf(IT_basal_core_norm, 'FaceAlpha', 1);
    shading interp
    set(gca, 'TickDir', 'out')
    title('IT Basal core NORM')
    zlim(z_scale_norm)
    caxis(c_scale_norm)
    ax4 = subplot(2,3,5);
    surf(IT_tbs_core_norm, 'FaceAlpha', 1);
    shading interp
    set(gca, 'TickDir', 'out')
    title('IT TBS core NORM')
    zlim(z_scale_norm)
    caxis(c_scale_norm)

    ax5 = subplot(2,3,3);
    surf(IT_basal_shell_norm, 'FaceAlpha', 1);
    shading interp
    set(gca, 'TickDir', 'out')
    title('IT Basal shell NORM')
    zlim(z_scale_norm)
    caxis(c_scale_norm)
    ax6 = subplot(2,3,6);
    surf(IT_tbs_shell_norm, 'FaceAlpha', 1);
    shading interp
    set(gca, 'TickDir', 'out')
    title('IT TBS shell NORM')
    zlim(z_scale_norm)
    caxis(c_scale_norm)

    Link = linkprop([ax1, ax2,ax3,ax4,ax5,ax6],{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
    setappdata(gcf, 'StoreTheLink', Link);

    %% CREATE MASK
     %Create normal mask
    mask = double(IT_tbs >= threshold);
    %nmask = IT_tbs ~= 0; %original
    nmask = IT_tbs >  0.5; %modificada para que incluya en nmask solo valores > 0.5
    nmask = nmask - mask;
    nmask(nmask==0) = NaN;
    mask(mask==0) = NaN;

    mask_core = double(IT_tbs_core >= threshold);
    %nmask = IT_tbs ~= 0; %original
    nmask_core = IT_tbs_core >  0.5; %modificada para que incluya en nmask solo valores > 0.5
    nmask_core = nmask_core - mask_core;
    nmask_core(nmask_core==0) = NaN;
    mask_core(mask_core==0) = NaN;

    mask_shell = double(IT_tbs_shell >= threshold);
    %nmask = IT_tbs ~= 0; %original
    nmask_shell = IT_tbs_shell >  0.5; %modificada para que incluya en nmask solo valores > 0.5
    nmask_shell = nmask_shell - mask_shell;
    nmask_shell(nmask_shell==0) = NaN;
    mask_shell(mask_shell==0) = NaN;
    %% Create mask norm
    mask_norm = double(IT_tbs_norm>= threshold_norm);
    nmask_norm = IT_tbs_norm >  0.5; %modificada para que incluya en nmask solo valores > 0.5
    nmask_norm = nmask_norm - mask_norm;
    nmask_norm(nmask_norm==0) = NaN;
    mask_norm(mask_norm==0) = NaN;

    mask_core_norm = double(IT_tbs_core_norm>= threshold_norm);
    nmask_core_norm = IT_tbs_core_norm >  0.5; %modificada para que incluya en nmask solo valores > 0.5
    nmask_core_norm = nmask_core_norm - mask_core_norm;
    nmask_core_norm(nmask_core_norm==0) = NaN;
    mask_core_norm(mask_core_norm==0) = NaN;

    mask_shell_norm = double(IT_tbs_shell_norm>= threshold_norm);
    nmask_shell_norm = IT_tbs_shell_norm >  0.5; %modificada para que incluya en nmask solo valores > 0.5
    nmask_shell_norm = nmask_shell_norm - mask_shell_norm;
    nmask_shell_norm(nmask_shell_norm==0) = NaN;
    mask_shell_norm(mask_shell_norm==0) = NaN;

    %Plot mask norm
    figure
    subplot(2,3,1)
    imshow(mask_norm)
    title('mask norm')
    subplot(2,3,4)
    imshow(nmask_norm)
    title('Negative mask norm')

    subplot(2,3,2)
    imshow(mask_core_norm)
    title('mask core norm')
    subplot(2,3,5)
    imshow(nmask_core_norm)
    title('Negative mask core norm')

    subplot(2,3,3)
    imshow(mask_shell_norm)
    title('mask shell norm')
    subplot(2,3,6)
    imshow(nmask_shell_norm)
    title('Negative mask shell norm')
    %% Compute values inside outside of mask
    mask_values = zeros(size(proj_images,2)/2,4);
    for ii = 1:size(proj_images,2)
        temp = proj_images(ii).proj.*mask;
        proj_images(ii).mask_values = nanmean(temp(:));
        temp = proj_images(ii).proj.*nmask;
        proj_images(ii).nmask_values = nanmean(temp(:));
        if rem(ii,2) == 1
            mask_values((ii+1)/2,1) = proj_images(ii).mask_values;
            mask_values((ii+1)/2,2) = proj_images(ii).nmask_values;
        else
            mask_values(ii/2,3) = proj_images(ii).mask_values;
            mask_values(ii/2,4) = proj_images(ii).nmask_values;
        end
    end
    %% Compute values inside outside of mask CORE
    mask_values_core = zeros(size(proj_images,2)/2,4);
    for ii = 1:size(proj_images,2)
        temp = proj_images(ii).proj_core.*mask_core;
        proj_images(ii).mask_core_values = nanmean(temp(:));
        temp = proj_images(ii).proj_core.*nmask_core;
        proj_images(ii).nmask_core_values = nanmean(temp(:));
        if rem(ii,2) == 1
            mask_values_core((ii+1)/2,1) = proj_images(ii).mask_core_values;
            mask_values_core((ii+1)/2,2) = proj_images(ii).nmask_core_values;
        else
            mask_values_core(ii/2,3) = proj_images(ii).mask_core_values;
            mask_values_core(ii/2,4) = proj_images(ii).nmask_core_values;
        end
    end
    %% Compute values inside outside of mask SHELL
    mask_values_shell = zeros(size(proj_images,2)/2,4);
    for ii = 1:size(proj_images,2)
        temp = proj_images(ii).proj_shell.*mask_shell;
        proj_images(ii).mask_shell_values = nanmean(temp(:));
        temp = proj_images(ii).proj_shell.*nmask_shell;
        proj_images(ii).nmask_shell_values = nanmean(temp(:));
        if rem(ii,2) == 1
            mask_values_shell((ii+1)/2,1) = proj_images(ii).mask_shell_values;
            mask_values_shell((ii+1)/2,2) = proj_images(ii).nmask_shell_values;
        else
            mask_values_shell(ii/2,3) = proj_images(ii).mask_shell_values;
            mask_values_shell(ii/2,4) = proj_images(ii).nmask_shell_values;
        end
    end
    %% Compute values inside outside of mask NORM
    mask_values_norm = zeros(size(proj_images,2)/2,4);
    for ii = 1:size(proj_images,2)
        temp = proj_images(ii).proj_norm.*mask_norm;
        proj_images(ii).mask_values_norm = nanmean(temp(:));
        temp = proj_images(ii).proj_norm.*nmask_norm;
        proj_images(ii).nmask_values_norm = nanmean(temp(:));
        if rem(ii,2) == 1
            mask_values_norm((ii+1)/2,1) = proj_images(ii).mask_values_norm;
            mask_values_norm((ii+1)/2,2) = proj_images(ii).nmask_values_norm;
        else
            mask_values_norm(ii/2,3) = proj_images(ii).mask_values_norm;
            mask_values_norm(ii/2,4) = proj_images(ii).nmask_values_norm;
        end
    end
    %% Compute values inside outside of mask CORE NORM
    mask_values_core_norm = zeros(size(proj_images,2)/2,4);
    for ii = 1:size(proj_images,2)
        temp = proj_images(ii).proj_core_norm.*mask_core_norm;
        proj_images(ii).mask_core_values_norm = nanmean(temp(:));
        temp = proj_images(ii).proj_core_norm.*nmask_core_norm;
        proj_images(ii).nmask_core_values_norm = nanmean(temp(:));
        if rem(ii,2) == 1
            mask_values_core_norm((ii+1)/2,1) = proj_images(ii).mask_core_values_norm;
            mask_values_core_norm((ii+1)/2,2) = proj_images(ii).nmask_core_values_norm;
        else
            mask_values_core_norm(ii/2,3) = proj_images(ii).mask_core_values_norm;
            mask_values_core_norm(ii/2,4) = proj_images(ii).nmask_core_values_norm;
        end
    end
    %% Compute values inside outside of mask SHELL NORM
    mask_values_shell_norm = zeros(size(proj_images,2)/2,4);
    for ii = 1:size(proj_images,2)
        temp = proj_images(ii).proj_shell_norm.*mask_shell_norm;
        proj_images(ii).mask_shell_values_norm = nanmean(temp(:));
        temp = proj_images(ii).proj_shell_norm.*nmask_shell_norm;
        proj_images(ii).nmask_shell_values_norm = nanmean(temp(:));
        if rem(ii,2) == 1
            mask_values_shell_norm((ii+1)/2,1) = proj_images(ii).mask_shell_values_norm;
            mask_values_shell_norm((ii+1)/2,2) = proj_images(ii).nmask_shell_values_norm;
        else
            mask_values_shell_norm(ii/2,3) = proj_images(ii).mask_shell_values_norm;
            mask_values_shell_norm(ii/2,4) = proj_images(ii).nmask_shell_values_norm;
        end
    end
    %% Total value
    mask_values_norm_TOT = zeros(size(proj_images,2)/2,2);
    for ii = 1:size(proj_images,2)
        proj_images(ii).mask_values_norm_TOT = nanmean(proj_images(ii).proj_norm(:));
        if rem(ii,2) == 1
            mask_values_norm_TOT((ii+1)/2,1) = proj_images(ii).mask_values_norm_TOT;
        else
            mask_values_norm_TOT(ii/2,2) = proj_images(ii).mask_values_norm_TOT;
        end
    end
    %% Total value CORE
    mask_values_core_norm_TOT = zeros(size(proj_images,2)/2,2);
    for ii = 1:size(proj_images,2)
        proj_images(ii).mask_core_values_norm_TOT = nanmean(proj_images(ii).proj_core_norm(:));
        if rem(ii,2) == 1
            mask_values_core_norm_TOT((ii+1)/2,1) = proj_images(ii).mask_core_values_norm_TOT;
        else
            mask_values_core_norm_TOT(ii/2,2) = proj_images(ii).mask_core_values_norm_TOT;
        end
    end
    %% Total value Shell
    mask_values_shell_norm_TOT = zeros(size(proj_images,2)/2,2);
    for ii = 1:size(proj_images,2)
        proj_images(ii).mask_shell_values_norm_TOT = nanmean(proj_images(ii).proj_shell_norm(:));
        if rem(ii,2) == 1
            mask_values_shell_norm_TOT((ii+1)/2,1) = proj_images(ii).mask_shell_values_norm_TOT;
        else
            mask_values_shell_norm_TOT(ii/2,2) = proj_images(ii).mask_shell_values_norm_TOT;
        end
    end
end
clear ax ax1 ax2 ax3 ax4 ax5 ax6 ii image_temp info_check pos s1 s2 tax temp info_temp
struct_output.mask_values = mask_values;
struct_output.mask_values_core = mask_values_core;
struct_output.mask_values_shell = mask_values_shell;

struct_output.mask_values_norm = mask_values_norm;
struct_output.mask_values_core_norm = mask_values_core_norm;
struct_output.mask_values_shell_norm = mask_values_shell_norm;

struct_output.mask_values_norm_TOT = mask_values_norm_TOT;
struct_output.mask_values_core_norm_TOT = mask_values_core_norm_TOT;
struct_output.mask_values_shell_norm_TOT = mask_values_shell_norm_TOT;

struct_output.proj_images = proj_images;
struct_output.core_delimit = core_mask;
struct_output.shell_delimit = shell_mask;

struct_output.threshold = threshold;
struct_output.threshold_norm = threshold_norm;
struct_output.mask = mask;
struct_output.mask_core = mask_core;
struct_output.mask_shell = mask_shell;
struct_output.nmask = nmask;
struct_output.nmask_core = nmask_core;
struct_output.nmask_shell = nmask_shell;

struct_output.IT_basal = IT_basal;
struct_output.IT_basal_core = IT_basal_core;
struct_output.IT_basal_shell = IT_basal_shell;
struct_output.IT_tbs = IT_tbs;
struct_output.IT_tbs_core = IT_tbs_core;
struct_output.IT_tbs_shell = IT_tbs_shell;

struct_output.mask_norm = mask_norm;
struct_output.mask_core_norm = mask_core_norm;
struct_output.mask_shell_norm = mask_shell_norm;
struct_output.nmask_norm = nmask_norm;
struct_output.nmask_core_norm = nmask_core_norm;
struct_output.nmask_shell_norm = nmask_shell_norm;


struct_output.IT_basal_norm = IT_basal_norm;
struct_output.IT_basal_core_norm = IT_basal_core_norm;
struct_output.IT_basal_shell_norm = IT_basal_shell_norm;
struct_output.IT_tbs_norm = IT_tbs_norm;
struct_output.IT_tbs_core_norm = IT_tbs_core_norm;
struct_output.IT_tbs_shell_norm = IT_tbs_shell_norm;

