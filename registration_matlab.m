function [Info_CS, Image_PREC] = registration_matlab(my_folder)
clearw = warning ('off','all');
folder = dir(my_folder);
del = [];
for ii = 1:size(folder,1) 
    if folder(ii).isdir == 1
        del(end+1) = ii;
    elseif ~contains(folder(ii).name, '_pos.tif')
        del(end+1) = ii;
    elseif contains(folder(ii).name, '_ch00')
        del(end+1) = ii;
    end
end
folder(del) = [];
if contains(lower(my_folder(end-5:end)),"tbs")
    cond = "TBS";
elseif contains(lower(my_folder(end-5:end)),"basal")
    cond = "Basal";
else
    cond = "Unkown";
end
if size(folder,1) ~= 5
    keyboard;
end
tic
%% LOAD TIFF AND GET INFO
InfoImage_NA=imfinfo(strcat(folder(5).folder,'\',folder(5).name));
fprintf('\nLoading File...       '); 
%Background ROI
TifLink=Tiff(strcat(folder(1).folder,'\',folder(1).name),'r');    
TifLink.setDirectory(1);
temp = double(TifLink.read());
Mask_bg = temp(:,:,1);
Mask_bg(Mask_bg~=0) = 1;
Mask_bg(Mask_bg==0) = NaN;
    
%Core mask
TifLink=Tiff(strcat(folder(2).folder,'\',folder(2).name),'r');    
TifLink.setDirectory(1);
temp = double(TifLink.read());
Mask_core = temp(:,:,1);
Mask_core(Mask_core~=0) = 1;
Mask_core(Mask_core==0) = NaN;
    
%Shell mask
TifLink=Tiff(strcat(folder(4).folder,'\',folder(4).name),'r');    
TifLink.setDirectory(1);
temp = double(TifLink.read());
Mask_shell = temp(:,:,1);
Mask_shell(Mask_shell~=0) = 1;
Mask_shell(Mask_shell==0) = NaN;

    
    
Info_CS.file = folder(5).name;
Info_CS.my=InfoImage_NA(1,1).Width;
Info_CS.mx=InfoImage_NA(1,1).Height;
Info_CS.xresolution = 0.75758; %um/pixel
Info_CS.yresolution = 0.75758; %um/pixel
Info_CS.zresolution = 9.99; %um/pixel
Info_CS.planes= numel(InfoImage_NA);
Info_CS.MaxValue = InfoImage_NA(1,1).MaxSampleValue(1);
    
s = warning;
warning off all;
%The tiff loading is done according to the instructions at http://www.matlabtips.com/how-to-load-tiff-stacks-fast-really-fast/
Image_r = nan(Info_CS.mx,Info_CS.my,Info_CS.planes);
        
%get red
TifLink=Tiff(strcat(folder(5).folder,'\',folder(5).name),'r');
for ii=1:Info_CS.planes
	TifLink.setDirectory(ii);
	temp = double(TifLink.read());
	Image_r(:,:,ii) = temp;
end
if Info_CS.MaxValue == 255
	Image_r = Image_r*(65535/255);
end
    
    %% Check size
    if size(Mask_core,1) < size(Image_r,1)
        Mask_core = [nan(size(Image_r,1) - size(Mask_core,1),size(Mask_core,2)); Mask_core];
    elseif size(Mask_core,1) > size(Image_r,1)
	    Mask_core = Mask_core(1:size(Image_r,1),:);
    end
	if size(Mask_core,2) < size(Image_r,2)
        Mask_core = [nan(size(Mask_core,1), size(Image_g,2) - size(Mask_core,2)), Mask_core];
	elseif size(Mask_core,2) > size(Image_r,2)
        Mask_core = Mask_core(:, 1:size(Image_r,2));
    end
	if size(Mask_shell,1) < size(Image_r,1)
        Mask_shell = [nan(size(Image_r,1) - size(Mask_shell,1),size(Mask_shell,2)); Mask_shell];
	elseif size(Mask_shell,1) > size(Image_r,1)
        Mask_shell = Mask_shell(1:size(Image_r,1),:);
    end
	if size(Mask_shell,2) < size(Image_r,2)
        Mask_shell = [nan(size(Mask_shell,1), size(Image_r,2) - size(Mask_shell,2)), Mask_shell];
	elseif size(Mask_shell,2) > size(Image_r,2)
        Mask_shell = Mask_shell(:, 1:size(Image_r,2));
    end
	if size(Mask_bg,1) < size(Image_r,1)
        Mask_bg = [nan(size(Image_r,1) - size(Mask_bg,1),size(Mask_bg,2)); Mask_bg];
	elseif size(Mask_bg,1) > size(Image_r,1)
        Mask_bg = Mask_bg(1:size(Image_r,1),:);
    end
	if size(Mask_bg,2) < size(Image_r,2)
        Mask_bg = [nan(size(Mask_bg,1), size(Image_r,2) - size(Mask_bg,2)), Mask_bg];
	elseif size(Mask_bg,2) > size(Image_r,2)
        Mask_bg = Mask_bg(:, 1:size(Image_r,2));
    end
    Mask_NA = Mask_core;
    Mask_NA(isnan(Mask_NA)) = 0;
    Mask_NA2 = Mask_shell;
    Mask_NA2(isnan(Mask_NA2)) = 0;
    Mask_NA = Mask_NA + Mask_NA2;
    Mask_NA(Mask_NA==0) = NaN;
    Mask_NA(~isnan(Mask_NA)) = 1;
    clear Mask_NA2;
    
    Irc = (Image_r.*Mask_core);
    Irs = (Image_r.*Mask_shell);
    warning(s);
    TifLink.close();
    Info_CS.r_bg = squeeze(sum(sum(Image_r.*Mask_bg,1,'omitnan'),2,'omitnan')./(sum(sum(Mask_bg,1,'omitnan'),2,'omitnan')));
    for ii = 1:size(Irc,3)
        Irc(:,:,ii) = Irc(:,:,ii)./Info_CS.r_bg(ii);
        Irs(:,:,ii) = Irs(:,:,ii)./Info_CS.r_bg(ii);
    end
    Irc(Irc<0) = 0.0001;
    Irs(Irs<0) = 0.0001;

    Info_CS.r_core = squeeze(sum(sum(Irc,1,'omitnan'),2,'omitnan')./(sum(sum(Mask_core,1,'omitnan'),2,'omitnan')));
    Info_CS.r_shell = squeeze(sum(sum(Irs,1,'omitnan'),2,'omitnan')./(sum(sum(Mask_shell,1,'omitnan'),2,'omitnan')));
    
    
    Info_CS.area_NAc = (sum(sum(Mask_core, 'omitnan'),'omitnan'))*Info_CS.xresolution*Info_CS.yresolution*(10^-6); %mm^2
    Info_CS.area_NAs = (sum(sum(Mask_shell, 'omitnan'),'omitnan'))*Info_CS.xresolution*Info_CS.yresolution*(10^-6); %mm^2
    Info_CS.area_NA = (sum(sum(Mask_NA, 'omitnan'),'omitnan'))*Info_CS.xresolution*Info_CS.yresolution*(10^-6); %mm^2
    Info_CS.prot = cond;
    Image_PREC.NAc_mask = Mask_core;
    Image_PREC.NAs_mask = Mask_shell;
    Image_PREC.NA_mask = Mask_NA;
    Image_PREC.Image_r = Image_r;
    fprintf('\b\b\b\b\bDone\n');
    %% SUBDIVIDE IN SQUARES AND GET RED/GREEN FOR EACH
    fprintf('Dividing image in squares...       ');

    PREC_core_r = nan(size(Irc));
    PREC_shell_r = nan(size(Irs));
    
    wx = round(50/Info_CS.xresolution);
    wy = round(50/Info_CS.yresolution);
    
    nwx = floor(Info_CS.mx/wx);
    nwy = floor(Info_CS.my/wy);
    Info_CS.wx = wx;
    Info_CS.wy = wy;
    Info_CS.nwx = nwx;
    Info_CS.nwy = nwy;
        %core
    r_core = nan(nwx, nwy,size(PREC_core_r,3));
    for xx=1:nwx
        for yy = 1:nwy
            wind_red = Irc(1+(xx-1)*wx:xx*wx,1+(yy-1)*wy:yy*wy,:);
            temp = unique(wind_red);
            temp(isnan(temp)) = [];
            if size(temp,1)>1              
                red_temp =   squeeze(sum(sum(wind_red,1,'omitnan'),2,'omitnan')./(size(wind_red,1)*size(wind_red,2)));
                r_core(xx,yy,:) = red_temp;
                PREC_core_r(1+(xx-1)*wx:xx*wx,1+(yy-1)*wy:yy*wy,:) = Mask_core(1+(xx-1)*wx:xx*wx,1+(yy-1)*wy:yy*wy).*r_core(xx,yy,:);
          
            end
        end
    end
    Info_CS.r_core = r_core;
    Image_PREC.PREC_core_r = PREC_core_r;
        %shell
    r_shell = nan(nwx, nwy,size(PREC_shell_r,3));
    for xx=1:nwx
        for yy = 1:nwy
            wind_red = Irs(1+(xx-1)*wx:xx*wx,1+(yy-1)*wy:yy*wy,:);
            temp = unique(wind_red);
            temp(isnan(temp)) = [];
            if size(temp,1)>1
                red_temp = squeeze(sum(sum(wind_red,1,'omitnan'),2,'omitnan')./(size(wind_red,1)*size(wind_red,2)));
                r_shell(xx,yy,:) = red_temp;
                PREC_shell_r(1+(xx-1)*wx:xx*wx,1+(yy-1)*wy:yy*wy,:) = Mask_shell(1+(xx-1)*wx:xx*wx,1+(yy-1)*wy:yy*wy).*r_shell(xx,yy,:);
            end
        end
    end
    Info_CS.r_shell = r_shell;
    Image_PREC.PREC_shell_r = PREC_shell_r;
    fprintf('\b\b\b\b\bDone\n');
    toc
    
end