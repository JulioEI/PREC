%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Step 0: S0_registration.m : create Info_IX.mat grid structures from input z-stack images.

	input: folder path (Go to Registration folder) containing two folders (Basal, TBS (optostim condition)) with TIF z-stack images after registration using bUnwarpJ ImageJ plugin:
		-'XX_XX_c01_pos.tif': TIF file containing red channel aligned z-stack image
		-'XX_XX_NAc_pos.tif': TIF file containing binary NAc mask after alignment. Used for grid construction. 
		-'XX_XX_CORE_pos.tif': TIF file containing binary core NAc mask after alignment. Used for subregion delimitation. 
		-'XX_XX_SHELL_pos.tif': TIF file containing binary shell NAc mask after alignment. Used for subregion delimitation. 
	output: -'Info_IX.mat': matlab structure containing resulting 50umx50um grid data. User will be prompt to specify planes of interest (in S0_example, planes = 6 to 10)
				from original z-stack. 
		-'Image_IX.mat': matlab structutre containing original data.
	dependencies: 'registration_matlab.m'
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%WARNING: MANUALLY CHANGE VARIABLE NAME: Info_IX & Image_IX IN ALL THE
%SCRIPT (use control+F to replace)
loadf = uigetdir();
loadf = dir(loadf);
loadf(1:2) = [];

for hh = 1:size(loadf,1)
    my_folder = strcat(loadf(hh).folder, '\', loadf(hh).name);
    [Info_CS(hh), Image_PREC(hh)] = registration_matlab(my_folder);
end

Image_I1 = Image_PREC;
Info_I1 = Info_CS;

st_plane = input('Select starting plane for posterior analysis: ');
en_plane = input('Select ending plane for posterior analysis: ');
Info_I1(1).select_planes = st_plane:en_plane;
Info_I1(2).select_planes = st_plane:en_plane;

%save('Image_I1.mat', 'Image_I1','-v7.3');
%save('Info_I1.mat', 'Info_I1');
clear Image_PREC Info_CS