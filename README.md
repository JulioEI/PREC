# PRQ

**Step 0:** S0_registration.m : create Info_IX.mat grid structures from input z-stack images.

	input: folder path (Go to Registration folder) containing two folders (Basal, TBS (optostim condition)) with TIF z-stack images 
	after registration using bUnwarpJ ImageJ plugin:
		-'XX_XX_c01_pos.tif': TIF file containing red channel aligned z-stack image
		-'XX_XX_NAc_pos.tif': TIF file containing binary NAc mask after alignment. Used for grid construction. 
		-'XX_XX_CORE_pos.tif': TIF file containing binary core NAc mask after alignment. Used for subregion delimitation. 
		-'XX_XX_SHELL_pos.tif': TIF file containing binary shell NAc mask after alignment. Used for subregion delimitation. 
	output: -'Info_IX.mat': matlab structure containing resulting 50umx50um grid data. User will be prompt to specify planes of
	interest (in S0_example, planes = 6 to 10)
				from original z-stack. 
		-'Image_IX.mat': matlab structutre containing original data (not included on the example output for file-size reasons, not 
		necessary for further processing).
	dependencies: 'registration_matlab.m'


**Step 1:** S1_create_autofluo_reference.m : code used to create a reference autoflorescence signal 
for each grid from control z-stack images with no viral infection.

	input:  -'Info_IX.mat': matlab structure containing control data (Info_I1, Info_I2, ..., Info_I9). Need
			to be loaded on the workspace.
	params: none
	output: -'Info_autofluo.mat': matlab structure containing selfdescriptive fields regarding autofluorescence signal. Info duplicated
			in two rows to ease later analysis (Basal & TBS). Key fields:
			-'r_bg': red background fluorescence outside NAc to account for microscope conditions
			-r_proj_og': mean red autofluorescence signal on each 50umx50um grid


**Step 2:** S2_fix_info_files.m: substract reference autofluorescence signal from experimental Info_IX structures.

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


**Step 3 (optional):** S3_plot_change_in_Info_files.m :code used to visualize change in images after substracting the autofluorescence signal
	input: 	-'Info_IX.mat' matlab structure containing experimental processed data (after Step 2)
		-'Info_autofluo.mat' (see Step 1). 
	params: -'orig_color_th': upper threshold colormap for original data
		-'norm_color_th': upper threshold colormap for normalized data
	output: figure for data visualization

**Step 4:** 'S4_k_cluster_threshold.m' : compute activation threshold through a k-mean clustering approach.
	input: -'Info_IX.mat': all matlab structures used to compute an activation threshold.
	output: -'k_XX': resulting clustering groups
			'k_mean': mean fluorescence signal of each clustering group
			'k_min': min fluorescence signal of each clustering group
		-'Ik': final image where each pixel values corresponds to that of the cluster it has been assigned.
	dependencies: 'k_cluster_seg.m'


**Step 5:** 'S5_create_mask_struct_from_info.m' : compute structure with containing mean data for a given experimental condition (e.g. amyg, mPFC). 
	input: 	-'Info_IX.mat': all matlab structures corresponding to a given experimental condition (see Step 2).
	params: -'threshold': threshold to determine active/inactive regions (default: 1.71 computed through a k-cluster approach, see Step4 
	and Methods for further information).
	output: -'struct_output.mat': structure output containing selfdescriptive fields with the results. Key fields include:
			-'proj_images': inner structure containing processed data for each individual pair of images of the condition. 
			-'IT_basal_norm': average matrix containing mean fluorescence signal for each grid on the experimental condition before
			optostimulation normalized to mean fluorescence basal condition.
	 		-'IT_tbs_norm': average matrix containing mean fluorescence signal for each grid on the experimental condition after 
			optostimulation normalized to mean fluorescence basal condition.
			-'mask_norm': average binary matrix containing active regions inside NAc (computed according to threshold param).


Note that this pipeline allows for data exploration in terms of stimuli condition, NAc subregions, and data normalization. In our case, we 
focused on normalized data considering the global NAc.
