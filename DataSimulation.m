% This code simulate gradient echo (GRE) data: Magnitude and phase ( total, positive, and negative) (see README.md
% for more details). 

%%%%%%%%%%%%%%%%%%%%% VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%
%%% Inputs: %%% 
% Chimap_file- 3D total susceptibility map
% Chiposmap_file- 3D positive susceptibility map
% Chinegmap_file- 3D negative susceptibility map
% R1map_file- 3D R1 map (1/T1).
% R2starmap_file- 3D R2star map (1/T2star).
% M0map_file- 3D Net magnetization map (M0).             
% Segmentation_file- 3D brain segmentation file.         
% BrainMask_file- 3D maks of ROI.
% R2map_file- 3D R2star map (1/T2).
% Drposmap_file- 3D positive Dr.
% Drnegmap_file- 3D negative Dr.
% TR- Repition Time in seconds.
% TE- Echo Time in seconds (array).
% FlipAngle- Flip angle in degrees above ernst angle (13) to improve gre contrast between Grey and white matter.
% B0- Main magnetic field strength in Tesla.
% B0_dir- B0 direction .
% PhaseOffset- Phase Offset (multiplier term of a quadratic phase over the brain)
% Shim- boolean 0 if no additional shimms are applied 1 if
% Res- resolution of output
%%% Output: %%%
%
% Output_dir- Path to results ( 4D Simulated GRE magnitude and phase)
%%%%%%%%%%%%%%%%%%%%% VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%
%% T2 simulation
seg_file = 'data/masks/SegmentedModel.nii.gz';
T2_star_file = 'data/maps/T2_star.nii.gz';
M0_file = 'data/maps/M0.nii.gz';
% Output:
T2_uniform_file = 'data/maps/T2_uniform.nii.gz';
R2_file = 'data/maps/R2.nii.gz';

% Load the NIfTI segmentation file
seg_img = load_untouch_nii(seg_file);
seg_data = double(seg_img.img);

% Define the region values lookup table
region_values = [57.46, 41.47, 50.44, 44.07, 71.71, 47.255, 56.62, 45.54, 84.71, 1029.6, 97.5];

% Assign values to each region in the segmentation
for region_num = 1:11
    region_voxels = (seg_data == region_num);
    region_value = region_values(region_num);
    seg_data(region_voxels) = region_value;
end

% Create a new NIfTI image with the modified segmentation data
modified_seg_img = make_nii(seg_data, seg_img.hdr.dime.pixdim(2:4), [], 64);

% Save the modified NIfTI segmentation file
save_nii(modified_seg_img, T2_uniform_file);

% Load NIfTI files
T2_star_img = load_untouch_nii(T2_star_file);
T2_img = load_untouch_nii(T2_uniform_file);
M0_img = load_untouch_nii(M0_file);
segmentation_img=load_untouch_nii(seg_file);
% Get the data arrays
T2_star_data = double(T2_star_img.img);
T2_data = double(T2_img.img);
M0_data = double(M0_img.img);
segmentation_data=double(segmentation_img.img);
% Preprocess T2_star and M0 maps
T2_star_data(isnan(T2_star_data)) = 0;
M0_data(isnan(M0_data)) = 0;

% Create empty arrays to store percentage masks for T2_star and M0 variations
percentage_masks_t2 = zeros(size(T2_star_data));
percentage_masks_m0 = zeros(size(M0_data));

% Loop over each region (from 1 to 11) to calculate variations
for region_id = 1:11
    % Create a mask for the current region based on the segmentation map
    region_mask = (segmentation_data == region_id);

    % Calculate variations for T2_star
    percentile_value_in_region_t2 = mean(T2_star_data(region_mask));
    percentage_mask_t2 = region_mask .* (T2_star_data / percentile_value_in_region_t2);
    percentage_masks_t2 = percentage_masks_t2 + percentage_mask_t2;

    % Calculate variations for M0
    percentile_value_in_region_m0 = mean(M0_data(region_mask));
    percentage_mask_m0 = region_mask .* (M0_data / percentile_value_in_region_m0);
    percentage_masks_m0 = percentage_masks_m0 + percentage_mask_m0;
end

% Combine the variations from T2_star and M0 into a single percentage mask
percentage_masks_combined = (percentage_masks_t2 + percentage_masks_m0) / 2.0;

% Apply the combined percentage masks to the T2 map to get T2_with_variations
T2_with_variations = T2_data .* percentage_masks_combined;

% Add a 0.2 Gaussian filter to T2_with_variations
T2_with_variations_smoothed = imgaussfilt(T2_with_variations, 0.2);
R2=(1/T2_with_variations_smoothed)*1000;
R2(isnan(R2)) = 0;
R2(isinf(R2)) = 0;
R2(R2>100)=0;
% Save the T2_with_variations_smoothed map as a new NIfTI file
T2_with_variations_smoothed_img = make_nii(R2, T2_img.hdr.dime.pixdim(2:4), [], 64);
save_nii(T2_with_variations_smoothed_img, R2_file);
Mask(R2_file)
%% Creating Dr+ and Dr-
input_nifti_path = 'data/masks/SegmentedModel.nii.gz';
angle_nifti_path = 'data/maps/theta.nii.gz';
% Output:
Dr_pos_file = 'data/maps//Dr_positive.nii.gz';
Dr_neg_variable_file = 'data/maps/Dr_negative_variable.nii.gz';
Dr_neg_constant_file = 'data/maps/Dr_negative_constant.nii.gz';

% Load the segmentation NIfTI image
image = load_untouch_nii(input_nifti_path);
data = double(image.img);

% Define the regions to modify (1 to 11, excluding 8)
regions_to_modify = [1, 2, 3, 4, 5, 6, 7, 9, 10, 11];

% Set the specified regions to the desired value for Dr_positive
value_positive = (2 * pi * 2 * pi * 42.58 * 7) / (9 * sqrt(3));
for region = regions_to_modify
    data(data == region) = value_positive;
end
data(data == 8) = 0;

% Create and save a new NIfTI image with the modified data for Dr_positive
modified_image_pos = make_nii(data, image.hdr.dime.pixdim(2:4), [], 64);
save_nii(modified_image_pos, Dr_pos_file);
Mask(Dr_pos_file)

% Now for Dr negative variable
% Load the angle NIfTI image 
angle_image = load_untouch_nii(angle_nifti_path);
angle_data = double(angle_image.img);

% Convert angle from degrees to radians
angle_data_radians = deg2rad(angle_data);

% Calculate Dr_negative using the formula you provided
Drneg_data = 0.5 * 42.58 * 2 * pi * 7 * sin(angle_data_radians).^2;
Drneg_data(isnan(Drneg_data)) = 0;
% Create a new NIfTI image for Dr_negative
Drneg_image = make_nii(Drneg_data, angle_image.hdr.dime.pixdim(2:4), [], 64);

% Save the Dr_negative NIfTI image
save_nii(Drneg_image, Dr_neg_variable_file);
Mask(Dr_neg_variable_file);

% Now for Dr negative constant
% Load the segmentation NIfTI image
image = load_untouch_nii(input_nifti_path);
data2 = double(image.img);
% Create a new NIfTI image with a constant value (700.8) for Dr_negative constant
value_constant = 700.8;
% Set the specified regions to the desired value for Dr_positive
for region = regions_to_modify
    data2(data2 == region) = 0;
end
data2(data2 == 8) = value_constant;
% Create a new NIfTI image for Dr_negative constant
Drneg_image_constant = make_nii(data2, angle_image.hdr.dime.pixdim(2:4), [], 64);

% Save the Dr_negative constant NIfTI image
save_nii(Drneg_image_constant, Dr_neg_constant_file);
Mask(Dr_neg_constant_file);

%% Sequence protocol for simulation

ModelParams.Chimap_file = 'data/chimodel/Chi.nii.gz';
ModelParams.Chiposmap_file='data/chimodel/Chi_positive.nii.gz';
ModelParams.Chinegmap_file='data/chimodel/Chi_negative.nii.gz';
ModelParams.R1map_file = 'data/maps/R1.nii.gz';
ModelParams.R2starmap_file = 'data/maps/R2star.nii.gz';
ModelParams.M0map_file = 'data/maps/M0.nii.gz';
ModelParams.Segmentation_file = 'data/masks/SegmentedModel.nii.gz';
ModelParams.BrainMask_file = 'data/masks/BrainMask.nii.gz';
ModelParams.R2map_file='data/maps/R2.nii.gz';
ModelParams.Drposmap_file='data/maps/Dr_positive.nii.gz';
ModelParams.Drnegmap_file='data/maps/Dr_negative_variable.nii.gz';

Protocol = 3 ;
SeqParams{Protocol}.TR = 50e-3;                               
SeqParams{Protocol}.TE = [4e-3 12e-3  20e-3  28e-3  ];     
SeqParams{Protocol}.FlipAngle = 15;                             
SimParams{Protocol}.B0 = 7;
SimParams{Protocol}.B0_dir = [0 0 1];
SimParams{Protocol}.PhaseOffset = 0 ;                           
% 0 no phase offset; pi phase difference inside brain mask
SimParams{Protocol}.Shimm = 1;                                  
% extracted output
SimParams{Protocol}.Res = [1 1 1];                              
SimParams{Protocol}.Output_dir = 'Simdata/Simulation_Results'; 

DataSimulationFunction(ModelParams,SeqParams{Protocol},SimParams{Protocol})




