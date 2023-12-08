% This code Create three susceptibility phantoms ( total, positive, and negative) (see README.md
% for more details). 

%%%%%%%%%%%%%%%%%%%%% VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%
%%% Inputs: %%% 
%
% R1map_file- 3D R1 map (1/T1).
% R2starmap_file- 3D R2star map (1/T2star).
% M0map_file- 3D Net magnetization map (M0).             
% Segmentation_file- 3D brain segmentation file.
% rawField_file- 3D raw field map.           
% BrainMask_file- 3D maks of ROI.
% highGradMask_file- 3D high gradient mask
% ChiModulation_file- Mat file contains all susceptibility and weighting values.
%%% Output: %%%
%
% OutputChiModel_file- Chi: total, positive, and negative ( Chi total and negative are created in CreateOwnRealisticPhantom file)
%%%%%%%%%%%%%%%%%%%%% VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%
%% parameters

ModelParams.R1map_file = 'data/maps/R1.nii.gz';
ModelParams.R2starmap_file = 'data/maps/R2star.nii.gz';
ModelParams.M0map_file = 'data/maps/M0.nii.gz';
ModelParams.Segmentation_file = 'data/masks/SegmentedModel.nii.gz';
ModelParams.rawField_file = 'data/raw/MP2RAGEME_total-field_ppm_head_r2s-header.nii.gz';
ModelParams.BrainMask_file = 'data/masks/BrainMask.nii.gz';
% file that defines regions where R2* values can and can't be trusted
ModelParams.highGradMask_file = 'data/masks/highgrad.nii.gz';
% File with the various parameters to create a susceptbility map
ModelParams.ChiModulation_file = 'data/chimodel/SusceptibilityValues.mat';
ModelParams.OutputChiModel_file = 'data/chimodel/Chi_positive.nii.gz';
PhantomCreationFunction(ModelParams)
%Mask('data/chimodel/Chi_positive_masked.nii.gz');
%Mask('data/chimodel/Chi.nii.gz_masked');
%Mask('data/chimodel/Chi_negative_masked.nii.gz');
%% Add anisotropy
%% Theta
v1_img_path = 'data/raw/3T/Derived_HighResSpace/V1.nii.gz';
segmentation_path = 'data/masks/SegmentedModel.nii.gz';

v1_img = load_untouch_nii(v1_img_path);
segmentation_img = load_untouch_nii(segmentation_path);

% Assuming you want to work with region 8
region_label = 8;

% Extract the V1 data and segmentation data as NumPy arrays
v1_data = double(v1_img.img);
segmentation_data = double(segmentation_img.img);

% Create a mask for Region 8
roi_mask = (segmentation_data == region_label);

% Apply the mask to the V1 data
v1_roi_data = v1_data .* repmat(roi_mask, [1, 1, 1, size(v1_data, 4)]);

% Calculate sin(theta) using your formula within Region 8
numerator = sum(v1_roi_data(:, :, :, 1:2) .^ 2, 4);
denominator = sum(v1_roi_data .^ 2, 4) .^ 2;
sin_theta_squared = numerator ./ denominator;

% Calculate theta within Region 8
theta_sqrt = sqrt(sin_theta_squared);
theta = (asin(theta_sqrt)) * (180 / pi);

% Create a new NIfTI image using the calculated theta within Region 8
theta_img = make_nii(theta, v1_img.hdr.dime.pixdim(2:4), [], 64);

% Save the theta image as a NIfTI file
theta_img_path = 'data/maps/theta.nii.gz';
save_nii(theta_img, theta_img_path);
%% Apparent susceptibility
segmentation_file = 'data/masks/white_matter_mask.nii.gz';
theta_file = 'data/maps/theta.nii.gz';
R1_file = 'data/maps/R1.nii.gz';
segmented_model_file = 'data/masks/SegmentedModel.nii.gz';
chi_negative_file = 'data/chimodel/Chi_negative.nii.gz'; 
Chi_positive_file='data/chimodel/Chi_positive.nii.gz';
% Output:
deltaX_file = 'data/chimodel/delta_X.nii.gz';
Xzero_file = 'data/chimodel/X_zero.nii.gz';
Xapp_file = 'data/chimodel/Xapp.nii.gz';
Xapp_weighted_file = 'data/chimodel/Xapp_weighted.nii.gz';
Chi_negative_with_anisotropy_file = 'data/chimodel/Chi_negative_with_anisotropy.nii.gz';
Chi_with_anisotropy_file = 'data/chimodel/Chi_with_anisotropy.nii.gz';

% Load the NIfTI segmentation file
img = load_untouch_nii(segmentation_file);
segmentation_data = double(img.img);

% Define the delta X and Xzero values for each region
deltaX_values = [0.032, 0.024, 0.014, 0.016, 0.016, 0.005, 0.008, 0.006, -0.015, -0.015, 0.0091];
Xzero_values = [-0.0512, -0.0522, -0.0382, -0.0512, -0.0592, -0.0442, -0.0542, -0.0462, -0.0382, -0.0372, -0.0385];

% Create delta X and Xzero maps based on region values
delta_X_map = zeros(size(segmentation_data));
Xzero_map = zeros(size(segmentation_data));

for region = 1:11
    delta_X_map(segmentation_data == region) = deltaX_values(region);
    Xzero_map(segmentation_data == region) = Xzero_values(region);
end

% Save delta X and Xzero maps as NIfTI files
delta_X_img = make_nii(delta_X_map, img.hdr.dime.pixdim(2:4), [], 64);
Xzero_img = make_nii(Xzero_map, img.hdr.dime.pixdim(2:4), [], 64);

save_nii(delta_X_img, deltaX_file);
save_nii(Xzero_img, Xzero_file);

% Load the delta X and Xzero maps
deltaX_img = load_untouch_nii(deltaX_file);
Xzero_img = load_untouch_nii(Xzero_file);
deltaX_data = double(deltaX_img.img);
Xzero_data = double(Xzero_img.img);

% Load the theta map
theta_img = load_untouch_nii(theta_file);
theta_data = double(theta_img.img);

% Calculate the Xapp map
Xapp_data = deltaX_data .* ((cos(theta_data * pi / 180)) .^ 2) + Xzero_data;

% Save the Xapp map as a NIfTI file
Xapp_img = make_nii(Xapp_data, deltaX_img.hdr.dime.pixdim(2:4), [], 64);
save_nii(Xapp_img, Xapp_file);

% Load the R1 map and segmented model
R1_img = load_untouch_nii(R1_file);
segmentation_img = load_untouch_nii(segmented_model_file);
R1_data = double(R1_img.img);
segmentation_data = double(segmentation_img.img);

% Define the region number
region_id = 8;

% Create a mask for region 8 based on the segmentation map
region_mask = (segmentation_data == region_id);

% Calculate variations for R1 in region 8
percentile_value_in_region_r1 = mean(R1_data(region_mask));
percentage_mask_r1 = region_mask .* (R1_data / percentile_value_in_region_r1);

% Define the mean and standard deviation values for Gaussian noise
mean_val = -0.04;
std_dev = 0.05;

% Initialize Xapp_masked_data as Xapp_data for the first iteration
Xapp_masked_data = Xapp_data;

for iteration = 1:3
    % Create a new noised mask for each iteration
    noise = normrnd(mean_val, std_dev, size(percentage_mask_r1));
    percentage_mask_r1_with_noise = percentage_mask_r1 + noise;

    % Multiply Xapp_masked_data with the noised mask
    Xapp_masked_data = Xapp_masked_data .* percentage_mask_r1_with_noise;
end

% Save the final Xapp image with the repeated noised mask applied
Xapp_masked_img = make_nii(Xapp_masked_data, deltaX_img.hdr.dime.pixdim(2:4), [], 64);
save_nii(Xapp_masked_img, Xapp_weighted_file);
Mask(Xapp_weighted_file);
% Load the Chi_negative image
chi_negative_img = load_untouch_nii(chi_negative_file);
chi_negative_data = double(chi_negative_img.img);

% Replace region 8 in Chi_negative with Xapp_with_R1_mask_in_region_8_3.nii.gz
chi_negative_data(segmentation_data == region_id) = Xapp_masked_data(segmentation_data == region_id);
chi_negative_data(isnan(chi_negative_data))=0;
% Save the modified Chi_negative image
chi_negative_modified_img = make_nii(chi_negative_data, chi_negative_img.hdr.dime.pixdim(2:4), [], 64);
save_nii(chi_negative_modified_img, Chi_negative_with_anisotropy_file);

% Load the Chi_negative_with_anisotropy image
chi_negative_with_anisotropy_img = load_untouch_nii(Chi_negative_with_anisotropy_file);
chi_negative_with_anisotropy_data = double(chi_negative_with_anisotropy_img.img);

% Load the Chi_positive image
chi_positive_img = load_untouch_nii(Chi_positive_file);
chi_positive_data = double(chi_positive_img.img);

% Add Chi_positive to Chi_negative_with_anisotropy
chi_with_anisotropy_data = chi_positive_data + chi_negative_with_anisotropy_data;

% Save the resulting Chi_with_anisotropy image
chi_with_anisotropy_img = make_nii(chi_with_anisotropy_data, chi_negative_with_anisotropy_img.hdr.dime.pixdim(2:4), [], 64);
save_nii(chi_with_anisotropy_img, Chi_with_anisotropy_file);

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
T2_with_variations_smoothed_img = make_nii(T2_with_variations_smoothed, T2_img.hdr.dime.pixdim(2:4), [], 64);
save_nii(T2_with_variations_smoothed_img, 'data/maps/T2.nii.gz');
R2=(1/T2_with_variations_smoothed)*1000;
R2(isnan(R2)) = 0;
R2(isinf(R2)) = 0;
R2(R2>100)=0;
% Save the T2_with_variations_smoothed map as a new NIfTI file
R2_with_variations_smoothed_img = make_nii(R2, T2_img.hdr.dime.pixdim(2:4), [], 64);
save_nii(R2_with_variations_smoothed_img, R2_file);
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