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
Mask('data/chimodel/Chi_positive.nii.gz');
Mask('data/chimodel/Chi.nii.gz');
Mask('data/chimodel/Chi_negative.nii.gz');
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
%Xapp_masked_img = make_nii(Xapp_masked_data, deltaX_img.hdr.dime.pixdim(2:4), [], 64);
%save_nii(Xapp_masked_img, Xapp_weighted_file);
%Mask(Xapp_weighted_file);
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
