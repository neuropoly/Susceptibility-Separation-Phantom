% This code Create an R2 prime map (used for X-separation) (see README.md for more details). 
%
%%%%%%%%%%%%%%%%%%%%% VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%
%%% Inputs: %%% 
% T2star_img- 3D T2* map
% T2_img- 3D T2 map (created using the CreateT2map.m script)
%%% Output: %%%
%
% R2prime_img-  3D R2' map
%%%%%%%%%%%%%%%%%%%%% VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%
% Load the input images
% Input:
T2star_img = load_untouch_nii('Path/to/T2_star.nii.gz');
T2_img = load_untouch_nii('Path/to/T2.nii.gz');
% Output:
R2star_file = 'Path/to/R2_star.nii.gz';
R2_file = 'Path/to/R2.nii.gz';
R2prime_file = 'Path/to/R2_prime.nii.gz';

% Calculate the inverse of each image
inv_T2star_img = 1./double(T2star_img.img);
inv_T2_img = 1./double(T2_img.img);

% Replace any inf or NaN values with zero
inv_T2star_img(~isfinite(inv_T2star_img)) = 0;
inv_T2_img(~isfinite(inv_T2_img)) = 0;

inv_T2star_img = squeeze(inv_T2star_img);
inv_T2_img = squeeze(inv_T2_img);

% Save R2_star and R2
R2star_img = make_nii(inv_T2star_img, T2star_img.hdr.dime.pixdim(2:4), [], 64);
save_nii(R2star_img, R2star_file);

R2_img = make_nii(inv_T2_img, T2_img.hdr.dime.pixdim(2:4), [], 64);
save_nii(R2_img, R2_file);

% Calculate the difference between the two images
R2prime_img = inv_T2star_img - inv_T2_img;
% Set any negative values to zero
R2prime_img(R2prime_img < 0) = 0;

% Save the result to a NIfTI file
result_img = make_nii(R2prime_img, T2star_img.hdr.dime.pixdim(2:4), [], 64);
save_nii(result_img, R2prime_file);
