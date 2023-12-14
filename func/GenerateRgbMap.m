% This code create an RGB directional color coded map  (see README.md for more details). 
%
%%%%%%%%%%%%%%%%%%%%% VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%
%%% Inputs: %%% 
% 
% img_fa- Fractional anisotropy map from DTI data
% data_v1- Principal vectore map from DTI data
%%% Output: %%%
%
% rgb- RGB map of the FA data ( color codded according to the direction of the white matter fibers)
%%%%%%%%%%%%%%%%%%%%% VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%
% File paths and names
path = 'Path/to';

img_fa = load_untouch_nii(fullfile(path, 'FA.nii'));
data_fa = double(img_fa.img);

data_v1 = load_untouch_nii(fullfile(path, 'V1.nii.gz'));
data_v1 = double(data_v1.img);

% Display shapes
disp(size(data_v1));
disp(size(data_fa));

% Color encoded
rgb = abs(data_v1) .* repmat(max(0, min(1, data_fa)), [1, 1, 1, 3]);
rgb_nii = make_nii(rgb, img_fa.hdr.dime.pixdim(2:4), [], 64);

% Save the RGB NIfTI image
save_nii(rgb_nii, fullfile(path, 'rgb.nii.gz'));
