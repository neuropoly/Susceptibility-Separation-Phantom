function Mask(input_brain_path)
    % Hardcoded mask image path
    mask_path = 'data/masks/BrainMask.nii.gz';

    % Load the brain image and its header
    brain_nii = load_nii(input_brain_path);
    brain_image = brain_nii.img;
    brain_image(isnan(brain_image)) = 0;
    
    % Load the mask image and its header
    mask_nii = load_nii(mask_path);
    mask_image = mask_nii.img;

    % Apply the mask to the brain image
    masked_image = brain_image .* mask_image;
    masked_image(isnan(masked_image)) = 0;
    
    % Create a new NIfTI structure for the masked image with the same header as the brain image
    masked_nii = mask_nii;
    masked_nii.img = masked_image;

    % Save the masked image as a new NIfTI file
    save_nii(masked_nii, input_brain_path);
end
