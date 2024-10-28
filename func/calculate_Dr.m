function [Dr_pos_file, Dr_neg_variable_file, Dr_neg_constant_file] = calculate_Dr(input_nifti_path, angle_nifti_path)
    % Define output file paths
    Dr_pos_file = 'data/maps/Dr_positive.nii.gz';
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

    % Calculate Dr_negative using the formula
    Drneg_data = 0.5 * 42.58 * 2 * pi * 7 * sin(angle_data_radians).^2;
    Drneg_data(isnan(Drneg_data)) = 0;

    % Create a new NIfTI image for Dr_negative
    Drneg_image = make_nii(Drneg_data, angle_image.hdr.dime.pixdim(2:4), [], 64);
    save_nii(Drneg_image, Dr_neg_variable_file);
    Mask(Dr_neg_variable_file);

    % Now for Dr negative constant
    % Load the segmentation NIfTI image again
    data2 = double(image.img);
    
    % Create a new NIfTI image with a constant value (700.8) for Dr_negative constant
    value_constant = 700.8;
    
    % Set the specified regions to the desired value for Dr_negative constant
    for region = regions_to_modify
        data2(data2 == region) = 0;
    end
    data2(data2 == 8) = value_constant;

    % Create a new NIfTI image for Dr_negative constant
    Drneg_image_constant = make_nii(data2, angle_image.hdr.dime.pixdim(2:4), [], 64);
    save_nii(Drneg_image_constant, Dr_neg_constant_file);
    Mask(Dr_neg_constant_file);
end
