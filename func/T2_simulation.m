function T2_simulation(seg_file, T2_star_file, M0_file, T2_uniform_file, R2_file)
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
    segmentation_img = load_untouch_nii(seg_file);

    % Get the data arrays
    T2_star_data = double(T2_star_img.img);
    T2_data = double(T2_img.img);
    M0_data = double(M0_img.img);
    segmentation_data = double(segmentation_img.img);

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

    % Calculate R2 from T2
    R2 = (1 ./ T2_with_variations_smoothed) * 1000;
    R2(isnan(R2)) = 0;
    R2(isinf(R2)) = 0;
    R2(R2 > 300) = 0;

    % Save the R2 map as a NIfTI file
    R2_with_variations_smoothed_img = make_nii(R2, T2_img.hdr.dime.pixdim(2:4), [], 64);
    save_nii(R2_with_variations_smoothed_img, R2_file);
        % Calculate T2_3T 
    T2_3T = T2_with_variations_smoothed / 0.65;

    % Save the R2_3T map as a NIfTI file
    T2_3T_img = make_nii(T2_3T, T2_img.hdr.dime.pixdim(2:4), [], 64);
    save_nii(T2_3T_img,'data/maps/T2_3T.nii.gz' );


    % Apply mask (assuming Mask function is provided)
    Mask(R2_file);
    Mask('data/maps/T2_3T.nii.gz');
end
