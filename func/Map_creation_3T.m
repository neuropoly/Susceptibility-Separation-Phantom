function Map_creation_3T(R2_file, R2star_file, R1_file, seg_file)
    % Load R2 and multiply by 0.65, then save as R2_3T
    R2_img = load_untouch_nii(R2_file);
    R2_data = double(R2_img.img);
    R2_data_3T = R2_data * 0.65;
    R2_3T_img = make_nii(R2_data_3T, R2_img.hdr.dime.pixdim(2:4), [], 64);
    save_nii(R2_3T_img, 'data/maps/R2_3T.nii.gz');

    % Load R2star and multiply by 0.5, then save as R2star_3T
    R2star_img = load_untouch_nii(R2star_file);
    R2star_data = double(R2star_img.img);
    R2star_data_3T = R2star_data * 0.5;
    R2star_3T_img = make_nii(R2star_data_3T, R2star_img.hdr.dime.pixdim(2:4), [], 64);
    save_nii(R2star_3T_img, 'data/maps/R2star_3T.nii.gz');

    % Load R1 and Segmented Model
    R1_img = load_untouch_nii(R1_file);
    R1_data = double(R1_img.img);
    seg_img = load_untouch_nii(seg_file);
    seg_data = double(seg_img.img);

    % Division factors for regions 1 to 11
    division_factors = [0.75929, 0.73274, 0.74212, 0.65, 0.65, 0.65, 0.73898, 0.72472, 0.73648, 1.0051, 0.75672];

    % Apply division to each region
    output_data = R1_data;
    for label = 1:11
        mask = (seg_data == label);
        output_data(mask) = R1_data(mask) / division_factors(label);
    end

    % Save R1_3T
    R1_3T_img = make_nii(output_data, R1_img.hdr.dime.pixdim(2:4), [], 64);
    save_nii(R1_3T_img, 'data/maps/R1_3T.nii.gz');
end
