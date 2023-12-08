function ModelParams = Anisotropy(ModelParams, AnisotropyFlag)
    % Check if AnisotropyFlag is true or false
    if AnisotropyFlag
        disp('Simulating GRE data with anisotropy');
        ModelParams.Chimap_file = 'data/chimodel/Chi_with_anisotropy.nii.gz';
        ModelParams.Chiposmap_file = 'data/chimodel/Chi_positive.nii.gz';
        ModelParams.Chinegmap_file = 'data/chimodel/Chi_negative_with_anisotropy.nii.gz';
        ModelParams.Drposmap_file = 'data/maps/Dr_positive.nii.gz';
        ModelParams.Drnegmap_file = 'data/maps/Dr_negative_variable.nii.gz';
    else
        disp('Simulating GRE data without anisotropy');
        ModelParams.Chimap_file = 'data/chimodel/Chi.nii.gz';
        ModelParams.Chiposmap_file = 'data/chimodel/Chi_positive.nii.gz';
        ModelParams.Chinegmap_file = 'data/chimodel/Chi_negative.nii.gz';
        ModelParams.Drposmap_file = 'data/maps/Dr_positive.nii.gz';
        ModelParams.Drnegmap_file = 'data/maps/Dr_negative_constant.nii.gz';
    end
end
