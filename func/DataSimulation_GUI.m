function DataSimulation_GUI
    % Fetch the B0 field strength, TR, TE, Flip Angle, and anisotropy values from the base workspace
    B0 = evalin('base', 'B0');
    anisotropy = evalin('base', 'anisotropy');
    TR = evalin('base', 'TR')/1000;  % Fetch TR from the GUI
    % Dynamically fetch all TEs based on the number of echoes
    numEchoes = evalin('base', 'numEchoes');
    TE = zeros(1, numEchoes);  % Initialize TE array
    for i = 1:numEchoes
        teVarName = sprintf('TE%d', i);  % Dynamically create TE variable name (TE1, TE2, TE3, ...)
        TE(i) = evalin('base', teVarName) / 1000;  % Fetch TE from the GUI and convert to seconds
    end
    %TE=[0.004,0.012,0.02, 0.028];
    FlipAngle = evalin('base', 'FlipAngle');  % Fetch Flip Angle from the GUI
    
    % Path to input NIfTI files (no change)
    input_nifti_path = 'data/masks/SegmentedModel.nii.gz';
    angle_nifti_path = 'data/maps/theta.nii.gz';
    
    % Calculate Dr (no change)
    calculate_Dr(input_nifti_path, angle_nifti_path);
    
    %% Sequence protocol for simulation
    
    % Set file paths for various maps (no change)
    ModelParams.R1map_file = 'data/maps/R1.nii.gz';
    ModelParams.M0map_file = 'data/maps/M0.nii.gz';
    ModelParams.Segmentation_file = 'data/masks/SegmentedModel.nii.gz';
    ModelParams.BrainMask_file = 'data/masks/BrainMask.nii.gz';
    ModelParams.R2map_file = 'data/maps/R2.nii.gz';
    
    % Set anisotropy based on the value fetched from the GUI
    ModelParams = Anisotropy(ModelParams, anisotropy);  % true= with anisotropy, false= without anisotropy
    
    % Protocol parameters with TR, TE, and Flip Angle fetched from GUI
    Protocol = 3;
    SeqParams{Protocol}.TR = TR;  % TR from GUI                               
    SeqParams{Protocol}.TE = TE;  % TE from GUI   
    SeqParams{Protocol}.FlipAngle = FlipAngle;  % Flip Angle from GUI
    
    % Set B0 in SimParams based on the value fetched from the GUI
    SimParams{Protocol}.B0 = B0;  % B0 fetched from the GUI (3T or 7T)
    
    % Remaining parameters (no change)
    SimParams{Protocol}.B0_dir = [0 0 1];
    SimParams{Protocol}.PhaseOffset = 0;                      
    SimParams{Protocol}.Shimm = 1;                                
    SimParams{Protocol}.Res = [1 1 1];                              
    SimParams{Protocol}.Output_dir = 'Simdata/Simulation_Results'; 
    
    % Call the simulation function (no change)
    DataSimulationFunction(ModelParams, SeqParams{Protocol}, SimParams{Protocol});
    
    % Display completion message
    disp('GRE data simulation completed successfully.');
end
