% This code simulate gradient echo (GRE) data: Magnitude and phase ( total, positive, and negative) (see README.md
% for more details). 

%%%%%%%%%%%%%%%%%%%%% VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%
%%% Inputs: %%% 
% Chimap_file- 3D total susceptibility map
% Chiposmap_file- 3D positive susceptibility map
% Chinegmap_file- 3D negative susceptibility map
% R1map_file- 3D R1 map (1/T1).
% M0map_file- 3D Net magnetization map (M0).             
% Segmentation_file- 3D brain segmentation file.         
% BrainMask_file- 3D maks of ROI.
% R2map_file- 3D R2star map (1/T2).
% Drposmap_file- 3D positive Dr.
% Drnegmap_file- 3D negative Dr.
% TR- Repition Time in seconds.
% TE- Echo Time in seconds (array).
% FlipAngle- Flip angle in degrees above ernst angle (13) to improve gre contrast between Grey and white matter.
% B0- Main magnetic field strength in Tesla.
% B0_dir- B0 direction .
% PhaseOffset- Phase Offset 0 no phase offset; pi phase difference inside brain mask. 
% Shim- boolean 0 if no additional shimms are applied 1 if
% Res- resolution of output
%%% Output: %%%
%
% Output_dir- Path to results ( 4D Simulated GRE magnitude and phase)
%%%%%%%%%%%%%%%%%%%%% VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%

%% Sequence protocol for simulation

ModelParams.R1map_file = 'data/maps/R1.nii.gz';
ModelParams.M0map_file = 'data/maps/M0.nii.gz';
ModelParams.Segmentation_file = 'data/masks/SegmentedModel.nii.gz';
ModelParams.BrainMask_file = 'data/masks/BrainMask.nii.gz';
ModelParams.R2map_file='data/maps/R2.nii.gz';

ModelParams = Anisotropy(ModelParams, false);%true= simulate data with anisotropy, false: simulate data without anisotropy

Protocol = 3 ;
SeqParams{Protocol}.TR = 50e-3;                               
SeqParams{Protocol}.TE = [4e-3 12e-3  20e-3  28e-3  ];     
SeqParams{Protocol}.FlipAngle = 15;                             
SimParams{Protocol}.B0 = 7;
SimParams{Protocol}.B0_dir = [0 0 1];
SimParams{Protocol}.PhaseOffset = 0 ;                      
SimParams{Protocol}.Shimm = 1;                                  
% extracted output
SimParams{Protocol}.Res = [1 1 1];                              
SimParams{Protocol}.Output_dir = 'Simdata/Simulation_Results'; 

DataSimulationFunction(ModelParams,SeqParams{Protocol},SimParams{Protocol})




