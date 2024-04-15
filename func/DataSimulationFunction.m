function DataSimulationFunction(ModelParams,SeqParams,SimParams)

Sufix(1).name='BrainExtracted';
Sufix(2).name='';
savedata = 1;
B0=7;
B0_dir=[0 0 1];
gyro = 42.57747892;
modelname = SimParams.Output_dir ;

    for BackGroundFieldRemoval = [1 2]

        M0=load_nii(ModelParams.M0map_file);
        R1=load_nii(ModelParams.R1map_file);
        chi0=load_nii(ModelParams.Chimap_file);
        Brain=load_nii(ModelParams.BrainMask_file);
        Brain=Brain.img;
        voxel_size = round(M0.hdr.dime.pixdim(2:4)*100)/100;
        R2=load_nii(ModelParams.R2map_file);
        Drpos=load_nii(ModelParams.Drposmap_file);
        Drneg=load_nii(ModelParams.Drnegmap_file);
        Chipos=load_nii(ModelParams.Chiposmap_file);
        Chineg=load_nii(ModelParams.Chinegmap_file);
        chi = chi0 ;
               
        if BackGroundFieldRemoval == 1
            chi.img = (chi.img - mean(chi.img(Brain == 1))) .* Brain ;
            Shimm = 0 ; 
        end
        
        if BackGroundFieldRemoval == 2
            Shimm = SimParams.Shimm ;
        end
        %% creates dipole Kernel
        dims = size(M0.img);       
        D = create_dipole_kernel(B0_dir , voxel_size, 2 * dims, 1);
        
        chitemp = ones( 2 * dims) * chi.img(end,end,end);
        chitemp (1:dims(1),1:dims(2),1:dims(3)) = chi.img;
        field= real(ifftn(fftn(chitemp).*(D)))  ;
        
        field = field(1:dims(1),1:dims(2),1:dims(3));
        clear chitemp
        clear D
        
        %% Brain shimming
        
        if Shimm==0
            fieldb.img=field;
        else
            [~,fieldb,~]=PolyFitShimLike(make_nii(field),make_nii(single(Brain)),2);
        end
        %% Brain Phase Offset
        
        if SimParams.PhaseOffset==0 
            PhaseOffset = 0;
        else
            [c , w ] = centerofmass(M0.img);
            [y,x,z] = meshgrid((1:dims(2))-c(2), (1:dims(1))-c(1), (1:dims(3))-c(3));
            temp = (x/w(1)).^2 + (y/w(2)).^2 + (z/w(3)).^2 ;
            PhaseOffset = - temp/(max(temp(Brain==1))-min(temp(Brain==1)))*pi*SimParams.PhaseOffset;
        end
            
        %%
        
        for TE = SeqParams.TE
            
            if savedata == 1
               
                FolderHR=[modelname,'/SimulatedHR/'];
               
                
                if exist(FolderHR,'dir')~=7
                    mkdir(FolderHR)
                end
                
            end
            
            SequenceParam.TR=SeqParams.TR;                     
            SequenceParam.TE=TE;                     
            SequenceParam.theta=SeqParams.FlipAngle; 
            TissueParam.B0=SimParams.B0;
            if length(SimParams.Res)==1
                SequenceParam.res=[1 1 1 ]*SimParams.Res;             
            elseif length(SimParams.Res)==3
                SequenceParam.res=SimParams.Res;             
            else
                disp('you have not entered the resolution in the correct format, we [ 1 1 1 ] will be assumed')
                SequenceParam.res=[1 1 1 ];             
            end
            TissueParam.M0=double(M0.img);              
            TissueParam.R1=double(R1.img);                                    
            TissueParam.field=field * SimParams.B0 * gyro*2*pi;                  
            TissueParam.R2=R2.img; 
            TissueParam.Drpos=Drpos.img; 
            TissueParam.Drneg=Drneg.img;
            TissueParam.Chipos=Chipos.img;
            TissueParam.Chineg=Chineg.img;
            TissueParam.PhaseOffset = PhaseOffset;             
            TissueParam.res = voxel_size;
           
            [sigHR]= GRESimulation(SequenceParam,TissueParam); 
             
            if savedata == 1
                save_nii(make_nii(abs(sigHR),voxel_size),[FolderHR,'Magnitude_TE',num2str(TE*1000,'%03.f'),Sufix(BackGroundFieldRemoval).name,'.nii.gz'])
                save_nii(make_nii(angle(sigHR),voxel_size),[FolderHR,'Phase_TE',num2str(TE*1000,'%03.f'),Sufix(BackGroundFieldRemoval).name,'.nii.gz'])
            end
            %Mask([FolderHR,'Magnitude_TE',num2str(TE*1000,'%03.f'),Sufix(BackGroundFieldRemoval).name,'.nii.gz']);
            %Mask([FolderHR,'Phase_TE',num2str(TE*1000,'%03.f'),Sufix(BackGroundFieldRemoval).name,'.nii.gz']);
        end
      
        
    end
   % Initialize cell arrays to store the file names
magnitude_files = cell(1, length(SeqParams.TE));
phase_files = cell(1, length(SeqParams.TE));

for i = 1:length(SeqParams.TE)
    TE_str = num2str(SeqParams.TE(i) * 1000, '%03.f');

    % Construct file names for each echo
    magnitude_files{i} = ['Simdata/Simulation_Results/SimulatedHR/Magnitude_TE', TE_str, '.nii.gz'];
    phase_files{i} = ['Simdata/Simulation_Results/SimulatedHR/Phase_TE', TE_str, '.nii.gz'];

    % Load NIfTI images for each echo
    magnitude_echoes{i} = load_untouch_nii(magnitude_files{i});
    phase_echoes{i} = load_untouch_nii(phase_files{i});
end

% Extract the image data and dimensions for magnitude echoes
data_magnitude = double(magnitude_echoes{1}.img);
% Concatenate the data by grouping the echoes
grouped_data_magnitude = cat(4, data_magnitude, double(magnitude_echoes{2}.img), double(magnitude_echoes{3}.img), double(magnitude_echoes{4}.img));

% Create a new NIfTI image using the grouped data
grouped_nii_magnitude = make_nii(grouped_data_magnitude, magnitude_echoes{1}.hdr.dime.pixdim(2:4), [], 64);

% Save the new NIfTI image to a file
save_nii(grouped_nii_magnitude, 'Simdata/Simulation_Results/SimulatedHR/Magnitude_data.nii.gz');

% Extract the image data and dimensions for phase echoes
data_phase = double(phase_echoes{1}.img);
% Concatenate the data by grouping the echoes
grouped_data_phase = cat(4, data_phase, double(phase_echoes{2}.img), double(phase_echoes{3}.img), double(phase_echoes{4}.img));

% Create a new NIfTI image using the grouped data
grouped_nii_phase = make_nii(grouped_data_phase, phase_echoes{1}.hdr.dime.pixdim(2:4), [], 64);

% Save the new NIfTI image to a file
save_nii(grouped_nii_phase, 'Simdata/Simulation_Results/SimulatedHR/Phase_data.nii.gz');

    save( [ modelname, '/SimulationParameters.mat'],'ModelParams','SeqParams','SimParams')
 
