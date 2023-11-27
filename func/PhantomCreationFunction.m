function PhantomCreationFunction(params)

load(params.ChiModulation_file,'label')
R12chipos = zeros(1, length(label));
R2star2chipos = zeros(1, length(label));
R12chineg = zeros(1, length(label));
R2star2chineg = zeros(1, length(label));
for k=1:length(label)
    R12chipos(k)=label(k).R12chipos;
    R2star2chipos(k)=label(k).R2star2chipos;
    R12chineg(k)=label(k).R12chineg;
    R2star2chineg(k)=label(k).R2star2chineg;
end
label(k).R2star2chi=0;
FinalSegment = load_nii(params.Segmentation_file); 
FinalSegm = FinalSegment.img;

R1 = load_nii(params.R1map_file);
tmp=Wavedec3Denoising(R1.img,25/1000,8,'db2','verysoft',double(R1.img~=0));
R1.img = (abs(real(tmp)));

R2star= load_nii(params.R2starmap_file);
tmp=Wavedec3Denoising(R2star.img,45/10,8,'db2','soft',double(R1.img~=0));
R2star.img =(abs(real(tmp)));

%% Create Chi positive
HighGradMask=load_nii(params.highGradMask_file);
RealFreqMap=load_nii(params.rawField_file);
[fieldgradient(:,:,:,1), fieldgradient(:,:,:,2), fieldgradient(:,:,:,3)]=gradient_3D(RealFreqMap.img,[],0);
fieldgradient=(sos(fieldgradient,4));
fieldgradient(FinalSegm>=11)=0;

fieldgradient=fieldgradient.*single(HighGradMask.img);
fieldgradient(and(fieldgradient~=0,fieldgradient>0.2))=0.2;
fieldgradient(and(fieldgradient~=0,fieldgradient<0.05))=0.05;
fieldgradient(fieldgradient~=0)=(fieldgradient(fieldgradient~=0)-0.05)/(0.2-0.05);
WeightingBetweenModels=cos(fieldgradient);

%%
Chimap3=0*zeros(size(FinalSegm));
Chimap4=0*zeros(size(FinalSegm));
FWHM_V=1.2;
LG_mask_index=find(or(HighGradMask.img==0,HighGradMask.img==1));

ProbabilityAccumulated=0*zeros(size(FinalSegm));

for k=1:length(label)
    
    Mask=FinalSegm==k;
    indexes=find(FinalSegm==k);
    
    if sum(Mask(:))>1
        
        Mask_smooth=real(smooth3D(double(Mask),FWHM_V,[1,1,1]));
        if k <= 10
            ProbabilityAccumulated=ProbabilityAccumulated+Mask_smooth;
            chitemp=...
                R2star2chipos(k)*(R2star.img(LG_mask_index)-mean(R2star.img(indexes)))+...
                R12chipos(k)*(R1.img(LG_mask_index)-mean(R1.img(indexes)));            
            
            NSTD=3;
            smoothThresh=0.9;
            smoothThresh2=0.05; 
            modulation_std = std(chitemp(Mask_smooth(LG_mask_index)>smoothThresh));
            modulation_mean = mean(chitemp(Mask_smooth(LG_mask_index)>smoothThresh));
           
            chitemp(and(Mask_smooth(LG_mask_index)>smoothThresh2,abs(chitemp-modulation_mean)>NSTD*modulation_std))=modulation_mean;            
            
            Chimap3(LG_mask_index)=Chimap3(LG_mask_index)+Mask_smooth(LG_mask_index)...
                .*(label(k).chipos+1*chitemp);
                       
            modulation_std = std(chitemp(Mask_smooth(LG_mask_index)>smoothThresh));
            modulation_mean = mean(chitemp(Mask_smooth(LG_mask_index)>smoothThresh));
            
            chitemp(and(Mask_smooth(LG_mask_index)>smoothThresh2,abs(chitemp-modulation_mean)>NSTD*modulation_std))=modulation_mean;
            Chimap4(LG_mask_index)=Chimap4(LG_mask_index)+Mask_smooth(LG_mask_index)...
                .*(label(k).chipos+chitemp);
            
        else
           
            chitemp=R2star2chipos(k)*(R2star.img(indexes)-mean(R2star.img(indexes)))+...
                R12chipos(k)*(R1.img(indexes)-mean(R1.img(indexes)));
            
            modulation_std = std(chitemp(Mask(indexes)>smoothThresh));
            modulation_mean = mean(chitemp(Mask(indexes)>smoothThresh));

            chitemp(and(Mask(indexes)>smoothThresh2,chitemp>modulation_mean+NSTD*modulation_std))=modulation_mean;
            chitemp(and(Mask(indexes)>smoothThresh2,chitemp<modulation_mean+NSTD*modulation_std))=modulation_mean;
            
            Chimap3(indexes)=Chimap3(indexes)+(Mask(indexes)-ProbabilityAccumulated(indexes))...
                .*(label(k).chipos+1*chitemp);
 
            modulation_std = std(chitemp(Mask(indexes)>smoothThresh));
            modulation_mean = mean(chitemp(Mask(indexes)>smoothThresh));
            chitemp(and(Mask(indexes)>smoothThresh2,chitemp>modulation_mean+NSTD*modulation_std))=modulation_mean;
            chitemp(and(Mask(indexes)>smoothThresh2,chitemp<modulation_mean+NSTD*modulation_std))=modulation_mean;
            Chimap4(indexes)=Chimap4(indexes)+(Mask(indexes)-ProbabilityAccumulated(indexes))...
                .*(label(k).chipos+chitemp);            
        end
    end
end

Chimap5=(WeightingBetweenModels.^4.*Chimap3+(1-WeightingBetweenModels.^4).*Chimap4);
chimodel=R1;
chimodel.img=Chimap5;
save_nii(chimodel,params.OutputChiModel_file)

%% Create Chi negative 
Chimap3=0*zeros(size(FinalSegm));
Chimap4=0*zeros(size(FinalSegm));
FWHM_V=1.2;
LG_mask_index=find(or(HighGradMask.img==0,HighGradMask.img==1));
ProbabilityAccumulated=0*zeros(size(FinalSegm));

for k=1:length(label)
    Mask=FinalSegm==k;
    indexes=find(FinalSegm==k);
    
    if sum(Mask(:))>1    
        Mask_smooth=real(smooth3D(double(Mask),FWHM_V,[1,1,1]));
        if k <= 10
            ProbabilityAccumulated=ProbabilityAccumulated+Mask_smooth;
            chitemp=...
                R2star2chineg(k)*(R2star.img(LG_mask_index)-mean(R2star.img(indexes)))+...
                R12chineg(k)*(R1.img(LG_mask_index)-mean(R1.img(indexes)));            
            
            NSTD=3;
            smoothThresh=0.9;
            smoothThresh2=0.05; 
            modulation_std = std(chitemp(Mask_smooth(LG_mask_index)>smoothThresh));
            modulation_mean = mean(chitemp(Mask_smooth(LG_mask_index)>smoothThresh));
           
            chitemp(and(Mask_smooth(LG_mask_index)>smoothThresh2,abs(chitemp-modulation_mean)>NSTD*modulation_std))=modulation_mean;            
            
            Chimap3(LG_mask_index)=Chimap3(LG_mask_index)+Mask_smooth(LG_mask_index)...
                .*(label(k).chineg+1*chitemp);
                       
            modulation_std = std(chitemp(Mask_smooth(LG_mask_index)>smoothThresh));
            modulation_mean = mean(chitemp(Mask_smooth(LG_mask_index)>smoothThresh));
            
            chitemp(and(Mask_smooth(LG_mask_index)>smoothThresh2,abs(chitemp-modulation_mean)>NSTD*modulation_std))=modulation_mean;
            Chimap4(LG_mask_index)=Chimap4(LG_mask_index)+Mask_smooth(LG_mask_index)...
                .*(label(k).chineg+chitemp);
            
        else
           
            chitemp=R2star2chineg(k)*(R2star.img(indexes)-mean(R2star.img(indexes)))+...
                R12chineg(k)*(R1.img(indexes)-mean(R1.img(indexes)));
            
            modulation_std = std(chitemp(Mask(indexes)>smoothThresh));
            modulation_mean = mean(chitemp(Mask(indexes)>smoothThresh));

            chitemp(and(Mask(indexes)>smoothThresh2,chitemp>modulation_mean+NSTD*modulation_std))=modulation_mean;
            chitemp(and(Mask(indexes)>smoothThresh2,chitemp<modulation_mean+NSTD*modulation_std))=modulation_mean;
            
            Chimap3(indexes)=Chimap3(indexes)+(Mask(indexes)-ProbabilityAccumulated(indexes))...
                .*(label(k).chineg+1*chitemp);
 
            modulation_std = std(chitemp(Mask(indexes)>smoothThresh));
            modulation_mean = mean(chitemp(Mask(indexes)>smoothThresh));
            chitemp(and(Mask(indexes)>smoothThresh2,chitemp>modulation_mean+NSTD*modulation_std))=modulation_mean;
            chitemp(and(Mask(indexes)>smoothThresh2,chitemp<modulation_mean+NSTD*modulation_std))=modulation_mean;
            Chimap4(indexes)=Chimap4(indexes)+(Mask(indexes)-ProbabilityAccumulated(indexes))...
                .*(label(k).chineg+chitemp);            
        end
    end
end

Chimap5=(WeightingBetweenModels.^4.*Chimap3+(1-WeightingBetweenModels.^4).*Chimap4);

chimodelneg=R1;
chimodelneg.img=Chimap5;
save_nii(chimodelneg,'data/chimodel/Chi_negative.nii.gz');


%% Create Chi total
%Calculation of ChiTotal from ChiNeg and ChiPos
% File paths of the two NIfTI images to be loaded
ChiNegPath = 'data/chimodel/Chi_negative.nii.gz';
ChiPosPath = 'data/chimodel/Chi_positive.nii.gz';

% Load the first NIfTI image
chineg = load_untouch_nii(ChiNegPath);
chinegdata = chineg.img;

% Load the second NIfTI image
chipos = load_untouch_nii(ChiPosPath);
chiposdata = chipos.img;

% Perform addition of the two images
Chitotal = chinegdata + chiposdata;

% Save the result as a NIfTI image
resultNii = chineg;  % Use image1 as a template for header information
resultNii.img = Chitotal;  % Assign the result to the image data

% File path to save the result NIfTI image
resultPath = 'data/chimodel/Chi.nii.gz';

% Save the result NIfTI image
save_untouch_nii(resultNii, resultPath);

% Display a message upon successful save
disp('Result NIfTI image saved successfully.');