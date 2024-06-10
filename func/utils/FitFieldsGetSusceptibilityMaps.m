function [FittedSuscep, FittedField]=FitFieldsGetSusceptibilityMaps(Field,SegmentationModel,ExtraModulation, B0dir,voxel_size,shimm,varargin)
% [FittedSuscep, FittedField]=FitFieldsGetSuisceptibilityMaps(Field,SegmentationModel,ExtraModulation, B0dir,voxel_size,shimm)
% FieldSegmentationModel
% ExtraModulation can be left empty or be a 4D volume with R1, R2*, T1w, MT and blabla
% B0dir of this model
% voxel_size
% shim, can be 0 or 1 meaning that other than the field it fits a second
% order polynominal
% if nargin = 7 structure with tissue with stepwise smooth assumptions (modulation is not used ) 
% and no stepwosesmooth assumption (Only the modulation models are used there)
%    TissueModelsNoStepWiseSmooth=   varargin{1}.ModulationOnly;
%    TissueModelsWithStepWiseSmooth= varargin{1}.StepWiseSmooth;

if isempty(ExtraModulation)
    models=SegmentationModel;
else
    models=SegmentationModel;
    for k=1:size(ExtraModulation,4)
        models=cat(4, models, bsxfun(@times,SegmentationModel,ExtraModulation(:,:,:,k)));
    end;
end;


%%

if nargin == 7
    NModulations = size(ExtraModulation,4);
    NTissues = size(SegmentationModel,4);
    SusceptibilitySources=size(models,4);

    if isfield(varargin{1},'ModulationOnly')
        TissueModelsNotStepWiseSmooth=   varargin{1}.ModulationOnly;
    else
        TissueModelsNotStepWiseSmooth= []
    end;

    if isfield(varargin{1},'StepWiseSmooth')
        temp=   varargin{1}.StepWiseSmooth;
        TissueModelsStepWiseSmooth=[];
        for k=1:NModulations
        TissueModelsStepWiseSmooth = cat(2,TissueModelsStepWiseSmooth, NTissues+temp);
        end
    else
        TissueModelsStepWiseSmooth= []
    end;
    SusceptibilitySources=size(models,4);
   Models2keep = setdiff(1:SusceptibilitySources,TissueModelsNotStepWiseSmooth);

   Models2keep = setdiff(Models2keep,TissueModelsStepWiseSmooth);

   models=models(:,:,:,Models2keep);
   
else
    
    
    
end;

SusceptibilitySources=size(models,4);

N=size(Field);


mask=Field~=0;
Indices=find(mask);



%%
D = create_dipole_kernel(B0dir, voxel_size,N, 1);

modelfields=models;
clear modelfieldvector
for k=1:size(models,4)
    temp=ifftn(fftn(modelfields(:,:,:,k)).*D);
    modelfields(:,:,:,k)=temp;
    modelfieldvector(:,k)=temp(Indices);
end
if shimm==0
    modelfields(:,:,:,end+1)=1;
    modelfieldvector(:,end+1)=1;
else
    
    [x1,y1,z1]=ndgrid((1:N(1))-(N(1)+1)/2,(1:N(2))-(N(2)+1)/2,(1:N(3))-(N(3)+1)/2);
    
    modelfields(:,:,:,end+1)=1;
    modelfieldvector(:,end+1)=1;
    modelfields(:,:,:,end+1)=x1;
    modelfieldvector(:,end+1)=x1(Indices);
    modelfields(:,:,:,end+1)=y1;
    modelfieldvector(:,end+1)=y1(Indices);
    modelfields(:,:,:,end+1)=z1;
    modelfieldvector(:,end+1)=z1(Indices);
    modelfields(:,:,:,end+1)=x1.^2;
    modelfieldvector(:,end+1)=x1(Indices).^2;
    modelfields(:,:,:,end+1)=y1.^2;
    modelfieldvector(:,end+1)=y1(Indices).^2;
    modelfields(:,:,:,end+1)=z1.^2;
    modelfieldvector(:,end+1)=z1(Indices).^2;
    modelfields(:,:,:,end+1)=x1.*y1;
    modelfieldvector(:,end+1)=x1(Indices).*y1(Indices);
    modelfields(:,:,:,end+1)=y1.*z1;
    modelfieldvector(:,end+1)=y1(Indices).*z1(Indices);
    modelfields(:,:,:,end+1)=z1.*x1;
    modelfieldvector(:,end+1)=z1(Indices).*x1(Indices);
    
    
end;

% try
%  MultFact=Field(Indices)'*pinv(modelfieldvector');
% catch

% if nargin<7
MultFact=lsqr(modelfieldvector,Field(Indices),1e-6,1000);
% else
%     in this case we put an weighting factor on the fit
% MultFact=lsqr(modelfieldvector,varargin{1}.*Field(Indices),1e-4,100);
%     
% % end

FittedField=real(sum(bsxfun(@times,reshape(MultFact,[1 1 1 length(MultFact)]),modelfields),4));
FittedSuscep=real(sum(bsxfun(@times,reshape(MultFact(1:(SusceptibilitySources+1)),[1 1 1 SusceptibilitySources+1]),cat(4,models,ones(N))),4));
% fv(real(FittedSuscep))
% fv(real(FittedField))

ShowResults=1;
if ShowResults==1
    dims=size(Field)
    pos=round(dims/2)+[+8 5 -9];
    figure()
%     subplot(321)
    subplot1(3,2,'Gap',[0.02 0.02])
    subplot1(1)
    Orthoview(FittedSuscep,pos)
    title('Susceptibility map ')
%     subplot(322)
    subplot1(2)

    % Orthoview(-FittedSuscep,pos,[ 50 600])
temp=FittedSuscep(Field~=0)-mean(FittedSuscep((Field~=0)));
    Orthoview(FittedSuscep-mean(FittedSuscep((Field~=0))),pos, prctile(temp(:),[4 96]));
    title('Susceptibility map zoomed')
%     subplot(323)
    subplot1(3)

    Orthoview(Field,pos,[-200 40]);
    title('Data to fit')
    subplot1(4)
%     subplot(324)
    Orthoview(FittedField.*single(Field~=0),pos,[-200 40]);
    title('fitted field')
%     subplot(325)
     subplot1(5)
    Orthoview(Field-FittedField.*single(Field~=0),pos,[-200 40]);
    title('residuals of fit')    
end;

