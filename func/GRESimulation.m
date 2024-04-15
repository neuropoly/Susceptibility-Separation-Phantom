function [sigHR] = GRESimulation(SequenceParam,TissueParam,SimParams)
%% reading sequence parameters

if isfield(SequenceParam,'TR')
    TR=SequenceParam.TR;
else
    TR=1;
end

if isfield(SequenceParam,'TE')
    TE=SequenceParam.TE;
else
    TE=30e-3;
end

if isfield(SequenceParam,'theta')
    theta=SequenceParam.theta;
else
    theta=90;
end

%% reading tissue model parameters

if isfield(TissueParam,'R1')
    R1=TissueParam.R1;
else
    R1=1;
end


if isfield(TissueParam,'M0')
    M0=TissueParam.M0;
else
    M0=1;
end


if isfield(TissueParam,'field')
    field=TissueParam.field;
else
    field=0;
end

if isfield(TissueParam,'PhaseOffset')
    PhaseOffset=TissueParam.PhaseOffset;
else
    PhaseOffset=0;
end

if isfield(TissueParam,'R2')
    R2=TissueParam.R2;
end

if isfield(TissueParam,'Drpos')
    Drpos=TissueParam.Drpos;
end

if isfield(TissueParam,'Drneg')
    Drneg=TissueParam.Drneg;
end

if isfield(TissueParam,'Chipos')
    Chipos=TissueParam.Chipos;
end

if isfield(TissueParam,'Chineg')
    Chineg=TissueParam.Chineg;
end
%% calculating signal 
if SimParams.B0==3 
    sigHR = (M0.*(1-exp(-TR.*R1*1.54)).*sind(theta)./(1-cosd(theta).*exp(-TR.*R1*1.54))...
    .*exp(1i.* (field * TE + PhaseOffset)) .*exp(-TE.*(R2*0.65+(Drpos*(3/7)).*abs(Chipos)+(Drneg*(3/7)).*abs(Chineg))));
    sigHR(isnan(sigHR))=0;
end 
if SimParams.B0==7 
    sigHR = (M0.*(1-exp(-TR.*R1)).*sind(theta)./(1-cosd(theta).*exp(-TR.*R1))...
    .*exp(1i.* (field * TE + PhaseOffset)) .*exp(-TE.*(R2+(Drpos).*abs(Chipos)+(Drneg).*abs(Chineg))));
    sigHR(isnan(sigHR))=0;
end 


% Parameters for Gaussian noise 
%mean_val = 0;     % Mean of the Gaussian noise
%std_dev = 20;    % Standard deviation of the Gaussian noise

% Generate Gaussian noise with the same style
%noise = normrnd(mean_val, std_dev, size(sigHR1));

% Add Gaussian noise to the signal
%sigHR= sigHR1 + noise


