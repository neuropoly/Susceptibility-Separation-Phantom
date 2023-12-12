function [sigHR] = GRESimulation(SequenceParam,TissueParam)
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
%% calculating signal at full resolution

sigHR = M0.*(1-exp(-TR.*R1)).*sind(theta)./(1-cosd(theta).*exp(-TR.*R1))...
.*exp(1i.* (field * TE + PhaseOffset)) .*exp(-TE.*(R2+(Drpos).*abs(Chipos)+(Drneg).*abs(Chineg)));
sigHR(isnan(sigHR))=0;





