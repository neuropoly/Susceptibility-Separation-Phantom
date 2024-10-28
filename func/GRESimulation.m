function [sigHR] = GRESimulation(SequenceParam, TissueParam, SimParams)
    %% Reading sequence parameters
    if isfield(SequenceParam, 'TR')
        TR = SequenceParam.TR;
    else
        TR = 1;
    end

    if isfield(SequenceParam, 'TE')
        TE = SequenceParam.TE;
    else
        TE = 30e-3;
    end

    if isfield(SequenceParam, 'theta')
        theta = SequenceParam.theta;
    else
        theta = 90;
    end

    %% Reading tissue model parameters
    if isfield(TissueParam, 'R1')
        R1 = TissueParam.R1;
    else
        R1 = 1;
    end

    if isfield(TissueParam, 'M0')
        M0 = TissueParam.M0;
    else
        M0 = 1;
    end

    if isfield(TissueParam, 'field')
        field = TissueParam.field;
    else
        field = 0;
    end

    if isfield(TissueParam, 'PhaseOffset')
        PhaseOffset = TissueParam.PhaseOffset;
    else
        PhaseOffset = 0;
    end

    if isfield(TissueParam, 'R2')
        R2 = TissueParam.R2;
    else
        R2 = 1;
    end

    if isfield(TissueParam, 'Drpos')
        Drpos = TissueParam.Drpos;
    else
        Drpos = 0;
    end

    if isfield(TissueParam, 'Drneg')
        Drneg = TissueParam.Drneg;
    else
        Drneg = 0;
    end

    if isfield(TissueParam, 'Chipos')
        Chipos = TissueParam.Chipos;
    else
        Chipos = 0;
    end

    if isfield(TissueParam, 'Chineg')
        Chineg = TissueParam.Chineg;
    else
        Chineg = 0;
    end

    %% Load R1_3T map if B0 == 3 T
    if SimParams.B0 == 3
        % Load R1_3T map
        R1nii = load_nii('data/maps/R1_3T.nii.gz');
        R1 = R1nii.img;
    end

    %% Calculating signal
    if SimParams.B0 == 3
        % Remove multiplication by 1.54 for R1
        sigHR1 = (M0 .* (1 - exp(-TR .* R1)) .* sind(theta) ./ (1 - cosd(theta) .* exp(-TR .* R1)) ...
            .* exp(1i * (field * TE + PhaseOffset)) ...
            .* exp(-TE .* (R2 * 0.65 + (Drpos * (3 / 7)) .* abs(Chipos) + (Drneg * (3 / 7)) .* abs(Chineg))));
        sigHR1(isnan(sigHR1)) = 0;
    elseif SimParams.B0 == 7
        sigHR1 = (M0 .* (1 - exp(-TR .* R1)) .* sind(theta) ./ (1 - cosd(theta) .* exp(-TR .* R1)) ...
            .* exp(1i * (field * TE + PhaseOffset)) ...
            .* exp(-TE .* (R2 + Drpos .* abs(Chipos) + Drneg .* abs(Chineg))));
        sigHR1(isnan(sigHR1)) = 0;
    else
        error('Unsupported B0 field strength.');
    end

    % Parameters for Gaussian noise
    mean_val = 0;    % Mean of the Gaussian noise
    std_dev = 1;     % Standard deviation of the Gaussian noise

    % Generate Gaussian noise with the same size as the signal
    noise = normrnd(mean_val, std_dev, size(sigHR1));

    % Add Gaussian noise to the signal
    sigHR = sigHR1 + noise;
end
