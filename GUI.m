function GUI
    % Create the main figure
    hFig = figure('Position', [100, 100, 900, 700], 'MenuBar', 'none', ...
                  'Name', 'MRI Phantom and Simulation', 'NumberTitle', 'off', 'Resize', 'on');
              
    % Create a tab group
    tabGroup = uitabgroup(hFig);
    
    % Create tabs for Phantom Creation and Data Simulation
    creationTab = uitab(tabGroup, 'Title', 'Phantom Creation');
    simulationTab = uitab(tabGroup, 'Title', 'Data Simulation');
    
    %% Adjusting the Phantom Creation Tab Layout
    generateMapsPanel = uipanel('Parent', creationTab, 'Title', 'Generate Maps', ...
                                'Position', [0.03, 0.88, 0.95, 0.1]);
    
    % Dropdown for field strength
    uicontrol('Parent', generateMapsPanel, 'Style', 'text', ...
              'Position', [10 20 100 20], 'String', 'Field Strength:', 'FontSize', 10);
    fieldStrengthPopup = uicontrol('Parent', generateMapsPanel, 'Style', 'popupmenu', 'Position', [120 20 70 20], ...
              'String', {'3T', '7T'}, 'FontSize', 10);
    
    % Checkbox for anisotropy and enlarged Generate Maps button
    anisotropyCheckbox = uicontrol('Parent', generateMapsPanel, 'Style', 'checkbox', ...
                                   'Position', [350 20 150 20], 'String', 'Include Anisotropy', ...
                                   'FontSize', 10);
    generateMapsButton = uicontrol('Parent', generateMapsPanel, 'Style', 'pushbutton', 'Position', [680 12 150 40], ...
              'String', 'Generate Maps', 'FontSize', 12, 'FontWeight', 'bold', 'Callback', @generateMaps);
    
    % Panels for Field Strength Dependent Maps and Susceptibility Maps
    fieldStrengthPanel = uipanel('Parent', creationTab, 'Title', 'Field Strength Dependent Maps', ...
                                 'Position', [0.03, 0.02, 0.47, 0.85]);
    
    susceptibilityPanel = uipanel('Parent', creationTab, 'Title', 'Susceptibility Maps', ...
                                  'Position', [0.51, 0.02, 0.47, 0.85]);
    
    % Initially no axes created until maps are generated
    r2Axes = [];
    r1Axes = [];
    r2StarAxes = [];
    chiPlusAxes = [];
    chiMinusAxes = [];
    chiTotalAxes = [];
    
    % Add Colormap Selection to fieldStrengthPanel
    uicontrol('Parent', fieldStrengthPanel, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.01, 0.95, 0.15, 0.04], 'String', 'Colormap:', 'FontSize', 10);

    colormapPopupField = uicontrol('Parent', fieldStrengthPanel, 'Style', 'popupmenu', 'Units', 'normalized', ...
        'Position', [0.15, 0.94, 0.2, 0.05], 'String', {'Jet', 'Grayscale', 'Blue to White', 'Red to White', 'Blue-White-Red'}, ...
        'FontSize', 10, 'Tag', 'ColormapPopup', 'Callback', @(src, event) colormapSelectionChanged(src, event, fieldStrengthPanel));

    % Add Colormap Selection to susceptibilityPanel
    uicontrol('Parent', susceptibilityPanel, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.01, 0.95, 0.15, 0.04], 'String', 'Colormap:', 'FontSize', 10);

    colormapPopupSusceptibility = uicontrol('Parent', susceptibilityPanel, 'Style', 'popupmenu', 'Units', 'normalized', ...
        'Position', [0.15, 0.94, 0.2, 0.05], 'String', {'Jet', 'Grayscale', 'Blue to White', 'Red to White', 'Blue-White-Red'}, ...
        'FontSize', 10, 'Tag', 'ColormapPopup', 'Callback', @(src, event) colormapSelectionChanged(src, event, susceptibilityPanel));

    %% Data Simulation Tab
    simulationPanel = uipanel('Parent', simulationTab, 'Position', [0.03, 0.88, 0.95, 0.1], ...
                              'Title', 'Simulation Inputs');

    % TR input
    uicontrol('Parent', simulationPanel, 'Style', 'text', 'Position', [10 20 55 20], ...
              'String', 'TR (ms):', 'FontSize', 10);
    trInput = uicontrol('Parent', simulationPanel, 'Style', 'edit', 'Position', [65 20 60 20], 'String', '50');

    % FA input
    uicontrol('Parent', simulationPanel, 'Style', 'text', 'Position', [125 20 100 20], ...
              'String', 'FA (degrees):', 'FontSize', 10);
    faInput = uicontrol('Parent', simulationPanel, 'Style', 'edit', 'Position', [225 20 60 20], 'String', '15');

    % Number of Echoes input
    uicontrol('Parent', simulationPanel, 'Style', 'text', 'Position', [285 20 130 20], ...
              'String', 'Number of Echoes:', 'FontSize', 10);  
    numEchoesInput = uicontrol('Parent', simulationPanel, 'Style', 'edit', 'Position', [415 20 60 20], 'String', '4', ...
                               'Callback', @numEchoesChanged);

    % TE selection droplist and value input
    uicontrol('Parent', simulationPanel, 'Style', 'text', 'Position', [485 20 55 20], ...
              'String', 'TE (ms):', 'FontSize', 10);
    teDropdown = uicontrol('Parent', simulationPanel, 'Style', 'popupmenu', ...
                           'Position', [540 20 70 20], 'String', {'TE1', 'TE2', 'TE3', 'TE4'}, ...
                           'Callback', @teDropdownCallback);
    teValueInput = uicontrol('Parent', simulationPanel, 'Style', 'edit', ...
                             'Position', [610 20 60 20], 'String', '');

    % Initialize TE values storage
    teValues = repmat({''}, 1, 4);  % Initialize for 4 TEs
    currentTEIndex = 1; % Track the currently selected TE

    % Simulate Button
    simulateButton = uicontrol('Parent', simulationPanel, 'Style', 'pushbutton', 'Position', [700, 15, 120, 30], ...
                               'String', 'Simulate', 'FontSize', 12, 'FontWeight', 'bold', 'Callback', @runSimulation);

    % Callback Function for the Simulate button
    function runSimulation(~, ~)
        % Start parallel pool if not already running
        if isempty(gcp('nocreate'))
            parpool('local', 8);
        end
        % Collect the values of TE from the workspace
        numEchoes = str2double(get(numEchoesInput, 'String'));
         % Save the number of echoes in the base workspace
        assignin('base', 'numEchoes', numEchoes);
        teValues = zeros(1, numEchoes);  % Initialize an array to store TE values
        
        % Retrieve the TE values saved in the base workspace
        for i = 1:numEchoes
            teVarName = sprintf('TE%d', i);  % TE variable name (TE1, TE2, ...)
            if evalin('base', ['exist(''', teVarName, ''', ''var'')'])  % Check if TE exists in the workspace
                teValues(i) = evalin('base', teVarName);  % Get TE value from workspace
            else
                errordlg(['TE' num2str(i) ' is not set. Please enter all TE values.']);
                return;
            end
        end
        
        % Save the TE values in the base workspace
        for i = 1:numEchoes
            assignin('base', sprintf('TE%d', i), teValues(i));
        end
        
        % Other inputs like TR, FA can be collected here
        trValue = str2double(get(trInput, 'String'));
        faValue = str2double(get(faInput, 'String'));
        
        % Make sure TR and FA are valid numbers
        if isnan(trValue) || isnan(faValue)
            errordlg('Invalid input for TR or FA. Please check your values.');
            return;
        end
        
        % Save TR and FA values to the base workspace
        assignin('base', 'TR', trValue);
        assignin('base', 'FlipAngle', faValue);
    
        % Extract B0 field strength from the dropdown
        B0_value = get(fieldStrengthPopup, 'Value');  % Get the index of the selected value
        if B0_value == 1
            B0_numeric = 3;  % '3T' corresponds to 3
        else
            B0_numeric = 7;  % '7T' corresponds to 7
        end
        assignin('base', 'B0', B0_numeric);  % Save the numeric B0 value to the workspace
    
        % Extract anisotropy value from checkbox (1 if checked, 0 if not)
        anisotropy = get(anisotropyCheckbox, 'Value');
        assignin('base', 'anisotropy', anisotropy);
    
        % Create the custom progress dialog
        progressFig = createProgressDialog('Simulation in progress...');
        
        % Ensure the progress dialog is drawn
        drawnow;
        
        try
            % Run the DataSimulation_GUI function
            DataSimulation_GUI();
            
            % Close the progress dialog when done
            if ishandle(progressFig)
                close(progressFig);
            end
            
            % Display success message
            msgbox('Simulation completed successfully!');
            
            % After the simulation completes, display the Magnitude and Phase images
            displaySimulationResults(); 
        catch ME
            % Close the progress dialog if an error occurs
            if ishandle(progressFig)
                close(progressFig);
            end
            
            % Display error message
            errordlg(['Error running simulation: ' ME.message], 'Simulation Error');
        end

        % Update the echo selection dropdown based on the number of echoes
        echoOptions = arrayfun(@(x) sprintf('Echo %d', x), 1:numEchoes, 'UniformOutput', false);
        set(echoDropdown, 'String', echoOptions, 'Value', 1);  % Reset to first echo
    end

    % Results panel for magnitude and phase images
    resultsPanel = uipanel('Parent', simulationTab, 'Position', [0.03, 0.05, 0.95, 0.83], 'Title', 'Results');
   
    % Add dropdown for selecting Echo number
    uicontrol('Parent', resultsPanel, 'Style', 'text', 'Position', [100, 50, 100, 20], ...
          'String', 'Select Echo:', 'FontSize', 10);
    echoDropdown = uicontrol('Parent', resultsPanel, 'Style', 'popupmenu', 'Position', [220, 50, 100, 20], ...
                         'String', {'Echo 1', 'Echo 2', 'Echo 3', 'Echo 4'}, 'Callback', @refreshEchoDisplay);

    % Add Colormap Selection to resultsPanel
    uicontrol('Parent', resultsPanel, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.01, 0.95, 0.1, 0.04], 'String', 'Colormap:', 'FontSize', 10);

    colormapPopupResults = uicontrol('Parent', resultsPanel, 'Style', 'popupmenu', 'Units', 'normalized', ...
        'Position', [0.12, 0.94, 0.2, 0.05], 'String', {'Jet', 'Grayscale', 'Blue to White', 'Red to White', 'Blue-White-Red'}, ...
        'FontSize', 10, 'Tag', 'ColormapPopup', 'Callback', @(src, event) colormapSelectionChanged(src, event, resultsPanel));

    %% Callback Functions
    function generateMaps(~, ~)
        % Start parallel pool if not already running
        if isempty(gcp('nocreate'))
            parpool('local', 8);
        end
        % Clear any previous axes to remove placeholders
        delete(r2Axes); delete(r1Axes); delete(r2StarAxes);
        delete(chiPlusAxes); delete(chiMinusAxes); delete(chiTotalAxes);
        
        % Create the custom progress dialog
        progressFig = createProgressDialog('Map generation in progress...');
        
        % Ensure the progress dialog is drawn
        drawnow;
        
        try
            % Run the PhantomCreation function
            evalin('base', 'PhantomCreation();');
            
            % Close the progress dialog when done
            if ishandle(progressFig)
                close(progressFig);
            end
            
            % Display success message
            msgbox('Map generation completed successfully!');
            
            % Call the function to display images once map generation is done
            displayGeneratedMaps();
        catch ME
            % Close the progress dialog if an error occurs
            if ishandle(progressFig)
                close(progressFig);
            end
            
            % Display error message
            errordlg(['Error: ' ME.message], 'Error');
        end
    end

    function displayGeneratedMaps()
        % Get the current selection of field strength and anisotropy
        fieldStrength = get(fieldStrengthPopup, 'Value');
        isAnisotropy = get(anisotropyCheckbox, 'Value');

        % Determine file paths based on field strength and anisotropy
        if fieldStrength == 1  % 3T
            if isAnisotropy
                R2File = 'data/maps/R2_3T.nii.gz';
                R1File = 'data/maps/R1_3T.nii.gz';
                R2StarFile = 'data/maps/R2star_3T.nii.gz';
                ChiPlusFile = 'data/chimodel/chi_positive.nii.gz';
                ChiMinusFile = 'data/chimodel/chi_negative_with_anisotropy.nii.gz';
                ChiTotalFile = 'data/chimodel/chi_with_anisotropy.nii.gz';
            else
                R2File = 'data/maps/R2_3T.nii.gz';
                R1File = 'data/maps/R1_3T.nii.gz';
                R2StarFile = 'data/maps/R2star_3T.nii.gz';
                ChiPlusFile = 'data/chimodel/chi_positive.nii.gz';
                ChiMinusFile = 'data/chimodel/chi_negative.nii.gz';
                ChiTotalFile = 'data/chimodel/chi.nii.gz';
            end
        else  % 7T
            R2File = 'data/maps/R2.nii.gz';
            R1File = 'data/maps/R1.nii.gz';
            R2StarFile = 'data/maps/R2star.nii.gz';
            if isAnisotropy
                ChiPlusFile = 'data/chimodel/chi_positive.nii.gz';
                ChiMinusFile = 'data/chimodel/chi_negative_with_anisotropy.nii.gz';
                ChiTotalFile = 'data/chimodel/chi_with_anisotropy.nii.gz';
            else
                ChiPlusFile = 'data/chimodel/chi_positive.nii.gz';
                ChiMinusFile = 'data/chimodel/chi_negative.nii.gz';
                ChiTotalFile = 'data/chimodel/chi.nii.gz';
            end
        end

        % Define control positions for each map
        % Adjust these positions to fine-tune the placement of the controls
        controlPositions_R2.maxLabelPosition = [0.77, 0.85, 0.1, 0.03];
        controlPositions_R2.maxEditPosition  = [0.87, 0.85, 0.1, 0.03];
        controlPositions_R2.minLabelPosition = [0.77, 0.75, 0.1, 0.03];
        controlPositions_R2.minEditPosition  = [0.87, 0.75, 0.1, 0.03];

        controlPositions_R1.maxLabelPosition = [0.77, 0.55, 0.1, 0.03];
        controlPositions_R1.maxEditPosition  = [0.87, 0.55, 0.1, 0.03];
        controlPositions_R1.minLabelPosition = [0.77, 0.45, 0.1, 0.03];
        controlPositions_R1.minEditPosition  = [0.87, 0.45, 0.1, 0.03];

        controlPositions_R2Star.maxLabelPosition = [0.77, 0.25, 0.1, 0.03];
        controlPositions_R2Star.maxEditPosition  = [0.87, 0.25, 0.1, 0.03];
        controlPositions_R2Star.minLabelPosition = [0.77, 0.15, 0.1, 0.03];
        controlPositions_R2Star.minEditPosition  = [0.87, 0.15, 0.1, 0.03];

        controlPositions_ChiPlus.maxLabelPosition = [0.77, 0.85, 0.1, 0.03];
        controlPositions_ChiPlus.maxEditPosition  = [0.87, 0.85, 0.1, 0.03];
        controlPositions_ChiPlus.minLabelPosition = [0.77, 0.75, 0.1, 0.03];
        controlPositions_ChiPlus.minEditPosition  = [0.87, 0.75, 0.1, 0.03];

        controlPositions_ChiMinus.maxLabelPosition = [0.77, 0.55, 0.1, 0.03];
        controlPositions_ChiMinus.maxEditPosition  = [0.87, 0.55, 0.1, 0.03];
        controlPositions_ChiMinus.minLabelPosition = [0.77, 0.45, 0.1, 0.03];
        controlPositions_ChiMinus.minEditPosition  = [0.87, 0.45, 0.1, 0.03];

        controlPositions_ChiTotal.maxLabelPosition = [0.77, 0.25, 0.1, 0.03];
        controlPositions_ChiTotal.maxEditPosition  = [0.87, 0.25, 0.1, 0.03];
        controlPositions_ChiTotal.minLabelPosition = [0.77, 0.15, 0.1, 0.03];
        controlPositions_ChiTotal.minEditPosition  = [0.87, 0.15, 0.1, 0.03];

        % Display images with sliders and custom control positions
        r2Axes = displayImageWithSlider(fieldStrengthPanel, [0.25, 0.7, 0.65, 0.25], R2File, 'R2 (1/s)', [0, 30], 'r2Slider', false, controlPositions_R2);
        r1Axes = displayImageWithSlider(fieldStrengthPanel, [0.25, 0.4, 0.65, 0.25], R1File, 'R1 (1/s)', [0, 1.5], 'r1Slider', false, controlPositions_R1);
        r2StarAxes = displayImageWithSlider(fieldStrengthPanel, [0.25, 0.1, 0.65, 0.25], R2StarFile, 'R2* (1/s)', [0, 100], 'r2StarSlider', false, controlPositions_R2Star);

        chiPlusAxes = displayImageWithSlider(susceptibilityPanel, [0.25, 0.7, 0.65, 0.25], ChiPlusFile, 'χ+ (ppm)', [0, 0.1], 'chiPlusSlider', false, controlPositions_ChiPlus);
        chiMinusAxes = displayImageWithSlider(susceptibilityPanel, [0.25, 0.4, 0.65, 0.25], ChiMinusFile, 'χ- (ppm)', [-0.1, 0], 'chiMinusSlider', false, controlPositions_ChiMinus);
        chiTotalAxes = displayImageWithSlider(susceptibilityPanel, [0.25, 0.1, 0.65, 0.25], ChiTotalFile, 'χ total (ppm)', [-0.1, 0.1], 'chiTotalSlider', false, controlPositions_ChiTotal);

        % Implement the dashed red box and message when fieldStrength == 2 (7T)
        if fieldStrength == 2  % 7T
            % For r1Axes
            axes(r1Axes);
            hold on;
            xlim = get(r1Axes, 'XLim');
            ylim = get(r1Axes, 'YLim');
            rect1 = rectangle('Position', [xlim(1), ylim(1), diff(xlim), diff(ylim)], ...
                'EdgeColor', 'red', 'LineStyle', '--', 'LineWidth', 2);
            hold off;

            % For r2StarAxes
            axes(r2StarAxes);
            hold on;
            xlim = get(r2StarAxes, 'XLim');
            ylim = get(r2StarAxes, 'YLim');
            rect2 = rectangle('Position', [xlim(1), ylim(1), diff(xlim), diff(ylim)], ...
                'EdgeColor', 'red', 'LineStyle', '--', 'LineWidth', 2);
            hold off;

            % Create text message in fieldStrengthPanel
            msgHandle = uicontrol('Parent', fieldStrengthPanel, 'Style', 'text', 'Units', 'normalized', ...
                'Position', [0.20, 0.02, 0.7, 0.05], 'String', '7T R1 and R2* from the QSM validation phantom', ...
                'ForegroundColor', 'red', 'BackgroundColor', get(fieldStrengthPanel, 'BackgroundColor'), ...
                'HorizontalAlignment', 'left', 'FontSize', 10);

            % Store handles to delete later if needed
            setappdata(hFig, 'rect1', rect1);
            setappdata(hFig, 'rect2', rect2);
            setappdata(hFig, 'msgHandle', msgHandle);
        else
            % Delete any existing rectangles and message
            if isappdata(hFig, 'rect1')
                delete(getappdata(hFig, 'rect1'));
                rmappdata(hFig, 'rect1');
            end
            if isappdata(hFig, 'rect2')
                delete(getappdata(hFig, 'rect2'));
                rmappdata(hFig, 'rect2');
            end
            if isappdata(hFig, 'msgHandle')
                delete(getappdata(hFig, 'msgHandle'));
                rmappdata(hFig, 'msgHandle');
            end
        end
    end

    function teDropdownCallback(~, ~)
        % Save the current value in the text box for the current TE
        currentTEValue = str2double(get(teValueInput, 'String'));
        
        % Make sure the value is valid before saving
        if isnan(currentTEValue)
            errordlg('Please enter a valid numeric value for TE.');
            return;
        end
        
        % Save the value to the base workspace (e.g., TE1, TE2, etc.)
        teVarName = sprintf('TE%d', currentTEIndex);  % Current TE (TE1, TE2, ...)
        assignin('base', teVarName, currentTEValue);
        
        % Update the current TE index based on the dropdown selection
        currentTEIndex = get(teDropdown, 'Value');
        
        % Clear the text box and show the value for the new TE selection
        if evalin('base', ['exist(''', sprintf('TE%d', currentTEIndex), ''', ''var'')'])
            set(teValueInput, 'String', evalin('base', sprintf('TE%d', currentTEIndex)));
        else
            set(teValueInput, 'String', '');  % If no value has been entered yet
        end
    end

    function numEchoesChanged(~, ~)
        % Get the new number of echoes
        numEchoes = str2double(get(numEchoesInput, 'String'));
        if isnan(numEchoes) || numEchoes <= 0
            numEchoes = 1; % Fallback to 1 if invalid
        end
        
        % Create the new dropdown options for TE values based on the number of echoes
        teOptions = arrayfun(@(x) sprintf('TE%d', x), 1:numEchoes, 'UniformOutput', false);

        % Update the existing dropdown (teDropdown) with new options
        set(teDropdown, 'String', teOptions, 'Value', 1);  % Reset to first TE
        
        % Reset the text field for entering TE values
        set(teValueInput, 'String', '');
    end

    % Function to display Magnitude and Phase images
    function displaySimulationResults()
        % Define paths for Magnitude and Phase images
        magnitudePath = 'Simdata/Simulation_Results/SimulatedHR/Magnitude_data.nii.gz';
        phasePath = 'Simdata/Simulation_Results/SimulatedHR/Phase_data.nii.gz';

        % Check if files exist
        if ~isfile(magnitudePath) || ~isfile(phasePath)
            errordlg('Magnitude or Phase image file not found!');
            return;
        end

        % Define control positions for Magnitude image
        controlPositions_Mag.maxLabelPosition = [0.15, 0.9, 0.05, 0.03];
        controlPositions_Mag.maxEditPosition  = [0.18, 0.9, 0.05, 0.03];
        controlPositions_Mag.minLabelPosition = [0.28, 0.9, 0.05, 0.03];
        controlPositions_Mag.minEditPosition  = [0.31, 0.9, 0.05, 0.03];

        % Create axes for Magnitude (placed to the left)
        magnitudeAxes = displayImageWithSlider(resultsPanel, [0.1, 0.1, 0.4, 0.8], ...
            magnitudePath, 'Magnitude (a.u.)', [0, 400], 'magnitudeSlider', true, controlPositions_Mag);

        % Define control positions for Phase image
        controlPositions_Phase.maxLabelPosition = [0.6, 0.9, 0.05, 0.03];
        controlPositions_Phase.maxEditPosition  = [0.63, 0.9, 0.05, 0.03];
        controlPositions_Phase.minLabelPosition = [0.73, 0.9, 0.05, 0.03];
        controlPositions_Phase.minEditPosition  = [0.78, 0.9, 0.05, 0.03];

        % Create axes for Phase (placed to the right)
        phaseAxes = displayImageWithSlider(resultsPanel, [0.55, 0.1, 0.4, 0.8], ...
            phasePath, 'Phase (radians)', [-pi, pi], 'phaseSlider', true, controlPositions_Phase);

        % Save the axes handles to the base workspace
        assignin('base', 'magnitudeAxes', magnitudeAxes);
        assignin('base', 'phaseAxes', phaseAxes);
    end

    % Define the refreshEchoDisplay function
    function refreshEchoDisplay(src, ~)
        % Get the currently selected echo from the dropdown
        selectedEcho = get(src, 'Value');
        
        % Update the displayed images for the selected echo
        updateDisplayedEcho(selectedEcho);
    end

    % Function to update the images for the selected echo's Magnitude and Phase
    function updateDisplayedEcho(selectedEcho)
        % Retrieve the axes handles from the base workspace
        magnitudeAxes = evalin('base', 'magnitudeAxes');
        phaseAxes = evalin('base', 'phaseAxes');

        % Load the Magnitude and Phase 4D data
        magnitudePath = 'Simdata/Simulation_Results/SimulatedHR/Magnitude_data.nii.gz';
        phasePath = 'Simdata/Simulation_Results/SimulatedHR/Phase_data.nii.gz';

        % Load the 4D NIfTI images
        magnitudeNii = load_untouch_nii(magnitudePath);
        phaseNii = load_untouch_nii(phasePath);

        % Extract the 3D volume for the selected echo
        magnitudeData = magnitudeNii.img(:, :, :, selectedEcho);
        phaseData = phaseNii.img(:, :, :, selectedEcho);

        % Update the image data in the sliders' UserData
        % For Magnitude
        magnitudeSlider = findobj(resultsPanel, 'Tag', 'magnitudeSlider');
        if isempty(magnitudeSlider)
            error('Magnitude slider not found.');
        end
        magnitudeSlider.UserData.imgData = magnitudeData;

        % For Phase
        phaseSlider = findobj(resultsPanel, 'Tag', 'phaseSlider');
        if isempty(phaseSlider)
            error('Phase slider not found.');
        end
        phaseSlider.UserData.imgData = phaseData;

        % Reset the sliders to the middle slice
        midSlice = round(size(magnitudeData, 3) / 2);
        set(magnitudeSlider, 'Min', 1, 'Max', size(magnitudeData, 3), 'Value', midSlice);
        set(phaseSlider, 'Min', 1, 'Max', size(phaseData, 3), 'Value', midSlice);

        % Update the displayed images
        % For Magnitude
        rotatedImg = rot90(magnitudeData(:, :, midSlice), 1);
        set(magnitudeSlider.UserData.imHandle, 'CData', rotatedImg);

        % For Phase
        rotatedImg = rot90(phaseData(:, :, midSlice), 1);
        set(phaseSlider.UserData.imHandle, 'CData', rotatedImg);
    end

    % Modified displayImageWithSlider function to accept custom control positions
    function axesHandle = displayImageWithSlider(parentPanel, axesPosition, filepath, mapTitle, contrastLimits, sliderTag, isSimulation, controlPositions)
        if exist(filepath, 'file')
            nii = load_nii(filepath);
            imgData = nii.img;  % 3D image data
            numSlices = size(imgData, 3);
            midSlice = round(numSlices / 2);

            % Create axes
            axesHandle = axes('Parent', parentPanel, 'Position', axesPosition);

            % Display the initial image (mid slice)
            rotatedImg = rot90(imgData(:, :, midSlice), 1);  % Rotate 90 degrees counterclockwise
            imHandle = imshow(rotatedImg, contrastLimits, 'Parent', axesHandle);  % Display the image in axes
            hold(axesHandle, 'on');  % Allow overlaying graphics
            title(axesHandle, mapTitle);  % Set the title for the map
            cbarHandle = colorbar('peer', axesHandle);  % Add colorbar and get handle

            % Create a slider
            sliderWidth = 0.03;  % Width of the slider
            sliderLeft = axesPosition(1) - sliderWidth - 0.02;  % Left position of the slider
            sliderPosition = [sliderLeft, axesPosition(2), sliderWidth, axesPosition(4)];  % Same height as axes

            sliderHandle = uicontrol('Parent', parentPanel, 'Style', 'slider', ...
                'Units', 'normalized', 'Position', sliderPosition, ...
                'Min', 1, 'Max', numSlices, 'Value', midSlice, ...
                'SliderStep', [1/(numSlices-1), 10/(numSlices-1)], ...
                'Callback', @sliceSliderCallback, ...
                'Tag', sliderTag);  % Assign the slider tag here

            % Store the image data and handles in the slider's UserData
            sliderHandle.UserData.imgData = imgData;
            sliderHandle.UserData.imHandle = imHandle;
            sliderHandle.UserData.axesHandle = axesHandle;
            sliderHandle.UserData.contrastLimits = contrastLimits;

            % Check if controlPositions is provided
            if nargin >= 8 && ~isempty(controlPositions)
                % Use the positions from controlPositions
                maxLabelPosition = controlPositions.maxLabelPosition;
                maxEditPosition  = controlPositions.maxEditPosition;
                minLabelPosition = controlPositions.minLabelPosition;
                minEditPosition  = controlPositions.minEditPosition;

                % Create labels and edit boxes
                maxLabel = uicontrol('Parent', parentPanel, 'Style', 'text', 'Units', 'normalized', ...
                    'Position', maxLabelPosition, 'String', 'Max:', 'HorizontalAlignment', 'left');
                maxEdit = uicontrol('Parent', parentPanel, 'Style', 'edit', 'Units', 'normalized', ...
                    'Position', maxEditPosition, 'String', num2str(contrastLimits(2)), ...
                    'Callback', @contrastEditCallback);

                minLabel = uicontrol('Parent', parentPanel, 'Style', 'text', 'Units', 'normalized', ...
                    'Position', minLabelPosition, 'String', 'Min:', 'HorizontalAlignment', 'left');
                minEdit = uicontrol('Parent', parentPanel, 'Style', 'edit', 'Units', 'normalized', ...
                    'Position', minEditPosition, 'String', num2str(contrastLimits(1)), ...
                    'Callback', @contrastEditCallback);

                % Store the handles in UserData for access in callbacks
                sliderHandle.UserData.minEdit = minEdit;
                sliderHandle.UserData.maxEdit = maxEdit;
            else
                % Default positions (you can adjust these if needed)
                controlLeft = axesPosition(1) + axesPosition(3) + 0.02;
                controlWidth = 0.05;
                controlHeight = 0.03;
                spacing = 0.005;

                maxLabelPosition = [controlLeft, axesPosition(2) + axesPosition(4) - controlHeight, controlWidth, controlHeight];
                maxEditPosition  = [controlLeft, axesPosition(2) + axesPosition(4) - 2*controlHeight - spacing, controlWidth, controlHeight];

                minLabelPosition = [controlLeft, axesPosition(2) + axesPosition(4) - 3*controlHeight - 2*spacing, controlWidth, controlHeight];
                minEditPosition  = [controlLeft, axesPosition(2) + axesPosition(4) - 4*controlHeight - 3*spacing, controlWidth, controlHeight];

                % Create labels and edit boxes
                maxLabel = uicontrol('Parent', parentPanel, 'Style', 'text', 'Units', 'normalized', ...
                    'Position', maxLabelPosition, 'String', 'Max:', 'HorizontalAlignment', 'left');
                maxEdit = uicontrol('Parent', parentPanel, 'Style', 'edit', 'Units', 'normalized', ...
                    'Position', maxEditPosition, 'String', num2str(contrastLimits(2)), ...
                    'Callback', @contrastEditCallback);

                minLabel = uicontrol('Parent', parentPanel, 'Style', 'text', 'Units', 'normalized', ...
                    'Position', minLabelPosition, 'String', 'Min:', 'HorizontalAlignment', 'left');
                minEdit = uicontrol('Parent', parentPanel, 'Style', 'edit', 'Units', 'normalized', ...
                    'Position', minEditPosition, 'String', num2str(contrastLimits(1)), ...
                    'Callback', @contrastEditCallback);

                % Store the handles in UserData for access in callbacks
                sliderHandle.UserData.minEdit = minEdit;
                sliderHandle.UserData.maxEdit = maxEdit;
            end

        else
            errordlg(['File not found: ' filepath]);
        end

        % Callback function for the slider
        function sliceSliderCallback(hObject, ~)
            % Get the new slice number
            sliceNum = round(get(hObject, 'Value'));

            % Get the stored image data and image handle
            imgData = hObject.UserData.imgData;
            imHandle = hObject.UserData.imHandle;
            axesHandle = hObject.UserData.axesHandle;

            % Get current contrast limits
            contrastLimits = hObject.UserData.contrastLimits;

            % Update the displayed image
            rotatedImg = rot90(imgData(:, :, sliceNum), 1);
            set(imHandle, 'CData', rotatedImg);

            % Update the contrast limits
            set(axesHandle, 'CLim', contrastLimits);
        end

        % Callback function for contrast edits
        function contrastEditCallback(~, ~)
            % Get the min and max values from the edit boxes
            minVal = str2double(get(sliderHandle.UserData.minEdit, 'String'));
            maxVal = str2double(get(sliderHandle.UserData.maxEdit, 'String'));

            % Validate the inputs
            if isnan(minVal) || isnan(maxVal) || minVal >= maxVal
                errordlg('Invalid contrast limits');
                return;
            end

            % Update the contrast limits in UserData
            sliderHandle.UserData.contrastLimits = [minVal, maxVal];

            % Update the displayed image
            sliceNum = round(get(sliderHandle, 'Value'));
            imgData = sliderHandle.UserData.imgData;
            imHandle = sliderHandle.UserData.imHandle;
            axesHandle = sliderHandle.UserData.axesHandle;

            rotatedImg = rot90(imgData(:, :, sliceNum), 1);
            set(imHandle, 'CData', rotatedImg);

            % Update the contrast limits
            set(axesHandle, 'CLim', [minVal, maxVal]);
        end
    end

    function cmap = getColormapByName(name)
        switch name
            case 'Jet'
                cmap = jet(256);
            case 'Grayscale'
                cmap = gray(256);
            case 'Blue to White'
                cmap = [linspace(0, 0, 256)', linspace(0, 1, 256)', ones(256, 1)];
            case 'Red to White'
                cmap = [ones(256, 1), linspace(0, 1, 256)', linspace(0, 0, 256)'];
            case 'Blue-White-Red'
                half = 128;
                blueToWhite = [linspace(0, 1, half)', linspace(0, 1, half)', ones(half, 1)];
                whiteToRed = [ones(half, 1), linspace(1, 0, half)', linspace(1, 0, half)'];
                cmap = [blueToWhite; whiteToRed];
            otherwise
                cmap = gray(256);  % Default colormap
        end
    end

    function colormapSelectionChanged(src, ~, panel)
        % Get selected colormap
        colormapOptions = get(src, 'String');
        selectedIndex = get(src, 'Value');
        selectedColormapName = colormapOptions{selectedIndex};
        selectedColormap = getColormapByName(selectedColormapName);

        % Apply the colormap to all axes in the panel
        axesHandles = findall(panel, 'Type', 'axes');
        for i = 1:length(axesHandles)
            colormap(axesHandles(i), selectedColormap);
        end
    end

    function progressFig = createProgressDialog(message)
        % Create a figure for the progress dialog
        progressFig = figure('MenuBar', 'none', 'ToolBar', 'none', 'NumberTitle', 'off', ...
            'Name', 'Progress', 'Resize', 'off', 'WindowStyle', 'modal');  % Changed to 'modal' to keep focus
        
        % Remove figure decorations (optional)
        set(progressFig, 'DockControls', 'off');
        
        % Center the figure on the screen
        movegui(progressFig, 'center');
        
        % Set the figure size
        set(progressFig, 'Position', [500, 500, 300, 100]);  % Adjust position and size as needed
        
        % Add a text message
        uicontrol('Parent', progressFig, 'Style', 'text', 'String', message, 'FontSize', 12, ...
            'Units', 'normalized', 'Position', [0.1, 0.4, 0.8, 0.2], 'HorizontalAlignment', 'center');
        
        drawnow;  % Ensure the figure renders immediately
    end
end
