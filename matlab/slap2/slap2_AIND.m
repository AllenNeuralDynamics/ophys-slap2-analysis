function slap2_AIND
    %starts or resumes a SLAP2 session and populates AIND-specific metadata

    localDataRoot = 'F:\slap2Data\';
    schemaUrl = 'https://neural-dynamics-dev.uw.r.appspot.com';

    answer = questdlg('AIND SLAP2 Startup','','New Session','Resume Session','Cancel','New Session');
    switch answer
        case 'Cancel'
            return
        case 'New Session'

            lastSessFn = [localDataRoot 'lastSessionMetadata.mat'];
            try
                lastSessionMetaData = load(lastSessFn);
            catch
                lastSessionMetaData = {'000000', 'default', '', '', '', '', ''};
            end

            %load the current rig description
            rigFn = [localDataRoot 'currentRigDescription.json'];
            if ~exist(rigFn, 'file')
                web(schemaUrl)
                rigFn = uigetfile([localDataRoot '*.json'], 'Select a Rig description');
            end
            if isempty(rigFn)
                return
            end

            answer = inputdlg({'Mouse ID; 000000 if in vitro', 'Session Description', ...
                               'FOV_number', 'Neuron_number', 'FOV_coordinate_ap', ...
                               'FOV_coordinate_ml', 'Targeted_structure'}, ...
                              'Session Metadata Grabber', ...
                              [1; 1; 1; 1; 1; 1; 1], ...
                              lastSessionMetaData);
            assert(~isempty(answer), 'Metadata cannot be empty!')

            date = datestr(now, 'YYYY-mm-DD_HH-MM-SS');
            sessionFolder = [localDataRoot 'slap2_' answer{1} '_' date];
            mkdir(sessionFolder)
            save([sessionFolder filesep 'DESC_ ' answer{2}], 'answer')
            copyfile(rigFn, [sessionFolder filesep 'rigDescription.json']); %copy the rig description to the data folder
            cd(sessionFolder)

            %pop up the webpage to generate a session .json file
            %this will be replaced by a custom app in future
            web(schemaUrl)

            %MAKE A SESSION JSON FILE
            %load a .json template
            %edit the template

        case 'Resume Session'
            sessionFolder = uigetdir(localDataRoot, 'select local session directory');
            cd(sessionFolder)
    end

    slap2
end