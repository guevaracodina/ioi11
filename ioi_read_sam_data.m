function [Frames, data] = ioi_read_sam_data(ExpeFolder)
% Memory maps data from Sam's setup
% SYNTAX
% [Frames, data] = ioi_read_sam_data(ExpeFolder)
% INPUT
% ExpeFolder        Experiment folder
% OUTPUTS
% Frames            Structure with frame information (lookup tables)
% data              Memory-mapped object containing the images
%_______________________________________________________________________________
% Copyright (C) 2014 LIOM Laboratoire d'Imagerie Optique et Moleculaire
%                    Ecole Polytechnique de Montreal
%_______________________________________________________________________________

try
    %% Standard processing
    NomSeq = [ExpeFolder filesep 'IOI_scan.seq'];
    SizeImage = memmapfile(NomSeq,'Offset',580,'Format','uint32','Repeat',1);
    NombreImage = memmapfile(NomSeq,'Offset',572,'Format','uint32','Repeat',1);
    ImRes_XY = memmapfile(NomSeq,'Offset',548,'Format','uint32','Repeat',2);
    SizeImage = double(SizeImage.Data)/2;
    ImRes_XY = double(ImRes_XY.Data);
    NombreImage = double(NombreImage.Data);
    data = memmapfile(NomSeq,'Offset',1024,'Format',{'uint16', [ImRes_XY(1) ImRes_XY(2)], 'framej';'uint16', SizeImage-ImRes_XY(1)*ImRes_XY(2), 'headerj'},'repeat',NombreImage);
    
    FrameSeq = load([ExpeFolder filesep 'IOI_scaninfo.mat']);
    AnalogIN = load([ExpeFolder filesep 'IOI_aux.mat']);
    
    [Path, FolderData] = fileparts(ExpeFolder);
    
    %% Unnecessary folder /EGC
%     FolderAnalyse = ['Analyse_' FolderData];
    
    %% Images a eliminer...
    figure(1);
    Frames.FrameToSkip_Start = 0;
    Frames.Stim = [];
    Frames.FreqR = [];
    Frames.FreqG = [];
    Frames.FreqY = [];
    Frames.FreqF = [];
    for indF = 1:10
        image = data.Data(indF).framej;
        imagesc(image)
        pause(0.01);
        Answer = questdlg('Is this frame OK?', ...
            'Bad frame manual detection', ...
            'Yes','No','Yes');
        if( strcmp(Answer,'Yes') )
            break;
        else
            Frames.FrameToSkip_Start = Frames.FrameToSkip_Start + 1;
        end
    end
    close(1);
    figure(1);
    Frames.FrameToSkip_End = 0;
    for indF = NombreImage:-1:NombreImage - 10
        image = data.Data(indF).framej;
        imagesc(image)
        pause(0.01);
        Answer = questdlg('Is this frame OK?', ...
            'Bad frame manual detection', ...
            'Yes','No','Yes');
        if( strcmp(Answer,'Yes') )
            break;
        else
            Frames.FrameToSkip_End = Frames.FrameToSkip_End + 1;
        end
    end
    close(1);
    clear Answer indF image;
    
    %% Determination des blocs de Stim
    FreqEch = 10000;
    StimOn = find(AnalogIN.aux(:,1)>7500);
    if( StimOn )
        TempsInterStim = diff(StimOn);
        LimiteStim = find(TempsInterStim > 1*FreqEch);
        BlocStim(1,1) = StimOn(1);
        BlocStim(1,2:length(LimiteStim)+1) = StimOn(LimiteStim+1);
        BlocStim(2,1:length(LimiteStim)) = StimOn(LimiteStim);
        BlocStim(2,end) = StimOn(end);
        
        clear LimiteStim TempsInterStim StimOn;
    else
        BlocStim(1,1) = 0;
        BlocStim(2,1) = 0;
    end

    %Separation des frames stim et non stim
    Frames.Temps = find(diff(AnalogIN.aux(:,2))>7500);
    %Test Nombre d'images
    if( (NombreImage - Frames.FrameToSkip_End - Frames.FrameToSkip_Start) ~= length(Frames.Temps) )
        warndlg('Number of frames differ from camera to NI acquisition files');
        disp(['Nombre images dans séquence: ' num2str((NombreImage - Frames.FrameToSkip_End - Frames.FrameToSkip_Start))]);
        disp(['Nombre images dans I/O:' num2str(length(Frames.Temps))]);
    else
        for indF = 1:length(Frames.Temps)
            if( sum((Frames.Temps(indF) > BlocStim(1,:)).*(Frames.Temps(indF) < BlocStim(2,:))) )
                Frames.Stim(indF) = 1;
            else
                Frames.Stim(indF) = 0;
            end
        end
        clear indF;
        
        %Separation des frames en couleurs (Red, Green, Yellow, Fluo)
        F = find( diff(FrameSeq.Signaux(:,1)) > 0.1);
        F = [1; F];
        Red = find(diff(FrameSeq.Signaux(:,2)) > 0.1);
        if( FrameSeq.Signaux(1,2) > 0 )
            Red = [1; Red];
        end
        Green = find(diff(FrameSeq.Signaux(:,3)) > 0.1);
        if( FrameSeq.Signaux(1,3) > 0 )
            Green = [1; Green];
        end
        Yellow = find(diff(FrameSeq.Signaux(:,4)) > 0.1);
        if( FrameSeq.Signaux(1,4) > 0 )
            Yellow = [1; Yellow];
        end
        Fluo = find(diff(FrameSeq.Signaux(:,5)) > 0.1);
        if( FrameSeq.Signaux(1,5) > 0 )
            Fluo = [1; Fluo];
        end
        
        for indF = 1:length(F)
            if( find(Red == F(indF)) )
                SequenceIllum(indF,1) = 'R';
            elseif( find(Green == F(indF)) )
                SequenceIllum(indF,1) = 'G';
            elseif( find(Yellow == F(indF)) )
                SequenceIllum(indF,1) = 'Y';
            elseif( find(Fluo == F(indF)) )
                SequenceIllum(indF,1) = 'F';
            end
        end
        ind_FR = 1;
        ind_FG = 1;
        ind_FY = 1;
        ind_FF = 1;
        for indF = 1:length(Frames.Temps)
            Frames.Color(indF) = SequenceIllum(mod(indF-1,length(SequenceIllum))+1);
            switch Frames.Color(indF)
                case 'R'
                    Frames.LUT_Red(ind_FR) = indF;
                    ind_FR = ind_FR + 1;
                case 'G'
                    Frames.LUT_Green(ind_FG) = indF;
                    ind_FG = ind_FG + 1;
                case 'Y'
                    Frames.LUT_Yellow(ind_FY) = indF;
                    ind_FY = ind_FY + 1;
                case 'F'
                    Frames.LUT_Fluo(ind_FF) = indF;
                    ind_FF = ind_FF + 1;
            end
        end
        
        Frames.FreqR = length(strfind(SequenceIllum','R'));
        Frames.FreqY = length(strfind(SequenceIllum','Y'));
        Frames.FreqG = length(strfind(SequenceIllum','G'));
        Frames.FreqF = length(strfind(SequenceIllum','F'));
        
        clear Fluo Yellow Green Red SequenceIllum FreqEch F BlocStim FrameSeq indF ind_FF ind_FG ind_FR ind_FY;
        
    end
    %% Resume de la manip:
    disp(' ');
    disp('**** INFO ****');
    disp(' ');
    disp(['Name: ' FolderData ]);
    if( sum(Frames.Stim) > 0 )
        disp('Stimulation?: Yes');
    else
        disp('Stimulation?: No');
    end
    Couleur = '';
    if( Frames.FreqR )
        Couleur = strcat(Couleur, ' Red,');
    end
    if( Frames.FreqG )
        Couleur = strcat(Couleur, ' Green,');
    end
    if( Frames.FreqY )
        Couleur = strcat(Couleur, ' Yellow,');
    end
    if( Frames.FreqF )
        Couleur = strcat(Couleur, ' Fluorescence');
    end
    disp(['Channels used?: ' Couleur]);
    
    %% Memory clean-up
    clear AnalogIN

catch
    %% Missing analog signal
    fprintf('Missing analog signal, retrying...\n');
    % opengl('software');
    clear data SizeImage NombreImage ImRes_XY
    
    %%
    % close all; clear all;
    % ExpeFolder = uigetdir;
    
%     NomSeq = [ExpeFolder filesep 'IOI_scan.seq'];
    SizeImage = memmapfile(NomSeq,'Offset',580,'Format','uint32','Repeat',1);
    NombreImage = memmapfile(NomSeq,'Offset',572,'Format','uint32','Repeat',1);
    ImRes_XY = memmapfile(NomSeq,'Offset',548,'Format','uint32','Repeat',2);
    SizeImage = double(SizeImage.Data)/2;
    ImRes_XY = double(ImRes_XY.Data);
    NombreImage = double(NombreImage.Data);
    data = memmapfile(NomSeq,'Offset',1024,'Format',{'uint16', [ImRes_XY(1) ImRes_XY(2)], 'framej';'uint16', SizeImage-ImRes_XY(1)*ImRes_XY(2), 'headerj'},'repeat',NombreImage);
    
    FrameSeq = load([ExpeFolder filesep 'IOI_scaninfo.mat']);
    AnalogIN = load([ExpeFolder filesep 'IOI_aux.mat']);
    
    [Path, FolderData] = fileparts(ExpeFolder);
%% Unnecessary folder // EGC    
%     FolderAnalyse = ['Analyse_' FolderData];
%     if( ~exist([Path filesep FolderAnalyse],'dir') )
%         mkdir([Path filesep FolderAnalyse]);
%         % else
%         %     rmdir([Path filesep FolderAnalyse],'s' );
%         %     mkdir([Path filesep FolderAnalyse]);
%     end
    
    %% Images a eliminer...
    figure(1);
    Frames.FrameToSkip_Start = 0;
    Frames.Stim = [];
    Frames.FreqR = [];
    Frames.FreqG = [];
    Frames.FreqY = [];
    Frames.FreqF = [];
    for indF = 1:10
        image = data.Data(indF).framej;
        imagesc(image)
        pause(0.01);
        Answer = questdlg('Is this frame OK?', ...
            'Bad frame manual detection', ...
            'Yes','No','Yes');
        if( strcmp(Answer,'Yes') )
            break;
        else
            Frames.FrameToSkip_Start = Frames.FrameToSkip_Start + 1;
        end
    end
    close(1);
    figure(1);
    Frames.FrameToSkip_End = 0;
    for indF = NombreImage:-1:NombreImage - 10
        image = data.Data(indF).framej;
        imagesc(image)
        pause(0.01);
        Answer = questdlg('Is this frame OK?', ...
            'Bad frame manual detection', ...
            'Yes','No','Yes');
        if( strcmp(Answer,'Yes') )
            break;
        else
            Frames.FrameToSkip_End = Frames.FrameToSkip_End + 1;
        end
    end
    close(1);
    clear Answer indF image;
    
    %% Determination des blocs de Stim
    FreqEch = 10000;
    BlocStim(1,1) = 0;
    BlocStim(2,1) = 0;

    %Separation des frames stim et non stim
    Frames.Temps = find(diff(AnalogIN.aux(:,1))>7500);
    %Test Nombre d'images
    if( (NombreImage - Frames.FrameToSkip_End - Frames.FrameToSkip_Start) ~= length(Frames.Temps) )
        warndlg('Number of frames differ from camera to NI acquisition files');
        disp(['Nombre images dans sï¿½quence: ' num2str((NombreImage - Frames.FrameToSkip_End - Frames.FrameToSkip_Start))]);
        disp(['Nombre images dans I/O:' num2str(length(Frames.Temps))]);
    else
        for indF = 1:length(Frames.Temps)
            if( sum((Frames.Temps(indF) > BlocStim(1,:)).*(Frames.Temps(indF) < BlocStim(2,:))) )
                Frames.Stim(indF) = 1;
            else
                Frames.Stim(indF) = 0;
            end
        end
        clear indF;
        
        %Separation des frames en couleurs (Red, Green, Yellow, Fluo)
        F = find( diff(FrameSeq.Signaux(:,1)) > 0.1);
        F = [1; F];
        Red = find(diff(FrameSeq.Signaux(:,2)) > 0.1);
        if( FrameSeq.Signaux(1,2) > 0 )
            Red = [1; Red];
        end
        Green = find(diff(FrameSeq.Signaux(:,3)) > 0.1);
        if( FrameSeq.Signaux(1,3) > 0 )
            Green = [1; Green];
        end
        Yellow = find(diff(FrameSeq.Signaux(:,4)) > 0.1);
        if( FrameSeq.Signaux(1,4) > 0 )
            Yellow = [1; Yellow];
        end
        Fluo = find(diff(FrameSeq.Signaux(:,5)) > 0.1);
        if( FrameSeq.Signaux(1,5) > 0 )
            Fluo = [1; Fluo];
        end
        
        for indF = 1:length(F)
            if( find(Red == F(indF)) )
                SequenceIllum(indF,1) = 'R';
            elseif( find(Green == F(indF)) )
                SequenceIllum(indF,1) = 'G';
            elseif( find(Yellow == F(indF)) )
                SequenceIllum(indF,1) = 'Y';
            elseif( find(Fluo == F(indF)) )
                SequenceIllum(indF,1) = 'F';
            end
        end
        ind_FR = 1;
        ind_FG = 1;
        ind_FY = 1;
        ind_FF = 1;
        for indF = 1:length(Frames.Temps)
            Frames.Color(indF) = SequenceIllum(mod(indF-1,length(SequenceIllum))+1);
            switch Frames.Color(indF)
                case 'R'
                    Frames.LUT_Red(ind_FR) = indF;
                    ind_FR = ind_FR + 1;
                case 'G'
                    Frames.LUT_Green(ind_FG) = indF;
                    ind_FG = ind_FG + 1;
                case 'Y'
                    Frames.LUT_Yellow(ind_FY) = indF;
                    ind_FY = ind_FY + 1;
                case 'F'
                    Frames.LUT_Fluo(ind_FF) = indF;
                    ind_FF = ind_FF + 1;
            end
        end
        
        Frames.FreqR = length(strfind(SequenceIllum','R'));
        Frames.FreqY = length(strfind(SequenceIllum','Y'));
        Frames.FreqG = length(strfind(SequenceIllum','G'));
        Frames.FreqF = length(strfind(SequenceIllum','F'));
        
        clear Fluo Yellow Green Red SequenceIllum FreqEch F BlocStim FrameSeq indF ind_FF ind_FG ind_FR ind_FY;
        
    end
    %% Resume de la manip:
    disp(' ');
    disp('**** INFO ****');
    disp(' ');
    disp(['Name: ' FolderData ]);
    if( sum(Frames.Stim) > 0 )
        disp('Stimulation?: Yes');
    else
        disp('Stimulation?: No');
    end
    Couleur = '';
    if( Frames.FreqR )
        Couleur = strcat(Couleur, ' Red,');
    end
    if( Frames.FreqG )
        Couleur = strcat(Couleur, ' Green,');
    end
    if( Frames.FreqY )
        Couleur = strcat(Couleur, ' Yellow,');
    end
    if( Frames.FreqF )
        Couleur = strcat(Couleur, ' Fluorescence');
    end
    disp(['Channels used?: ' Couleur]);
    
    %% Memory clean-up
    
    clear AnalogIN

    disp(exception.identifier)
    disp(exception.stack(1))
%     out.IOImat{scanIdx} = job.IOImat{scanIdx};
end
% EOF
