function ioi_HDM_display(M)
plot_algebraic_CMRO2 =1;
save_figures = M.save_figures;
HDMdisplay = M.generate_figures;
dir1 = M.dir1;
Ep = M.Ep;
Cp = M.Cp;
m = M.m;
U = M.U;
pE = M.PS.pE;
M0 = M.M0;
M1 = M.M1;
H1 = M.H1;
K1 = M.K1;
K2 = M.K2;

%-display results
%==========================================================================
t       = [1:M.N]*M.dt;
if HDMdisplay || save_figures
    Fhdm    = spm_figure;
    header = get(Fhdm,'Name');
    set(Fhdm,'name','Hemodynamic Modeling')
    
    
    % display input parameters
    %--------------------------------------------------------------------------
    subplot(2,2,1)
    switch M.PS.PhysioModel_Choice
        case 0 %Buxton-Friston
            P     = Ep(7:end);
            C     = diag(Cp(7:end,7:end));
        case 1 %Zheng-Mayhew
            P     = Ep(9:end);                  %MODIFIER
            C     = diag(Cp(9:end,9:end));      %MODIFIER
        case 2 %Huppert1
            P     = Ep(11:end);                 %MODIFIER
            C     = diag(Cp(11:end,11:end));    %MODIFIER
        otherwise
            P     = Ep(7:end);
            C     = diag(Cp(7:end,7:end));
    end
    
    [dummy, j] = max(abs(P));
    spm_barh(P,C)
    axis square
    title({'stimulus efficacy'; 'with 90% confidence intervals'},'FontSize',10)
    switch M.PS.PhysioModel_Choice
        case {0,1} %Buxton-Friston
            set(gca,'Ytick',[1:m],'YTickLabel',U.name,'FontSize',8)
            str = {};
            for i = 1:m
                str{end + 1} = U.name{i};
                str{end + 1} = sprintf('mean = %0.2f',P(i));
                str{end + 1} = '';
            end
            set(gca,'Ytick',[1:m*3]/3 + 1/2,'YTickLabel',str)
        case 2
            m2 = 2*m;
            for m3=1:m
                temp_labels{m3} = ['CMRO2' U.name{m3}];
            end
            labels2 = [U.name temp_labels];
            set(gca,'Ytick',[1:m2],'YTickLabel',labels2,'FontSize',8)
            str = {};
            for i = 1:m2
                str{end + 1} = labels2{i}; %U.name{i};
                str{end + 1} = sprintf('mean = %0.2f',P(i));
                str{end + 1} = '';
            end
            set(gca,'Ytick',[1:m2*3]/3 + 1/2,'YTickLabel',str)
        otherwise
    end
    xlabel('relative efficacy per event/sec')
    
    
    % display hemodynamic parameters
    %---------------------------------------------------------------------------
    subplot(2,2,3)
    switch M.PS.PhysioModel_Choice
        case 0 %Buxton-Friston
            P     = Ep(1:6);
            pE    = pE(1:6);
            C     = diag(Cp(1:6,1:6));
            spm_barh(P,C,pE)
            title({ 'hemodynamic parameters'},'FontSize',10)
            set(gca,'Ytick',[1:18]/3 + 1/2)
            set(gca,'YTickLabel',{  'SIGNAL decay',...
                sprintf('%0.2f per sec',P(1)),'',...
                'FEEDBACK',...
                sprintf('%0.2f per sec',P(2)),'',...
                'TRANSIT TIME',...
                sprintf('%0.2f seconds',P(3)),'',...
                'EXPONENT',...
                sprintf('%0.2f',P(4)),'',...
                'EXTRACTION',...
                sprintf('%0.0f %s',P(5)*100,'%'),'',...
                'log SIGNAL RATIO',...
                sprintf('%0.2f %s',P(6),'%'),''},'FontSize',8)
        case 1 %Zheng-Mayhew
            P     = Ep(1:8);            %MODIFIER
            pE    = pE(1:8);            %MODIFIER
            C     = diag(Cp(1:8,1:8));  %MODIFIER
            spm_barh(P,C,pE)
            title({ 'hemodynamic parameters'},'FontSize',10)
            set(gca,'Ytick',[1:18+6]/3 + 1/2)   %MODIFIER? was 18
            set(gca,'YTickLabel',{  'SIGNAL decay',...
                sprintf('%0.2f per sec',P(1)),'',...
                'FEEDBACK',...
                sprintf('%0.2f per sec',P(2)),'',...
                'TRANSIT TIME',...
                sprintf('%0.2f seconds',P(3)),'',...
                'EXPONENT',...
                sprintf('%0.2f',P(4)),'',...
                'EXTRACTION',...
                sprintf('%0.0f %s',P(5)*100,'%'),'',...
                'VASCULAR TONE',...                                %MODIFIER
                sprintf('%0.2f *10 seconds',P(6)),'',...              %MODIFIER
                'GAIN PARAMETER',...                               %MODIFIER
                sprintf('%0.3f *10 seconds',P(7)),'',...              %MODIFIER
                'log SIGNAL RATIO',...
                sprintf('%0.2f %s',P(8),'%'),''},'FontSize',8)   %MODIFIER
        case 2 %Huppert1
            P     = Ep(1:10);
            pE    = pE(1:10);
            C     = diag(Cp(1:10,1:10));
            spm_barh(P,C,pE)
            title({ 'hemodynamic parameters'},'FontSize',10)
            set(gca,'Ytick',[1:18+12]/3 + 1/2)
            set(gca,'YTickLabel',{  'SIGNAL decay',...
                sprintf('%0.2f per sec',P(1)),'',...
                'FEEDBACK',...
                sprintf('%0.2f per sec',P(2)),'',...
                'TRANSIT TIME',...
                sprintf('%0.2f seconds',P(3)),'',...
                'EXPONENT',...
                sprintf('%0.2f',P(4)),'',...
                'EXTRACTION',...
                sprintf('%0.0f %s',P(5)*100,'%'),'',...
                'log SIGNAL RATIO',...
                sprintf('%0.2f %s',P(6),'%'),''},'',...
                'CMRO2 SIGNAL decay',...
                sprintf('%0.2f per sec',P(7)),'',...
                'CMRO2 FEEDBACK',...
                sprintf('%0.2f per sec',P(8)),'',...
                'Volume fraction',...
                sprintf('%0.2f',P(9)),'',...
                'O2 Diffusion K',...
                sprintf('%0.2f per sec',P(10)),'','FontSize',8)
        otherwise
            P     = Ep(1:6);
            pE    = pE(1:6);
            C     = diag(Cp(1:6,1:6));
            spm_barh(P,C,pE)
            title({ 'hemodynamic parameters'},'FontSize',10)
            set(gca,'Ytick',[1:18]/3 + 1/2)
            set(gca,'YTickLabel',{  'SIGNAL decay',...
                sprintf('%0.2f per sec',P(1)),'',...
                'FEEDBACK',...
                sprintf('%0.2f per sec',P(2)),'',...
                'TRANSIT TIME',...
                sprintf('%0.2f seconds',P(3)),'',...
                'EXPONENT',...
                sprintf('%0.2f',P(4)),'',...
                'EXTRACTION',...
                sprintf('%0.0f %s',P(5)*100,'%'),'',...
                'log SIGNAL RATIO',...
                sprintf('%0.2f %s',P(6),'%'),''},'FontSize',8)
    end
    

% get display state kernels (i.e. state dynamics)
%==========================================================================

% Volterra kernels of states
%--------------------------------------------------------------------------

    subplot(3,2,2)
    if plot_algebraic_CMRO2
        tmp_H1 = exp(H1(:,:,j));
        %Algebraic relation for m = CMRO2, in arbitrary units
        %m = f * HbR /HbT; assuming gamma_R and gamma_T = 1; f: flow
        tmp_ma = tmp_H1(:,2) .* tmp_H1(:,4) ./ tmp_H1(:,3);
        tmp_H1 = [tmp_H1 tmp_ma];
        plot(t,tmp_H1)
        axis square
        title({['1st order kernels for ' U.name{j}];...
            'state variables'},'FontSize',9)
        ylabel('normalized values')
        switch M.PS.PhysioModel_Choice
            case 0 %Buxton-Friston
                legend('s','f','v','q','ma',0);
            case 1 %Zheng-Mayhew
                legend('s','f','v','q','w','ma',0); %MODIFIER
            case 2 %Huppert1
                legend('s','f','v','q','s2','m','Ct','Cv','ma',0)
            otherwise
                legend('s','f','v','q','ma',0);
        end
        grid on
    else
        plot(t,exp(H1(:,:,j)))
        axis square
        title({['1st order kernels for ' U.name{j}];...
            'state variables'},'FontSize',9)
        ylabel('normalized values')
        switch M.PS.PhysioModel_Choice
            case 0 %Buxton-Friston
                legend('s','f','v','q',0);
            case 1 %Zheng-Mayhew
                legend('s','f','v','q','w', 0); %MODIFIER
            case 2 %Huppert1
                legend('s','f','v','q','s2','m','Ct','Cv',0)
            otherwise
                legend('s','f','v','q',0);
        end
        grid on
    end
    
    
    
    % display output kernels (i.e. BOLD response)
    %--------------------------------------------------------------------------
    subplot(3,2,4)
    plot(t,K1(:,:,j))
    axis square
    modalities = [];
    xY = M.PS.xY;
    if xY.includeHbR
        modalities = [modalities 'HbR; '];
    end
    if xY.includeHbT
        modalities = [modalities 'HbT; '];
    end
    if xY.includeFlow
        modalities = [modalities 'Flow'];
    end
    title({ '1st order kernel';['output: ' modalities]},'FontSize',9)
 
    %ylabel('normalized flow signal')
            ylabel('normalized measure response')
    grid on
    
    subplot(3,2,6)
    axis square
    imagesc(t,t,K2(:,:,1,j,j))
    title({ '2nd order kernel';'output: first modality'},'FontSize',9)
    xlabel({'time {seconds} for'; U.name{j}})
    grid on
    
    %-Reset title
    %--------------------------------------------------------------------------
    spm('FigName',header);
    spm('Pointer','Arrow')
    spm_input('Thank you',1,'d')
    if save_figures
        %Save figure
        filen1 = fullfile(dir1,['HDM_kernels.fig']);
        filen2 = fullfile(dir1,['HDM_kernels.tiff']);
        saveas(Fhdm,filen1,'fig');
        print(Fhdm, '-dtiffn', filen2);
        Fsi = spm_figure('GetWin','SI');
        filen2 = fullfile(dir1,['HDM_fit.fig']);
        filen4 = fullfile(dir1,['HDM_fit.tiff']);
        saveas(Fsi,filen2,'fig');
        print(Fsi, '-dtiffn', filen4);
        if HDMdisplay
            try, close(Fhdm); end
        end
    end
end
end