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
            Np = 6;
        case 1 %Zheng-Mayhew
            Np = 8;
        case 2 %Huppert1
            Np = 8;
        otherwise
            Np = 6;
    end
    P     = Ep(Np:end);
    C     = diag(Cp(Np:end,Np:end));
    [dummy, j] = max(abs(P));
    %this variable j is the onset type
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
    
    Np = Np -1;
    % display hemodynamic parameters
    %---------------------------------------------------------------------------
    subplot(2,2,3)
    P     = Ep(1:Np);
    pE    = pE(1:Np);
    C     = diag(Cp(1:Np,1:Np));
    spm_barh(P,C,pE)
    title({ 'hemodynamic parameters'},'FontSize',10)
    set(gca,'Ytick',[1:3*Np]/3 + 1/2)
    switch M.PS.PhysioModel_Choice
        case 0 %Buxton-Friston
            set(gca,'YTickLabel',{  'Signal decay',...
                sprintf('%0.2f per sec',P(1)),'',...
                'Feedback',...
                sprintf('%0.2f per sec',P(2)),'',...
                'Transit time',...
                sprintf('%0.2f seconds',P(3)),'',...
                'Exponent',...
                sprintf('%0.2f',P(4)),'',...
                'Extraction',...
                sprintf('%0.0f %s',P(5)*100,'%'),''},'FontSize',8)
        case 1 %Zheng-Mayhew
            set(gca,'YTickLabel',{  'Signal decay',...
                sprintf('%0.2f per sec',P(1)),'',...
                'Feedback',...
                sprintf('%0.2f per sec',P(2)),'',...
                'Transit time',...
                sprintf('%0.2f seconds',P(3)),'',...
                'Exponent',...
                sprintf('%0.2f',P(4)),'',...
                'Extraction',...
                sprintf('%0.0f %s',P(5)*100,'%'),'',...
                'Vascular tone',...
                sprintf('%0.2f *10 seconds',P(6)),'',...
                'Gain parameter',...
                sprintf('%0.3f *10 seconds',P(7)),''},'FontSize',8)
        case 2 %Huppert1
            %{'ks'    'kr'    'kx'    'km'    'Vw0'    'beta'    'Ra0'    'effCMRO'    'effFlow'};
            set(gca,'YTickLabel',{  'Flow signal decay',...
                sprintf('%0.2f per sec',P(1)),'',...
                'Flow feedback',...
                sprintf('%0.2f per sec',P(2)),'',...
                'CMRO signal decay',...
                sprintf('%0.2f per sec',P(3)),'',...
                'CMRO feedback',...
                sprintf('%0.2f per sec',P(4)),'',...
                'Windkessel volume',...
                sprintf('%0.0f %s',P(5)*100,'%'),'',...
                'beta',...
                sprintf('%0.2f %s',P(6),'%'),'',...
                'Ra0',...
                sprintf('%0.2f',P(7)),''},'FontSize',8)
        otherwise
    end
    
    % get display state kernels (i.e. state dynamics)
    %==========================================================================
    
    % Volterra kernels of states
    %--------------------------------------------------------------------------
    
    subplot(3,2,2)
    leg_str = {'s';'f';'v';'q'};
    switch M.PS.PhysioModel_Choice
        case 0 %Buxton-Friston
        case 1 %Zheng-Mayhew
            leg_str = [leg_str; 'w'];
        case 2 %Huppert1
            leg_str = [leg_str; 's2';'m';'Ct';'Cv'];
        otherwise
    end
    if plot_algebraic_CMRO2
        leg_str = [leg_str; 'ma'];
    end
    if plot_algebraic_CMRO2
        tmp_H1 = exp(H1(:,:,j));
        %Algebraic relation for m = CMRO2, in arbitrary units
        %m = f * HbR /HbT; assuming gamma_R and gamma_T = 1; f: flow
        tmp_ma = tmp_H1(:,2) .* tmp_H1(:,4) ./ tmp_H1(:,3);
        tmp_H1 = [tmp_H1 tmp_ma];
        plot(t,tmp_H1)
    else
        plot(t,exp(H1(:,:,j)))
    end
    axis square
    title({['1st order kernels for ' U.name{j}];...
        'state variables'},'FontSize',9)
    ylabel('normalized values')
    
    grid on
    legend(leg_str,0);
    
    
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
    ylabel('normalized measure response')
    grid on
    
    subplot(3,2,6)
    axis square
    imagesc(t,t,K2(:,:,1,j,j))
    title({ '2nd order kernel';'output: first modality'},'FontSize',9)
    xlabel({'time {seconds} for'; U.name{j}})
    grid on
    
    
    %Add output of prediction of other variables
    Fhdm_fits    = figure; %spm_figure;
    title('Hemodynamic Fits');
    %header = get(Fhdm_fits,'Name');
    %set(Fhdm_fits,'name','Hemodynamic Fits')
    x    = (1:size(M.Y,1))*M.dt;
    u0 = full(M.U.u);
    u0(u0==0) = NaN;
    stem(x,u0,'k'); hold on
    %Loop over the measures -- this depends on the model
    cH1 = [];
    for i0=1:size(H1,2)
        tH1 = conv(full(M.U.u),exp(H1(:,i0,j)));
        cH1 = [cH1 tH1(1:length(U.u))];
    end
    if plot_algebraic_CMRO2
        %Algebraic relation for m = CMRO2, in arbitrary units
        %m = f * HbR /HbT; assuming gamma_R and gamma_T = 1; f: flow
        tH1 = cH1(:,2) .* cH1(:,4) ./cH1(:,3);
        cH1 = [cH1 tH1];
    end
    plot(x,cH1); hold on
    axis square
    title({['1st order kernels for ' U.name{j}];...
        'state variables'},'FontSize',9)
    ylabel('normalized values')
    leg_str = ['Stim'; leg_str];
    legend(leg_str);
    grid on
    xlabel('time (seconds)')
    
    %-Reset title
    %--------------------------------------------------------------------------
    spm('FigName',header);
    spm('Pointer','Arrow')
    spm_input('Thank you',1,'d')
    if save_figures
        HDM_str = ['_' M.HDM_str];
        %Save figure
        filen1 = fullfile(dir1,['HDM_kernels' HDM_str '.fig']);
        filen2 = fullfile(dir1,['HDM_kernels' HDM_str '.tiff']);
        saveas(Fhdm,filen1,'fig');
        print(Fhdm, '-dtiffn', filen2);
        Fsi = spm_figure('GetWin','SI');
        filen2 = fullfile(dir1,['HDM_fit' HDM_str '.fig']);
        filen4 = fullfile(dir1,['HDM_fit' HDM_str '.tiff']);
        saveas(Fsi,filen2,'fig');
        print(Fsi, '-dtiffn', filen4);
        filen5 = fullfile(dir1,['HDM_fits' HDM_str '.fig']);
        filen6 = fullfile(dir1,['HDM_fits' HDM_str '.tiff']);
        saveas(Fhdm_fits,filen5,'fig');
        print(Fhdm_fits, '-dtiffn', filen6);
        if ~HDMdisplay
            try close(Fhdm); end
            try close(Fhdm_fits); end
        end
    end
end
end