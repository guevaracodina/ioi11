function ioi_HDM_display(M)
plot_algebraic_CMRO2 =M.DO.plot_algebraic_CMRO2;
save_figures = M.DO.save_figures;
HDMdisplay = M.DO.generate_figures;
dir1 = M.dir1;
height_stems = 0.1; %for plotting inputs
Ep = M.Ep;
Cp = M.Cp;
pE = M.pE;
if M.DO.show_normalized_parameters
    Ep = Ep./pE;
    pE = ones(size(pE));
    sd1 = diag(M.pC).^0.5;
    Cp = Cp./(sd1 * sd1');
    norm_str = '(Normalized) ';
else
    norm_str = '';
end
m = M.m;
U = M.U;
M0 = M.M0;
M1 = M.M1;
H1 = M.H1;
K1 = M.K1;
K2 = M.K2;
M.plot_separate_figures = 1;

%-display results
%==========================================================================
t = [1:M.N]*M.dt;
if HDMdisplay || save_figures
    Fhdm    = spm_figure;
    header = get(Fhdm,'Name');
    set(Fhdm,'name','Hemodynamic Modeling')
    
    % display input parameters
    %--------------------------------------------------------------------------
    subplot(2,2,1)
    switch M.O.PhysioModel_Choice
        case 0 %Buxton-Friston
            Np = 6;
        case 1 %Zheng-Mayhew
            Np = 8;
        case 2 %Huppert1
            Np = 12;
        otherwise
            Np = 6;
    end
    P     = Ep(Np:end);
    C     = diag(Cp(Np:end,Np:end));
    %[dummy, j] = max(abs(P)); ?
    j=1;
    %this variable j is the onset type
    spm_barh(P,C)
    axis square
    title({'stimulus efficacy'; 'with 90% confidence intervals'},'FontSize',10)
    switch M.O.PhysioModel_Choice
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
    xlabel(['relative efficacy per event/sec' norm_str])
    %Careful here! this could be confusing later...
    Np = Np -1;
    % display hemodynamic parameters
    %---------------------------------------------------------------------------
    subplot(2,2,3)
    P     = Ep(1:Np);
    pE    = pE(1:Np);
    C     = diag(Cp(1:Np,1:Np));
    spm_barh(P,C,pE)
    title({ [norm_str 'hemodynamic parameters']},'FontSize',10)
    set(gca,'Ytick',[1:3*Np]/3 + 1/2)
    switch M.O.PhysioModel_Choice
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
            set(gca,'YTickLabel',{  'Flow signal decay',...
                sprintf('%0.2f per sec',P(1)),'',...
                'Flow feedback',...
                sprintf('%0.2f per sec',P(2)),'',...
                'CMRO signal decay',...
                sprintf('%0.2f per sec',P(3)),'',...
                'CMRO feedback',...
                sprintf('%0.2f per sec',P(4)),'',...
                'P0',...
                sprintf('%0.0f',P(5)),'',...
                'Windkessel volume',...
                sprintf('%0.2f %s',P(6)*100,'%'),'',...
                'beta',...
                sprintf('%0.2f',P(7)),'',...
                'delta',...
                sprintf('%0.2f',P(8)),'',...
                'HGB',...
                sprintf('%0.2f',P(2)),'',...
                'Ra0',...
                sprintf('%0.2f',P(3)),'',...
                'M0',...
                sprintf('%0.2f',P(4)),''},'FontSize',8)
        otherwise
    end
    
    % get display state kernels (i.e. state dynamics)
    %==========================================================================
    
    % Volterra kernels of states
    %--------------------------------------------------------------------------
    
    subplot(3,2,2)
    leg_str = {'s';'f';'v';'q'};
    switch M.O.PhysioModel_Choice
        case 0 %Buxton-Friston
        case 1 %Zheng-Mayhew
            leg_str = [leg_str; 'w'];
        case 2 %Huppert1
            leg_str = {'sf';'f';'v';'sm';'m';'t';'p'};
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
    if M.O.includeHbR
        modalities = [modalities 'HbR; '];
    end
    if M.O.includeHbO
        modalities = [modalities 'HbO; '];
    end
    if M.O.includeHbT
        modalities = [modalities 'HbT; '];
    end
    if M.O.includeFlow
        modalities = [modalities 'Flow'];
    end
    title({ '1st order kernel';['output: ' modalities]},'FontSize',9)
    ylabel('normalized measure response')
    grid on
    
    subplot(3,2,6)
    axis square
    if isreal(K2)
        imagesc(t,t,K2(:,:,1,j,j))
    else
        disp('K2 is not real, only plotting the real part')
        imagesc(t,t,real(K2(:,:,1,j,j)))
    end
    
    title({ '2nd order kernel';'output: first modality'},'FontSize',9)
    xlabel({'time {seconds} for'; U.name{j}})
    grid on
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %REPEATED CODE, to generate separate figures
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h3 = figure;
    leg_str = {'s';'f';'v';'q'};
    switch M.O.PhysioModel_Choice
        case 0 %Buxton-Friston
        case 1 %Zheng-Mayhew
            leg_str = [leg_str; 'w'];
        case 2 %Huppert1
            leg_str = {'sf';'f';'v';'sm';'m';'t';'p'};
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
    %axis square
    title({['1st order kernels for ' U.name{j}];...
        'state variables'},'FontSize',9)
    ylabel('normalized values')
    %grid on
    legend(leg_str,0);
    xlabel('time (s)')
    
    
    % display output kernels (i.e. BOLD response)
    %--------------------------------------------------------------------------
    h4 = figure;
    tY = K1(:,:,j);
    show_HbO = M.O.includeHbO;
    if show_HbO
        if M.O.includeHbR && M.O.includeHbT
            tY = [tY(:,1:2) -tY(:,1)+tY(:,2) tY(:,3:end)];
        end
    end
    
    plot(t,tY);
    %axis square
    modalities = [];
    leg_str0 = {};
    if M.O.includeHbR
        modalities = [modalities 'HbR; '];
        leg_str0 = [leg_str0 'HbR'];
    end
    if M.O.includeHbT
        modalities = [modalities 'HbT; '];
        leg_str0 = [leg_str0 'HbT'];
    end
    if show_HbO
        if M.O.includeHbR && M.O.includeHbT
            modalities = [modalities 'HbO; '];
            leg_str0 = [leg_str0 'HbO'];
        end
    end
    if M.O.includeFlow
        modalities = [modalities 'Flow'];
        leg_str0 = [leg_str0 'Flow'];
    end
    title({ '1st order kernel';['output: ' modalities]},'FontSize',9)
    ylabel('normalized measure response')
    legend(leg_str0)
    xlabel('time (s)')
    %grid on
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %Add output of prediction of other variables
    
    %Loop over the measures -- this depends on the model
    M.g     = 'ioi_gx_get_states';
    %generate forward model with estimated values
    M.pE = M.Ep;
    M.pC = zeros(size(M.Cp)); %should not be used
    M.l = 4; %?
    cH1 = ioi_direct_model(M,M.U,M.Y);
    
    if plot_algebraic_CMRO2
        %Algebraic relation for m = CMRO2, in arbitrary units
        %m = f * HbR /HbT; assuming gamma_R and gamma_T = 1; f: flow
        tH1 = cH1(:,2) .* cH1(:,4) ./cH1(:,3);
        cH1 = [cH1 tH1];
    end
    cH1 = cH1-ones(size(cH1));
    x    = (1:size(M.Y.y,1))*M.dt;
    u0 = full(M.U.u);
    u0(u0==0) = NaN;
    leg_str = ['Stim'; leg_str];
    for i0=0:size(cH1,2)
        Fhdm_pred{i0+1} = figure;
        stem(x,height_stems*u0,'k'); hold on
        if i0==0
            title('Hemodynamic predictions');
            plot(x,cH1); hold on
            legend(leg_str);
        else
            title(['Hemodynamic prediction: ' leg_str{i0+1}]);
            plot(x,cH1(:,i0)); hold on
        end
        %axis square
        ylabel('normalized values')
        grid on
        xlabel('time (seconds)')
    end
    
    %-Reset title
    %--------------------------------------------------------------------------
    spm('FigName',header);
    spm('Pointer','Arrow')
    spm_input('Thank you',1,'d')
    HDM_str = ['_' M.HDM_str];
    ioi_save_figures(save_figures,HDMdisplay,Fhdm,['HDM' HDM_str '_kernels'],dir1)
    if ~M.O.only_display
        Fsi = spm_figure('GetWin','SI');
        ioi_save_figures(save_figures,1,Fsi,['HDM' HDM_str 'fit'],dir1)
    end
    for i0=0:size(cH1,2)
        if i0 == 0
            strA = ['HDM' HDM_str '_predictions_less1_all'];
        else
            strA = ['HDM' HDM_str '_' leg_str{i0+1} '_predictions_less1'];
        end
        ioi_save_figures(save_figures,HDMdisplay,Fhdm_pred{i0+1},strA,dir1)
    end
    %extract some of the subpblots and put in full plots
    ioi_save_figures(save_figures,HDMdisplay,h3,['HDM' HDM_str '_kernels_large'],dir1)
    ioi_save_figures(save_figures,HDMdisplay,h4,['HDM' HDM_str '_output_large'],dir1)
end