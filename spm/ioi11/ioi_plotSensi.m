function ioi_plotSensi(varargin)
% plot sensitivity graph
%speci (e.g. infoSB): 
%vctDtc (e.g. vctB)
%stim
%S
%name (e.g. {'Stim','Bwd'}): names for plot
if nargin == 6 %comparing 2 series
   speci = varargin{1};
   vctDtc = varargin{2};
   stim = varargin{3};
   S = varargin{4};
   names = varargin{5}; 
   O = varargin{6};
else
    if nargin == 8 %comparing 3 series
        infoA = varargin{1};
        infoB = varargin{2};
        gld_series = varargin{3};
        vctDtc = varargin{4};
        stim = varargin{5};
        S = varargin{6};
        names = varargin{7};
        O = varargin{8};
    end
end

if nargin == 6
    gld = speci.sensi.vct.gld2;
    dtc = speci.sensi.vct.dtc2;
    TP = gld~=0 & dtc~=0;
    TN = gld==0 & dtc==0;
    FP = gld==0 & dtc~=0;
    FN = gld~=0 & dtc==0;
else 
    if nargin == 8
        vctA = infoA.sensi.vct;
        vctB = infoB.sensi.vct;
        %check: vctA.gld2 == vctB.gld2, vctA.temps2 == vctB.temps2 (otherwise, inputs are inconsistent)
        bck = vctA.gld2==0; %background
        dtcA = vctA.dtc2;
        dtcB = vctB.dtc2;
        gld = dtcA;
        %consider conjunction with gld==0: what is happening outside of
        %stimulations

        %Considering series A as the truth, to which B is compared
        TP = bck & dtcA~=0 & dtcB~=0;
        TN = bck & dtcA==0 & dtcB==0;
        FP = bck & dtcA==0 & dtcB~=0;
        FN = bck & dtcA~=0 & dtcB==0;
    end
end
%output values:
TP1 = sum(TP); TN1 = sum(TN); FP1 = sum(FP); FN1 = sum(FN);
disp(['TP: ' int2str(TP1) ' FP: ' int2str(FP1) ' FN: ' int2str(FN1) ' TN: ' int2str(TN1)]);
disp(['Sensitivity: ' num2str(TP1/(TP1+FN1),'%1.2f') ' Specificity: ' ...
    num2str(1-FP1/(FP1+TN1),'%1.2f') ' Positive Predictive value: ' num2str(TP1/(TP1+FP1),'%1.2f') ...
    ' Negative Predictive value: ' num2str(TN1/(TN1+FN1),'%1.2f')])
hMain = figure; %set(gcf,'Position',[0 0 900 1300])
%plot data series
offset=S.offset;
temps=(1:length(vctDtc))/S.freq+offset;
h0 = plot(temps,vctDtc);
if nargin == 6
    ylabel([names{2}]);
else
    if nargin == 8
        ylabel([names{2} ' against ' names{1}])
        hold on
        temps3=(1:length(gld_series))/S.freq+offset;
        h00 = plot(temps3,gld_series,'k.-');
    end
end
if S.set_limits
    set(gca,'xlim',S.xlim)
end
xlabel('Time (s)')
hold on
%get vertical limits, to use to set TP, FP, FN, TN labels and stim bars
yy=get(gca,'ylim');
in_vpos = yy(1)+0.05*diff(yy);

%set 4s windows boxes
temps2 = S.dtwindow*(1:length(gld))-S.dtwindow/2; %center of the window
%plot(temps2(TP),in_vpos,'bo',temps2(FP),in_vpos,'rs',temps2(FN),in_vpos,'gd',temps2(TN),in_vpos,'kx')
try h1 = plot(temps2(TP),in_vpos,'bo'); catch, h1 = plot(temps2,in_vpos,''); end
try h2 = plot(temps2(FP),in_vpos,'rs'); catch, h2 = plot(temps2,in_vpos,''); end
try h3 = plot(temps2(FN),in_vpos,'gd'); catch, h3 = plot(temps2,in_vpos,''); end
try h4 = plot(temps2(TN),in_vpos,'kx'); catch, h4 = plot(temps2,in_vpos,''); end
if nargin == 6
    h5 = plot(temps2,speci.ratioSeuil*ones(size(temps2)),'k-'); %
    legend([h0 h1(1) h2(1) h3(1) h4(1) h5(1)],names{2},'TP','FP','FN','TN','Thld')   
else
    if nargin == 8
        h5 = plot(temps2,infoB.ratioSeuil*ones(size(temps2)),'b-'); %
        h55 = plot(temps2,infoA.ratioSeuil*ones(size(temps2)),'k.-'); %
    
        legend([h0 h00 h1(1) h2(1) h3(1) h4(1) h5(1) h55(1)],names{2},names{1},'TP','FP','FN','TN',['Thld ' names{2}],['Thld ' names{1}])
    end
end
hold on
%set(gcf,'tag','doScale')

in_vpos2 = yy(2)-0.05*diff(yy);
t=offset+ stim.indGld/S.freq;
%t=offset + S.stimInd/S.freq; %center of the window - there is a problem in
%the timing of these measures
plot(t,in_vpos2,'k.','markersize',10); %ones(size(t))*
title(O.tit2)
hold off
ioi_save_figures(O.save_figures,O.generate_figures,hMain,O.fname2,O.dir);