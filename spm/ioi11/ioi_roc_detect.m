function [rocIB out vctDtc] = ioi_roc_detect(varargin)
%vctDtc (e.g. vctB): time series to search
%gld (e.g. stim): time series considered as gold standard
%S: structure of ROC 
%ratioseuil (e.g. S.s_th): threshold, used only after ROC for single point output
vctDtc = varargin{1};
S = varargin{3};
ratioseuil = varargin{4};
if nargin == 4
    gld = varargin{2};
else
    if nargin == 6
        gld_series = varargin{2};
        stim = varargin{5};
        infoGld = varargin{6};
        gld = infoGld.sensi.vct.dtc2;
        mask = stim.indGld;
        calc_corr = 0;
    end   
end

rocIB=[];
ii=1;
found_th = 0; %Boolean to search for point along ROC curve near ratioseuil
S.winpts = round(S.dtwindow*S.freq);
if S.ShiftedPower
    vctDtc = vctDtc-min(vctDtc); %shift to positive values
    vctDtc = vctDtc.^2; %power
end
addn = round(S.shift*S.freq);
if S.shift < 0
    vctDtc = vctDtc(abs(addn):end);
else
    vctDtc = [zeros(addn,1); vctDtc];
end
switch S.ROC_method
    case {2,3}
        %find if vctDtc exceeds threshold in each interval
        max2 = max(vctDtc); 
        min2 = min(vctDtc);    
        ratiobmax = max2;
        ratiobdelta = (max2-min2)/100;
        ratiob = min2-ratiobdelta;        
end

while ratiob<=ratiobmax+ratiobdelta
    ratiob=ratiob+ratiobdelta;
    %[temp.ampDtc temp.indDtc]=findpeaks(vctDtc,'minpeakdistance',dtpeak,'minpeakheight',max1*ratiob);
    
    %find if vctDtc exceeds threshold in each interval
    S.wind2s = S.sk_npts2:S.winpts:S.n; %array of starting points
    S.nwin2 = length(S.wind2s)-1; %number of windows
    tmp2amp = []; %zeros(1,S.nwin2);
    tmp2ind = []; %zeros(1,S.nwin2);
    for j1=1:S.nwin2 %much of this could be taken out of the while loop to save time
        window = S.wind2s(j1):(S.wind2s(j1)+S.winpts-1);
        switch S.ROC_method
            case 2
                [t2amp t2ind] = max(vctDtc(window));
            case 3
                %compute signed area
                t2amp = sum(vctDtc(window))/S.winpts;
                t2ind = S.winpts/2; %just put index in center of window
        end
        if t2amp >= ratiob
            %add to list of detected points
            tmp2ind = [tmp2ind t2ind + (j1-1)*S.winpts];
            tmp2amp = [tmp2amp t2amp];
        end                
    end
    if nargin == 4
        [stimB]=ioi_sensi(gld.indGld,tmp2ind,S);
    else
        if nargin == 6
            [stimB]=ioi_sensi(gld,tmp2ind,S,mask);
        end
    end
    rocIB(ii,1)=0; %stimB.sensi;
    rocIB(ii,2)=0;%1-stimB.speci;
    rocIB(ii,3)=stimB.sensi2;
    rocIB(ii,4)=1-stimB.speci2;
    out.ratio(ii)=ratiob;
    switch S.ROC_method
        case {2,3}
            if rocIB(ii,3) == 0
                ratiob = ratiobmax+2*ratiobdelta; %to exit the while
            end
    end
    if nargin==4 || nargin == 6
        if ~found_th && rocIB(ii,3)<=ratioseuil
            found_th = 1;
            out.ratioSeuil=ratiob;
            if nargin == 4
                out.indGld=gld.indGld;
                out.ampGld=gld.ampGld;
            else 
                if nargin == 6
                    out.indGld=gld;
                    out.ampGld=[];
                end
            end
            out.indDtc = tmp2ind;
            out.ampDtc = tmp2amp;
            out.sensi = stimB;
        end
        if nargin == 6 && ~calc_corr
            calc_corr = 1;
            %find correlation between 2 series outside of mask
            corr_amp = []; corr_ampGld = []; corr_ind = []; corr_indGld = [];
            mask1 = stimB.vct.mask;
            for j1=1:S.nwin2 %much of this could be taken out of the while loop to save time
                if ~mask1(j1)
                    window = S.wind2s(j1):(S.wind2s(j1)+S.winpts-1);
                    switch S.ROC_method
                        case 2
                            [t2amp t2ind] = max(vctDtc(window));
                            [t2ampGld t2indGld] = max(gld_series(window));
                        case 3
                            %compute signed area
                            t2amp = sum(vctDtc(window))/S.winpts;
                            t2ind = S.winpts/2; %just put index in center of window
                            t2ampGld = sum(gld_series(window))/S.winpts;
                            t2indGld = S.winpts/2; 
                    end                
                    corr_amp = [corr_amp t2amp]; corr_ampGld = [corr_ampGld t2ampGld];
                    corr_ind = [corr_ind t2ind]; corr_indGld = [corr_indGld t2indGld];
                end
            end
            out.corr_amp = corr_amp; out.corr_ampGld = corr_ampGld;
            out.corr_ind = corr_ind; out.corr_indGld = corr_indGld;
            out.corr = corr(corr_amp',corr_ampGld');
        end
    end
    ii=ii+1;    
end
%Compute area under ROC curve
meanROC=(rocIB(1:end-1,3)+rocIB(2:end,3))/2;
out.ROCarea = -sum(meanROC.*diff(rocIB(:,4)));