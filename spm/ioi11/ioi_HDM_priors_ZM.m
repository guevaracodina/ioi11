function M = ioi_HDM_priors_ZM(M,h)
% returns priors for a hemodynamic dynamic causal model

%
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_hdm_priors.m 2050 2008-09-05 19:15:50Z klaas $


[pE,pC] = nirs_hdm_priors_ZM(M.m,h);

if ~isfield(M,'name')
    M.name={'SIGNAL_decay'  'FEEDBACK', 'TRANSIT_TIME','EXPONENT',...
        'EXTRACTION', 'VASCULAR_TONE' 'GAIN_PARAMETER', 'log_SIGNAL_RATIO' 'effFlow'};
end
namePresent={'SIGNAL_decay'  'FEEDBACK', 'TRANSIT_TIME','EXPONENT',...
    'EXTRACTION', 'VASCULAR_TONE' 'GAIN_PARAMETER', 'log_SIGNAL_RATIO' 'effFlow' };

if M.Model_Choice==3 %decay
    dd=8;
    pE(end+5)=pE(end); pC(end+5,end+5)= pC(end,end);
    pE1=[.65 .41 .1 1 1]; pC1=[1 2 1 1 1]; %pE(9) à pE(13)
    for i1=1:length(pE1)
        pE(i1+dd)=pE1(i1);
        pC(i1+dd,i1+dd)=pC1(i1);
    end
    if ~isfield(M,'name')
        M.name={'SIGNAL_decay'  'FEEDBACK', 'TRANSIT_TIME','EXPONENT',...
            'EXTRACTION', 'VASCULAR_TONE' 'GAIN_PARAMETER', 'log_SIGNAL_RATIO' 'decay' 'd1' 'd2' 'd3' 'd4'   'effFlow'};
    end
    namePresent={'SIGNAL_decay'  'FEEDBACK', 'TRANSIT_TIME','EXPONENT',...
        'EXTRACTION', 'VASCULAR_TONE' 'GAIN_PARAMETER', 'log_SIGNAL_RATIO' 'decay' 'd1' 'd2' 'd3' 'd4'  'effFlow'};
    
end



%add scaling factor
namePresent=[ namePresent {'scFlow','scHbR'}];
pC=blkdiag(pC,eye(2,2));
pE=[pE ; 1; 1];
if isfield(M,'scaleGx') && isstruct(M.scaleGx)
    pE(end-1:end)=pE(end-1:end).*[M.scaleGx.Flow; M.scaleGx.HbR];
end





pE(pE==0)=min(pE(pE~=0));%obligatoire que cela vaille pas zéro

for i1=1:length(M.name)
    index(i1)= strmatch(M.name{i1},namePresent);
    pE2(i1,1)=pE(index(i1));
end
for i1=1:length(namePresent)
    M.var.(namePresent{i1})=pE(i1);
end


%normalise
pE=ones(length(M.name),1 );
pC=pC(index,index);
pC2= pC./ (pE2*pE2');
pC(isfinite(pC2))=pC2(isfinite(pC2));
pC=min(pC,1);

if isfield(M,'covarianceIni') && M.covarianceIni>0
    pC=eye(length(pE))*M.covarianceIni;
end
%remove nul variance
dd=diag(pC);
pC=pC-diag(dd);
dd(dd==0)=0.04;
pC=pC+diag(dd);

M.pE=pE;
M.pC=pC;
M.nameIni=namePresent;

1;