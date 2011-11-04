function M = ioi_HDM_priors_Buxton(M,h)
% returns priors for a hemodynamic dynamic causal model
%select Specific prior to do

[pE,pC] = ioi_HDM_priors_base_case(M.m,5);
pE(end)=0.025; % put a better approximation around effFlow to extract good estimate of covariance

if ~isfield(M,'name')
    M.name={'SIGNAL_decay', 'FEEDBACK', 'TRANSIT_TIME',...
        'EXPONENT', 'EXTRACTION',  'log_SIGNAL_RATIO' 'effFlow'};
    % M.name={'SIGNAL_decay',  'TRANSIT_TIME',...
    %     'EXPONENT', 'EXTRACTION' 'effFlow'};
end

namePresent={'SIGNAL_decay', 'FEEDBACK', 'TRANSIT_TIME',...
    'EXPONENT', 'EXTRACTION',  'log_SIGNAL_RATIO' 'effFlow'};

 
 if any(M.Model_Choice==[4 5 6 8])%decay
    dd=6;
    pE(end+5)=pE(end); pC(end+5,end+5)= pC(end,end);
    pE1=[.65 .281 .1 .01 1]; pC1=[1 2 1 1 1]; %pE(9) à pE(13)
    for i1=1:length(pE1)
        pE(i1+dd)=pE1(i1);
        pC(i1+dd,i1+dd)=pC1(i1);
    end
    if ~isfield(M,'name')
        M.name={'SIGNAL_decay', 'FEEDBACK', 'TRANSIT_TIME',...
    'EXPONENT', 'EXTRACTION',  'log_SIGNAL_RATIO'  'd0' 'd1' 'd2' 'd3' 'd4'   'effFlow'};
    end
    namePresent={'SIGNAL_decay', 'FEEDBACK', 'TRANSIT_TIME',...
    'EXPONENT', 'EXTRACTION',  'log_SIGNAL_RATIO'  'd0' 'd1' 'd2' 'd3' 'd4'  'effFlow'};
    
end
 
 %add scaleGx
namePresent=[ namePresent {'scFlow','scHbR'}];
 pC=blkdiag(pC,zeros(2,2));
 pE=[pE ; 1; 1];
 if isfield(M,'scaleGx') && isstruct(M.scaleGx)
    pE(end-1:end)=pE(end-1:end).*[M.scaleGx.Flow; M.scaleGx.HbR];
 end
 
 
for i1=1:length(M.name)
    try
        index(i1)= strmatch(M.name{i1},namePresent);
        pE2(i1,1)=pE(index(i1));
    catch
%          warning(['variable pas présente : ' M.name{i1}]) %let warning if
 %         variables are not all present in convolution performed.
% will take default value   
        error(['variable pas présente : ' M.name{i1}])
    end   
end






pE(pE==0)=min(pE(pE~=0));%obligatoire que cela vaille aps zéro
for i1=1:length(namePresent)
    M.var.(namePresent{i1})=pE(i1);
end

 

%
%normalise
 pE=ones(length(M.name),1 );
 pC=pC(index,index);
 pC2= pC./ (pE2*pE2');
 pC(isfinite(pC2))=pC2(isfinite(pC2));
 pC=min(pC,100);
 
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

