function M=ioi_HDM_set_model(M)
% Model choice
M.m=1;
M.N=64; %default 64

if ~isfield(M,'IS')
    M.IS='spm_int';
end

if ~isfield(M,'dispDo')
    M.dispDo=[0 0 0 0 0]; %brut/ fitInitial/ doFit / dataFit/ diffPrior /spmresults
end

%%%M.TE=0.03; % optional param for BOLD
%--------------------------------------------------------------------------
switch M.Model_Choice
    case 0 %Buxton-Friston
        M.f     = 'ioi_HDM_fx';
        M.g     = 'ioi_HDM_gx';
        M.x     = zeros(4,1);       
    case 1 %Zheng-Mayhew
        M.f     = 'ioi_HDM_fx_ZM';  %MODIFIER
        M.g     = 'ioi_HDM_gx';  %MODIFIER
        M.x     = zeros(5,1);    %MODIFIER
    case 2 %Huppert1
%         M.f     = 'IOI_fx_hdm_Hu';  %MODIFIER
        M.f     = 'ioi_HDM_fx_Hu2';  %MODIFIER
        M.g     = 'ioi_HDM_gx';  %MODIFIER
        M.x     = zeros(7,1);      %MODIFIER
    case 3 %zheng decay couplage
        M.f     = 'ioi_HDM_fx_ZM_Decay';  %MODIFIER
        M.g     = 'ioi_HDM_gx';  %MODIFIER
        M.x     =zeros(7,1);      %MODIFIER
%     case 4 %Buxton decay couplage
%         M.f     = 'ioi_HDM_fx_Decay';  %MODIFIER
%         M.g     = 'ioi_HDM_gx';  %MODIFIER
%         M.x     =zeros(6,1);      %MODIFIER
%     case 5 %Buxton decay couplage
%         M.f     = 'ioi_HDM_fx_test';  %MODIFIER
%         M.g     = 'ioi_HDM_gx';  %MODIFIER
%         try
%             [a n]=feval(M.f,[],[],[],M)
%         end
%         M.x     =zeros(n,1);      %MODIFIER
%     case 6 %Buxton decay couplage
%         M.f     = 'ioi_HDM_fx_BuDcy';  %MODIFIER
%         M.g     = 'ioi_HDM_gx';  %MODIFIER
%         M.x     =zeros(5,1);      %MODIFIER
%     case 7 %Hu decay couplage
%         M.f     = 'ioi_HDM_fx_HuDcy';  %MODIFIER
%         M.g     = 'ioi_HDM_gx';  %MODIFIER
%         M.x     =zeros(8,1);      %MODIFIER
% 
%     case 8 %Hu decay couplage
%         M.f     = 'ioi_HDM_fx_Zheng2010';  %MODIFIER
%         M.g     = 'ioi_HDM_gx';  %MODIFIER
%         M.x     =zeros(7,1);      %MODIFIER
%  case 9 %Hu decay couplage
%         M.f     = 'ioi_HDM_fx_noFlow';  %MODIFIER
%         M.g     = 'ioi_HDM_gx';  %MODIFIER
%         M.x     =zeros(2,1);      %MODIFIER
end

 M.n     = length(M.x);

if isfield(M,'curveToFitName')
    name={'HbO','HbR','HbT','Flow','CMRO'};
    M.YName={};
    for i1=1:length(name)
        if strfind(M.courbeToFitName,name{i1})
            M.YName{end+1}=name{i1};
        end
    end
else
    switch M.curveToFit
        case 1
            M.YName={'Flow'};
        case 2
            M.YName={'HbR','HbT'};
        case 3
            M.YName={'HbR','HbT','Flow'};
        case 4
            M.YName={'HbT','Flow'};
        case 5
            M.YName={'HbO','HbT'};
        case 6
            M.YName={'HbO','HbT','Flow'};
        case 7
            M.YName={'CMRO'};
        case 8
            M.YName={'HbR','HbT','HbO','Flow'};
    end
end
