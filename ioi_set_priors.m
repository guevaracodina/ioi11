function M=ioi_set_priors(M)
%load priors
switch M.PS.PhysioModel_Choice
    case 0 %Buxton-Friston
        M = ioi_SCKS_priors_Buxton(M,5);
%     case 1 %Zheng-Mayhew
%         SCKS = ioi_HDM_priors_ZM(SCKS,9); 
%     case 2 %Huppert1
%         SCKS = ioi_HDM_priors_Hu(SCKS,9);
%     otherwise
end