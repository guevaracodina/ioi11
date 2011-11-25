function M = ioi_set_physiomodel(M)
% Physiological Model choice
switch M.PS.PhysioModel_Choice
    case 0 %Buxton-Friston
        M.f     = 'ioi_fx';
        M.g     = 'ioi_gx';
        M.x     = zeros(4,1);       
    case 1 %Zheng-Mayhew
        M.f     = 'ioi_HDM_fx_ZM';  
        M.g     = 'ioi_HDM_gx';  
        M.x     = zeros(5,1);    
    case 2 %Huppert1
        M.f     = 'ioi_HDM_fx_Hu2';  
        M.g     = 'ioi_HDM_gx';  
        M.x     = zeros(7,1);           
end
M.n     = length(M.x);
%Number of inputs of direct model
M.m=1;

