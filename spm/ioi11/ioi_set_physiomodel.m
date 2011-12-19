function M = ioi_set_physiomodel(M)
% Physiological Model choice
M.f     = 'ioi_fx';
M.g     = 'ioi_gx';
switch M.O.PhysioModel_Choice
    case 0 %Buxton-Friston
        M.x     = zeros(4,1);
    case 1 %Zheng-Mayhew
        M.x     = zeros(5,1);    
    case 2 %Huppert1
        M.x     = zeros(7,1);           
end
M.n     = length(M.x);
%Number of inputs of direct model
M.m=1;
%Number of outputs
M.l = M.O.includeHbR + M.O.includeHbT + M.O.includeFlow;


