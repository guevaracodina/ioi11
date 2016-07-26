function fcli = ioi_lateralization_idx (LL, LR, RL,  RR)
% ioi_laterality_idx Computes the functional connectivity lateralization index
% From:
% SYNTAX
% fcli = ioi_lateralization_idx (LL, LR, RL,  RR)
% INPUTS
% LL, LR, RL,  RR
% OUTPUT
% fcli
fcli = ((LL - RL) - (RR - RL)) ./ (abs(LL) + abs(LR) + abs(RL) + abs(RR));

end

