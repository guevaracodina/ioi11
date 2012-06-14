function [ons amp IOI] = ioi_remove_onsets(ons, amp, rmi, ust, IOI,s1,k0)
%remove onsets
if ~isempty(rmi)
    IOI.rmi = rmi;
    r_onsets = round(ons);
    k_onsets = [];
    k_amp = [];
    rem_onsets = [];
    rem_amp = [];
    for u0=1:length(r_onsets)
        if any(r_onsets(u0) == rmi) || any(r_onsets(u0) == rmi + 1) ||  any(r_onsets(u0) == rmi - 1)
            %skip this onset
            rem_onsets = [rem_onsets ons(u0)];
            try rem_amp = [rem_amp amp(u0)]; end
        else
            k_onsets = [k_onsets ons(u0)];
            try k_amp = [k_amp amp(u0)]; end
        end
    end
    IOI.onsets_kept{s1}{k0} = k_onsets;
    IOI.onsets_removed{s1}{k0} = rem_onsets;
    ons = k_onsets;
    try
        IOI.amp_kept{s1}{k0} = k_amp;
        IOI.amp_removed{s1}{k0} = rem_amp;
        amp = k_amp;
    end
else
    IOI.onsets_kept{s1}{k0} = ons;
    IOI.onsets_removed{s1}{k0} = '';
    try
        IOI.amp_kept{s1}{k0} = amp;
        IOI.amp_removed{s1}{k0} = '';
    end
end
if ~isempty(ust)
    try
        ons = ons(ust);
    catch
        disp('Problem restricting number of onsets to user-specified -- probably not enough onsets');
    end
    try
        if ~any(k0 == job.which_onset_type)
            ons{k0} = ''; %remove these onset types
        end
    end
end