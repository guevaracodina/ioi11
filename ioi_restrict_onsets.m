function [IOI onsets_list pars_list] = ioi_restrict_onsets(IOI,job,rmi,ust)
%select a subset of sessions
[all_sessions selected_sessions] = ioi_get_sessions(job);

%loop over sessions
for s1=1:length(IOI.sess_res)
    if all_sessions || sum(s1==selected_sessions)
        tmp_onsets = IOI.sess_res{s1}.onsets;
        tmp_pars = IOI.sess_res{s1}.parameters;
        %remove onsets
        if ~isempty(rmi)
            IOI.rmi = rmi;
            for k0=1:length(tmp_onsets)
                r_onsets = round(tmp_onsets{k0});
                k_onsets = [];
                k_pars = [];
                rem_onsets = [];
                for u0=1:length(r_onsets)
                    if any(r_onsets(u0) == rmi) || any(r_onsets(u0) == rmi + 1) ||  any(r_onsets(u0) == rmi - 1)
                        %skip this onset
                        rem_onsets = [rem_onsets tmp_onsets{k0}(u0)];
                    else
                        k_onsets = [k_onsets tmp_onsets{k0}(u0)];
                        k_pars = [k_pars tmp_pars{k0}(u0)];
                    end
                end
                IOI.onsets_kept{s1}{k0} = k_onsets;
                IOI.onsets_removed{s1}{k0} = rem_onsets;
                IOI.pars_kept{s1}{k0} = k_pars;
                tmp_onsets{k0} = k_onsets;
                tmp_pars{k0} = k_pars;
            end
        else
            IOI.onsets_kept{s1} = tmp_onsets;
            IOI.pars_kept{s1} = tmp_pars;
            IOI.onsets_removed{s1} = '';
        end
        if ~isempty(ust)
            for k0=1:length(tmp_onsets)
                try
                    tmp_onsets{k0} = tmp_onsets{k0}(ust);
                    tmp_pars{k0} = tmp_pars{k0}(ust);
                catch
                    disp('Problem restricting number of onsets to user-specified -- probably not enough onsets');
                end
                try
                    if ~any(k0 == job.which_onset_type)
                        tmp_onsets{k0} = ''; %remove these onset types
                        tmp_pars{k0} = '';
                    end
                end
            end
        end
        onsets_list{s1} = tmp_onsets;
        pars_list{s1} = tmp_pars;
    end
end