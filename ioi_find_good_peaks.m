function [npkh npk dur] = ioi_find_good_peaks(pk,pkh,dP,el,E,MN,eSD)
npk = [];
npkh = [];
tmpdur = [];
dur = [];
%select fast and slow rising onsets
select_fast = 1; %peak reached within 40 ms of onsets
fast_pk = 400; %data points at 10000 Hz
select_slow = 1;
cf = 0.5; %?
%Narrow down to good events, by making sure it goes above then below baseline
%if E.spkOn
for i=1:length(pk)
    good = 0;
    %check that next peak is at least dP points later
    if (i>1 && pk(i)-pk(i-1) > dP) || (i==1)
        if i== 1
            pkm = -Inf;
        else
            %minimum over interval of 2*dP points
            pkm = min(el(pk(i-1):pk(i)));
        end
        if pkm < MN-cf*eSD
            %went below baseline, found new peak
            
            %calculate max
            [dummy pkM]= max(el(pk(i):min(length(el),(pk(i)+3*fast_pk))));
            
            if select_fast
                if pkM <= fast_pk
                    %if pkM > 1
                    good =1;
                    %end
                end
            end
            if select_slow
                if  pkM> fast_pk
                    good = 1;
                end
            end
            %exclude if previous spike was very close and current spike
            %is not deep -- to do
            
            if good
                npk =[npk pk(i)];
                npkh = [npkh pkh(i)];
                lpeak = pk(i);
                dur = [dur tmpdur]; %this cannot give the duration of the last seizure...            
            end
        end
    else
        tmpdur = (pk(i)-lpeak)/E.sf; %in seconds
    end
end
%last peak
try
    tmpdur = (pk(end)-lpeak)/E.sf; %in seconds
    dur = [dur tmpdur];
end