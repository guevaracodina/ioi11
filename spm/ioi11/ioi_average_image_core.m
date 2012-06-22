function [IOI,D,Z] = ioi_average_image_core(IOI,y,Z)
%y: filtered images
do_another_stim = 0; %Boolean, only turned on in interactive mode
onset_first_pass = 1;
[nt nxny] = size(y);
kb = 0; %counter of segments before onsets
ka = 0; %counter of segments after onsets
kb2 = 0; %skipped
ka2 = 0; %counter of skipped segments after onsets
%Ma = zeros(1,nxny); %mean
Sa = [];
%Mb = zeros(1,nxny);

try
    while onset_first_pass || Z.interactive_mode || do_another_stim
        onset_first_pass = 0;
        if ~onset_first_pass && Z.interactive_mode && do_another_stim
            %display interactive figure
        end
        
        %tmp_array_before = zeros(Z.window_before,nxny);
        %tmp_array_after = zeros(Z.window_after,nxny);
        %loop over onsets for that session
        if ~isempty(Z.ons)
            U = round(Z.ons/IOI.dev.TR)-Z.window_offset; %in data points
            %U0{s1} = ioi_get_U(IOI,[],U,0,s1); %only for plotting stims
            for u1=1:length(U)
                clear tmp_median;
                try
                    %Baseline
                    tmp1 = y((U(u1)-Z.window_before):(U(u1)-1),:);
                    switch Z.normalize_choice
                        case 1
                            tmp_median = median(tmp1,1);
                        case 2
                            tmp_median = tmp1(end,:);
                        case 3
                            tmp_median = mean(tmp1,1);
                    end
                    kb = kb+1;
                    %Mb = Mb + mean(tmp1,1) - tmp_median;
                catch
                        kb2 = kb2+1;
                end
                try
                    %window after
                    tmp1 = y((U(u1)+Z.window_start_delay):(U(u1)+Z.window_start_delay+Z.window_after-1),:);
                    tmp0 = mean(tmp1,1);
                    if ~exist('tmp_median','var') || Z.normalize_choice == 2
                        tmp_median = tmp1(1);
                    end
                    if Z.remove_segment_drift
%                         switch Z.normalize_choice
%                             case {1,3}
%                                 %use the same length as window before to estimate end value of segment
%                                 tmp_end = mean(tmp1(end-Z.window_before:end,:),1);
%                             case 2
%                                 tmp_end = tmp1(end,:);
%                         end
%                         slope = tmp_median + linspace(0,1,size(tmp1,1))*(tmp_end-tmp_median);
%                         tmp1 = tmp1 - slope;
                    else
                        tmp0 = tmp0-tmp_median;
                    end
                    %tmp_array_after = tmp_array_after + tmp1;
                    ka = ka+1;
                    %Ma = Ma + tmp0;
                    Sa = [Sa;tmp0];
                catch
                        ka2=ka2+1;
                end
            end
            
            %Mb{r2,m1}{c1,s1} = tmp_array_before/kb; %global mean before
            Ma = mean(Sa,1); %Ma/ka; %global mean after
            Da = std(Sa,0,1);
        end
    end
    
    if (ka2>0 || kb2 > 0)
        disp(['Skipped ' int2str(ka2) ' segments after onsets and ' int2str(kb2) ' segments before onsets']);
    end
    disp(['ka: ' int2str(ka) ', kb: ' int2str(kb)])
    disp(['ka2: ' int2str(ka2) ', kb2: ' int2str(kb2)])
    D.Ma = Ma;
    D.Da = Da;
    D.ka = ka;
    D.kb = kb;
    D.ka2 = ka2;
    D.kb2= kb2;
    
catch exception
    disp(exception.identifier)
    disp(exception.stack(1))
end