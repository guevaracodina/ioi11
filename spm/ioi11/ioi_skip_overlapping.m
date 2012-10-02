function  [onsets_list skipped] =ioi_skip_overlapping(skip_overlap,onsets_list,...
    which_onset_type,window_after,window_before)
%Function to skip_overlapping onsets
skipped = 0;
if skip_overlap
    %loop over onset types
    for m1=1:length(onsets_list)
        if any(m1==which_onset_type)
            if length(onsets_list{m1}) > 1
                tmp = [];
                for o1=1:(length(onsets_list{m1})-1)
                    if onsets_list{s1}{m1}(o1+1)-onsets_list{m1}(o1) > window_after+window_before %in seconds
                        tmp = [tmp onsets_list{m1}(o1)];
                    else
                        skipped = skipped + 1;
                    end
                end
                %always keep the last one
                tmp = [tmp onsets_list{m1}(end)];
                onsets_list{m1} = tmp;
            end
        end
    end
end