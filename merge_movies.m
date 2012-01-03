function merge_movies
%This is a utility to combine movies created for each onset type into
%the same onset type by simple averaging, since the onsets were really of
%the same type
colors = ['D' 'F' 'O' 'T'];
sessions = [1 2 4];
stims = 3;

for c1=1:length(colors)
    for s1=1:length(sessions)
        for m1=1:stims
            fname = ['cine_S' gen_num_str(sessions(s1),2) '_' colors(c1) '_onset' int2str(m1) '.Mdat'];
            load(fname,'-mat');
            if m1 == 1
                d0 = zeros(size(d));
            end
            d0 = d+d0;
        end
        d = d0/stims;
        fname = ['cine_S' gen_num_str(sessions(s1),2) '_' colors(c1) '.Mdat'];
        save(fname,'d');
        disp(['Session ' int2str(sessions(s1)) ', Color ' colors(c1) ' complete']);
    end
end