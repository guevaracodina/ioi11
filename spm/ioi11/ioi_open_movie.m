function [obj F Y] = ioi_open_movie(file,pathname)
file = fullfile(pathname,file);
[dir0 fil0 ext0] = fileparts(file);
switch ext0
    case '.avi'
        obj = VideoReader(file);
        nF = obj.NumberOfFrames;
        % Preallocate movie structure.
        F(1:nF) = struct('cdata',...
            zeros(obj.Height,obj.Width,3,'uint8'),'colormap',[]);
        % Read one frame at a time.
        for k=1:nF
            F(k).cdata = read(obj,k);
        end
    case '.mat'
        open(file);
end
[X dummy] = frame2im(F(1));
[nx ny] = size(X);
Y = zeros(nx,ny,nF);
for k=1:nF
    y = frame2im(F(k));
    Y(:,:,k) = y;
end