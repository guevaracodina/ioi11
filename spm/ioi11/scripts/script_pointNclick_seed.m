%% image setup
load mri;
currentimage = D(:,:,1,15);
figure;
map = colormap('gray');
imagesc(currentimage);
axis square;

% roisetup
t = 0:pi/100:2*pi;
% R0 = 5; x0 = 50; y0 = 50;
R0 = 3;
% Prompt user to point & cl
p = ginput(1);
x0 = p(1);
y0 = p(2);
% Parametric function for a circle
xi = R0*cos(t)+x0;
yi = R0*sin(t)+y0;
LineHandler = line(xi,yi,'LineWidth',3,'Color',[.8 0 0]);

% calc. roi stat.
roimask = poly2mask(xi,yi, size(currentimage,1), size(currentimage,2));
pr_r = find(roimask);
roimean = mean(currentimage(pr_r));
roistd = std(double(currentimage(pr_r)));
[roimean roistd]

% Display ROI
figure; imagesc(roimask)
axis image
colormap(map)
