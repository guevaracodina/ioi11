function eps_pathlength = ioi_epsilon_pathlength(lambda1,lambda2,npoints,whichSystem,whichCurve,debug)

%	This function estimates epsilon * D, it takes into account the camera
%	response, the leds spectra and uses a pathlength factor either set from
%	Kohl or Dunn in the literature.

%   This module is dependent on this file which contains all hardware info
%   for the setup, needs to specify the leds and the camera response.
%   we are still a bit dependent on the RGY but we could abstract this
%   (however the lambdas would need to be registered to specific hardware
%   still
if whichSystem
    load hardware_newsystem.mat;
else
    load hardware.mat;
end

% Rough baseline concentrations (in uM) : 100 uM (in the brain)
c_tot = 100e-6;

lambda_vec= linspace(lambda1,lambda2,npoints);
c_camera = private_reinterpolate_lambda(lambda1, lambda2, npoints, hardware.camera.lambda, hardware.camera.response);
c_led(1,:) = private_reinterpolate_lambda(lambda1, lambda2, npoints, hardware.leds.r.lambda, hardware.leds.r.intensity);
c_led(2,:) = private_reinterpolate_lambda(lambda1, lambda2, npoints, hardware.leds.g.lambda, hardware.leds.g.intensity);
c_led(3,:) = private_reinterpolate_lambda(lambda1, lambda2, npoints, hardware.leds.y.lambda, hardware.leds.y.intensity);
c_pathlength = ioi_path_length_factor(lambda1, lambda2, npoints, c_tot*1000, whichCurve);
[c_ext_hbo,c_ext_hbr] = ioi_get_extinctions(lambda1,lambda2,npoints);

if nargin>4
    figure;
    subplot(2,2,1)
    plot(lambda_vec,c_led(1,:),'r')
    hold on
    plot(lambda_vec,c_led(2,:),'g')
    hold on
    plot(lambda_vec,c_led(3,:),'y')
end

%Input basehbt1 not used...

% Create vectors of values for the fits
CHbO = 0.6*c_tot*(.5:.1:1.5);
CHbR = 0.4*c_tot*(.5:.1:1.5);

% In this computation below we neglect the fact that pathlength changes
% with total concentration (it is fixed for a Ctot of 100e-6)
eps_pathlength=zeros(3,2);

for iled=1:3
    for iconc = 1:length(CHbO)
        IHbO(iconc) = sum(c_camera .* c_led(iled,:) .* exp(-c_ext_hbo .* c_pathlength * CHbO(iconc)),2) ; %	Measured intensity for different concentrations
        IHbR(iconc) = sum(c_camera .* c_led(iled,:) .* exp(-c_ext_hbr .* c_pathlength * CHbR(iconc)),2) ;
    end
    IHbO = IHbO/max(IHbO);
    IHbR = IHbR/max(IHbR);

    % Compute effective eps
    p1 = polyfit(CHbO,-log(IHbO),1);
    p2 = polyfit(CHbR,-log(IHbR),1);
    HbRL = p2(1); %epsilon*D HbR effectif
    HbOL = p1(1);%epsilon*D HbO effectif
    eps_pathlength(iled,1)=HbOL;
    eps_pathlength(iled,2)=HbRL;
end


function interpolated=private_reinterpolate_lambda(lambda1, lambda2, npoints, i_lambda, i_intensity)

% Wanted values
xi = linspace(lambda1,lambda2,npoints);
% Actual values we have
x = i_lambda; 
y = i_intensity;

% watch for boundaries (extrapolation) set to zero in all cases
if x(1)>lambda1
x=[lambda1 ;x(1)*.9999 ;x];
y=[0;0; y];
end

if x(end)<lambda2
x=[x ; x(end)*1.0001 ;lambda2];
y=[y;0;0];
end

% perform interpolation
interpolated = interp1(x,y,xi); 

