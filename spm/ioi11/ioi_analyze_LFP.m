function A = ioi_analyze_LFP(E,el,ons,IOI,s,fdir)
n = length(ons);
sf = E.sf; %Hz, sampling frequency
tb = E.tb; %0.2; %time before, in seconds
ta = E.ta; %0.5; %time after, in seconds
Pf = 10; %rescaling factor for power, to put it on a similar scale in plots
lb = ceil(sf*tb);
la = ceil(sf*ta);
m = zeros(n,lb+la);
lp = linspace(-tb,ta,(ta+tb)*sf);
if E.use_epilepsy_convention
    el = -el;
end
for i=1:n
    st = ceil(ons(i)*sf);
    try
        %store electrophysiology data
        m(i,:) = el(st-lb+1:st+la);
    catch
        %data after last onset or before first onset might be incomplete
    end
end
h = figure; plot(lp,m','k')
xlabel('Time (s)')
ylabel('LFP (mV)')
title('Detected spikes, superposed')
filen1 = fullfile(fdir,['Onset_profile_S' gen_num_str(s,2) '.fig']);
saveas(h,filen1,'fig'); %save as .fig
filen2 = fullfile(fdir,['Onset_profile_S' gen_num_str(s,2) '.tiff']);
print(h, '-dtiffn', filen2);
close(h);

h = figure; plot(lp,m') %with colors
xlabel('Time(s)')
ylabel('LFP (mV)')
title('Detected spikes, superposed')
filen1 = fullfile(fdir,['Onset_profile_colors_S' gen_num_str(s,2) '.fig']);
saveas(h,filen1,'fig'); %save as .fig
filen2 = fullfile(fdir,['Onset_profile_colors_S' gen_num_str(s,2) '.tiff']);
print(h, '-dtiffn', filen2);
close(h);

%Calculate some properties of the spike
%max amplitude
aM = max(m(:,lb:lb+la),[],2);
%min amplitude
am = min(m(:,lb:lb+la),[],2);
%difference between max and min amplitude
ad = aM-am;
%area
A0 = mean(abs(m(:,lb:lb+la)),2);
%signed area
sA0 = mean(m(:,lb:lb+la),2);
%power
P = Pf*mean(m(:,lb:lb+la).*m(:,lb:lb+la),2);
%baseline
b = mean(m(:,1:lb),2);
%previous quantities normalized to baseline:
mb = m-repmat(b,1,lb+la);
%max amplitude wrt baseline
aMb = max(mb(:,lb:lb+la),[],2);
%min amplitude wrt baseline
amb = min(mb(:,lb:lb+la),[],2);
%area wrt baseline
Ab = mean(abs(mb(:,lb:lb+la)),2);
%signed area wrt baseline
sAb = mean(mb(:,lb:lb+la),2);
%power wrt baseline
Pb = Pf*mean(mb(:,lb:lb+la).*mb(:,lb:lb+la),2);

%Store
A{s}.n = n;
A{s}.m = m;
A{s}.b = b;
A{s}.aM = aM;
A{s}.am = am;
A{s}.ad = ad;
A{s}.A = A0;
A{s}.P = P;
A{s}.mb = mb;
A{s}.aMb = aMb;
A{s}.amb = amb;
A{s}.Ab = Ab;
A{s}.sAb = sAb;
A{s}.Pb = Pb;
A{s}.Y = [aM am ad A0 sA0 P b aMb amb Ab sAb Pb];
%Protocole
try
    A{s}.Pt = IOI.sess_res{s}.scan_info.protocoleStimulation;
end
A{s}.ons = ons;
%clear aM am ad A sA P b mb aMb amb Ab sAb Pb