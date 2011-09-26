%Explore onsets (LFPs): shape of unfiltered response
%load('IOI.mat')
Ns = 7;
Pf = 10; %rescaling factor for power, to put it on a similar scale in plots
for s = 1:Ns   
    load(['el_S0' int2str(s) '.mat']);
    ons = IOI.Sess(s).U.ons;
    n = length(ons);
    sf = 10000; %Hz, sampling frequency
    tb = 1; %time before, in seconds
    ta = 1; %time after, in seconds
    lb = ceil(sf*tb);
    la = ceil(sf*ta);
    m = zeros(n,lb+la);
    for i=1:n
        st = ceil(ons(i)*sf);
        try
            %store electrophysiology data
            m(i,:) = el(st-lb+1:st+la);
        catch
            %data after last onset or before first onset might be incomplete
        end
    end
    figure; plot(m')
    %Calculate some properties of the spike
    %max amplitude
    aM = max(m(:,lb:lb+la),[],2);
    %min amplitude
    am = min(m(:,lb:lb+la),[],2);
    %difference between max and min amplitude
    ad = aM-am;
    %area
    A = mean(abs(m(:,lb:lb+la)),2);
    %signed area
    sA = mean(m(:,lb:lb+la),2);
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
    E{s}.n = n;
    E{s}.m = m;
    E{s}.b = b;
    E{s}.aM = aM;
    E{s}.am = am;
    E{s}.ad = ad;
    E{s}.A = A;
    E{s}.P = P;
    E{s}.mb = mb;
    E{s}.aMb = aMb;
    E{s}.amb = amb;
    E{s}.Ab = Ab;
    E{s}.sAb = sAb;
    E{s}.Pb = Pb;
    E{s}.Y = [aM am ad A sA P b aMb amb Ab sAb Pb];
    %Protocole
    E{s}.Pt = IOI.sess_res{s}.scan_info.protocoleStimulation;
    E{s}.ons = ons;
    %clear aM am ad A sA P b mb aMb amb Ab sAb Pb
end
save('E','E');
Zm = []; Zs = []; Zsk = []; Zk = []; Pt = []; Zn = [];
for s = 1:Ns
    Y = E{s}.Y;
    Ym = mean(Y,1);
    Ys = std(Y,0,1);
    Ysk = skewness(Y,0,1);
    Yk = kurtosis(Y,0,1);
    Zm = [Zm Ym']; 
    Zs = [Zs Ys'];
    Zsk = [Zsk Ysk'];
    Zk = [Zk Yk'];
    Pt = [Pt; E{s}.Pt];
    Zn = [Zn; E{s}.n];
end
figure; plot(Zm(1:7,:)');
figure; plot(Zm([3 11 12],:)');
figure; plot(Zs([3 11 12],:)');
figure; plot(Zsk([3 11 12],:)');
figure; plot(Zk([3 11 12],:)');
figure; plot(m')
p = sum(m.*m,2);
Y = [aM am ad A sA P b aMb amb Ab sAb Pb];
Y = [aM am ad A sA b aMb amb Ab sAb];
figure; plot(Y);
Ym = mean(Y,1);
Ys = std(Y,0,1);
Ysk = skewness(Y,0,1);
Yk = kurtosis(Y,0,1);
figure; plot([Ym' Ys' Ysk' Yk']);

