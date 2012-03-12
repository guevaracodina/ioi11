function [score scoGld scoDtc ]= ioi_sensi(varargin)
%Gld: gold standard
%Dtc: data series to search 
%dst (e.g. dt): window size for detection, in seconds
%S: ROC structure
Gld = varargin{1};
Dtc = varargin{2};
S = varargin{3};
if nargin == 4
    mask0 = varargin{4};
    indexmask = ceil(mask0/S.freq/S.dtwindow);    
end

%find sensitivity matrix
nBox=ceil(S.n/(S.freq*S.dtwindow));

aa=zeros(nBox,1);
bb=aa;
indexdtc=ceil(Dtc/S.freq/S.dtwindow);

if nargin == 4
    mask = aa;
    mask(indexmask) = 1;
    score.mask0 = mask0;
    score.vct.mask = mask;
end
if nargin == 3
    indexgld=ceil(Gld/S.freq/S.dtwindow);
    aa(indexgld)=1;
else
    if nargin == 4
        aa = Gld';
    end
end
bb(indexdtc)=1;
score.gldp2=length(Gld);
score.dtcp2=length(Dtc);
if nargin == 3
    score.tp2=sum(aa&bb);
    score.fn2=sum(aa&~bb);
    score.fp2=sum(~aa&bb);
    score.tn2=sum(~aa&~bb);
else
    score.tp2=sum(aa&bb&~mask); %exclusive mask
    score.fn2=sum(aa&~bb&~mask);
    score.fp2=sum(~aa&bb&~mask);
    score.tn2=sum(~aa&~bb&~mask);
end
score.sensi2=score.tp2/(score.tp2+score.fn2);
score.speci2=score.tn2/(score.tn2+score.fp2);
score.gldpct2=score.tp2/length(Gld);
score.dtcPct2=score.tp2/length(Dtc);
score.vct.gld2=aa';
score.vct.dtc2=bb';
score.vct.temps2=(1:nBox)*S.dtwindow';
