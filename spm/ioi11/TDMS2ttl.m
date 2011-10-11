function out=TDMS2ttl(in)
% convert the TDMS to ttl for led camera and 


% same thing as num2bin, but go around user license problem
ttl1=rem(in.Data.MeasuredData(7).Data/1,2)>=1; %camera trigger 
ttl2=rem(in.Data.MeasuredData(7).Data/2,2)>=1; %stimulator trigger
ttl3=rem(in.Data.MeasuredData(7).Data/4,2)>=1; %color trigger
out=[ttl1 ttl2 ttl3]==1;
%plot([ttl1 ttl2 ttl3])
