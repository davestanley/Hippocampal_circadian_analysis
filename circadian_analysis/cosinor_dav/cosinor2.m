
function [MESOR,AMPL,PH,P]=cosinor2(y,t,T)
%[MESOR,AMPL,PH]=cosinor(y,t,T)
%
% Applies a Cosine fitting technique to the series y, taken
% at times t. T indicates the expected period of the Cosine.

% Matthias Fruehwirth 1.3.1996
% overall test added, mf 21.9.1998
% see -> F. Halberg et. al., Autorhythmometry - procedures for
% physiologic self-measurements and their analysis,
% The Physiology Teacher, Vol. I, No. 4, Jan 1972

w=2*pi/T;
x1=cos(w*t);
x2=sin(w*t);
sx1=sum(x1);
sx2=sum(x2);
sx1x2=sum(x1.*x2);

b=[ sum(y); ...
sum(x1.*y); ...
sum(x2.*y)];

a=zeros(3,3);
a(1,:)=[length(y) sx1 sx2];
a(2,:)=[sx1 sum(x1.^2) sx1x2];
a(3,:)=[sx2 sx1x2 sum(x2.^2)];

sol=a\b;
MESOR=sol(1);
AMPL=sqrt(sol(2)^2+sol(3)^2);
PH=atan2(-sol(3),sol(2));

yest=MESOR+AMPL*cos(w*t+PH);

% plot(t,y,'ob',t,yest,'-r',[min(t) max(t)],[MESOR MESOR],'--k')

f=sum((yest-mean(y)).^2)/sum((y-yest).^2)*(length(y)-3)/2;
P=1-fcdf(f,2,length(y)-3);

