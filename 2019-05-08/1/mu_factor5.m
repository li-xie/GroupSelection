function y=mu_factor5(n)
s1=0.05;
s2=0.05;
u=rand(n,1);
y=zeros(size(u));
div_value = s2*(1-exp(-1/s2))/(s1+s2*(1-exp(-1/s2)));
idx=find(u<=div_value);
y(idx)=log((s1+s2*(1-exp(-1/s2))) * u(idx)/s2+exp(-1/s2)) * s2;
idx=find(u>div_value);
y(idx)=-s1*log((1-u(idx))*(s1+s2*(1-exp(-1/s2)))/s1);
% y(y<-1)=-1;










