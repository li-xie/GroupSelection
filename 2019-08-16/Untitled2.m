
xc = (1e-8:1e-8:1e-3)';
v = 3.471e-6*ones(size(xc));
y = GoodEqs([v, xc]);
[m,p]=min(y(:,1))
figure(1)
% plot(xc, log10(abs(y(:,1))))

axis([0 1e-3 -10,10])
figure(2)
plot(xc, y(:,2))