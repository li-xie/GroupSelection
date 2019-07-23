function y = mu_specZeroNull(n)
u=rand(n, 1);
y=zeros(size(u));


% y(u <= 0.1) = -1;

idx = find(u > 0.9);
if ~isempty(idx)
    y(idx)=mu_factorDunham(length(idx));
end

