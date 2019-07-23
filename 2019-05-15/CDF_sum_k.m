function y = CDF_max_k(x, a, b, k)
y = zeros(size(x));
for i = 1:k
    y = y + exp(-exp(-(x-b)/a)-(x-b)*(i-1)/a)/factorial(i-1)*(k-i+1)/k;
end
