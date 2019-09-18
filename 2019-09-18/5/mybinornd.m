function [result] = mybinornd(N, p)
%MYBINORND Summary of this function goes here
%   Detailed explanation goes here
len = length(N);
result = zeros(len,1);
for i = 1 : len
    x = int32(0);
    for j = 1 : N(i)
        x = x + int32(rand < p(i));
    end
    result(i) = x;
end
end