% N is a column  vector of interger numbers. 
% p is the probability of drawing
% different from Alex's original code: changed the thresholds for
% approximations
function [result] = fastbinorv(N, p)
% CHECK TO MAKE SURE EVERYTHING IS RIGHT!!!
pcs = 1e-12;
if isempty(N)
    error('you are calling fastbniorv on an empty vector. This is a waste of everyone''s time.');
end

[~,n] = size(N);
if n > 1
    error('N must be a column vector')
end

% the minimal n for approximating the binomial distribution with a Poisson distribution 
thresh_b2p = 1e2;
% the minimal n for approximating the binomial distribution with a normal distribution 
thresh_p2n = 1e3;
% when p is not too close to 0 or 1, normal approximation is appropriate
normcond = @(N,p) N > 10 ./ p & N > 9 * p ./ (1-p);
% generate an outcome column vector
result = zeros(length(N),1);
if length(p) == 1
    p = ones(length(N),1) * p;
    result = N;
end
% when the normal approximation condition is satisfied for all elements in N,  use a normal approximation
if and(min(N) > thresh_p2n , normcond(N,p))
    result = round(normrnd(N .* p, sqrt(N .* p .* (1 - p))));
% decide whether one should use normal/Poisson approximation
else
    % class_ind = 1, normal approximation
    % class_ind = 2, Poisson approximation
    % class_ind = 3, Binomial random number generator
    class_ind = zeros(size(N), 'uint8');
    norm_ind = find(N >= thresh_p2n & normcond(N, p));
    class_ind(norm_ind) = 1;
%     pois_ind = find(N >= thresh_b2p & N <= 10 ./ p);
    pois_ind = find(N >= thresh_b2p & (N < thresh_p2n | ~normcond(N,p)));
    class_ind(pois_ind) = 2;
    bino_ind = find(N > pcs & N < thresh_b2p);
    class_ind(bino_ind) = 3;
    if nnz((class_ind == 0) & (N > pcs))>0
        error('missed a category')
    end
    if ~isempty(bino_ind)
        randy = randi(10^8);
        result(bino_ind) = (mybinornd_mex(int32(randy),(N(bino_ind))', (p(bino_ind))'))';
    end
    result(pois_ind) = poissrnd(N(pois_ind) .* p(pois_ind));
    result(norm_ind) = round(normrnd(N(norm_ind) .* p(norm_ind), sqrt(N(norm_ind) .*...
        p(norm_ind) .* (1 - p(norm_ind)))));
end
if sum(sum(result<0 + result>N)) > 0
    error('the Poisson/Normal approximation might not be good')
end
result(result < 0) = result(result < 0) * 0;
result(result > N) = result(result > N) * 0 + N(result > N);
end