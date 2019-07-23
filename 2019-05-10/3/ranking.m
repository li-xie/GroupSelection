function [result] = ranking(vec,dir)
if ~ or(isequal(dir, 'ascend'), isequal(dir, 'descend'))
    error('dir must be either ascend or descend')
end
[~,p] = sort(vec,dir);
result = 1:length(vec);
result(p) = result;
result = result';
end