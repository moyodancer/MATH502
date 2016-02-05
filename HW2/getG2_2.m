function [ cj ] = getG2_2( col, m, h, x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if col == 0
    cj = x-1;
elseif col == m+1
    cj = ones(length(x), 1);
else
    for i = 1:length(x)
        if i<=col
            cj(i) = h*(x(col)-1);
        else
            cj(i) = h*(x(i)-1);
        end
    end
end

end

