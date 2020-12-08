function [x_n,e_n] = next_x(x,B,u,p_n,r)
%NEXT_X Summary of this function goes here
%   Detailed explanation goes here
    x_n=((1+r)*x)-sum((B.*u));
    e_n=x_n-p_n;
end

