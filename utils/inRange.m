function [idx] = inRange(idx, args)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
coordinates = args{1};
limits = args{2};

idx = idx(all(coordinates(idx,:) >=limits(1,:),2) & all(coordinates(idx, :) <= limits(2,:),2));

end

