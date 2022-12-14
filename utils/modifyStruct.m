function [data] = modifyStruct(data,my_function, args)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if isstruct(data)
    names = fieldnames(data);
    for i=1:numel(names)
        data.(names{i}) = modifyStruct(data.(names{i}), my_function, args);
    end
end

if iscell(data)
    for i=1:numel(data)
        data{i} = modifyStruct(data{i}, my_function, args);
    end
end

if isnumeric(data)
    data = my_function(data, args);
end
end


