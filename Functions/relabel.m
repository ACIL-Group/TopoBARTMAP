%This file is for converting string array labels to numeric labels.

function output = relabel(input) 
    unique_elements = unique(input);
    shape_unique = size(unique_elements,2);
    output = zeros(size(input));
    for iter = 1:1:shape_unique
        output = output + iter*(input==unique_elements(iter));
    end
end