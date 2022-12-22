%{
This function is a part of implementation of Clustering Error (CE) an external
bicluster evaluation metric introduced by Patrikainen and Meila (2006).

This implimentation is based on Biclustlib, by Victor Alexandre Padilha

%}

function value = calculate_size(predicted, reference, num_rows, num_cols,operation)
    
    pred_count = zeros(num_rows,num_cols);
    ref_count = zeros(num_rows,num_cols);

    for itr=1:1:size(predicted,1)
        pred_count(predicted{itr,1},predicted{itr,2})=1;
    end

    for itr=1:1:size(reference,1)
        ref_count(reference{itr,1},reference{itr,2})=1;
    end
    if operation=="union" || operation=="Union"
        value = sum(max(pred_count,ref_count),'all');
    elseif operation=="intersection" || operation=="Intersection"
        value = sum(min(pred_count,ref_count),'all');
    else
        error("Incorrect operation specified, valid operations are: union and intersection");
    end

end