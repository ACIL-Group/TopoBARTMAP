%{
This function is an implementation of Clustering Error (CE) an external
bicluster evaluation metric introduced by Patrikainen and Meila (2006).

This implimentation is based on Biclustlib, by Victor Alexandre Padilha

Param: predicted_biclustering, reference_biclustering: cell array with
{[Bi rows],[Bi columns]; [Bj rows],[Bj columns]} form
%}


function value = clustering_error(predicted_biclustering, reference_biclustering, num_rows, num_cols)
    
    if size(predicted_biclustering,1)==0 && size(reference_biclustering,1)==0
        value = 1.0;
    elseif size(predicted_biclustering,1)==0 || size(reference_biclustering,1)==0
        value = 0.0;
    else
        union_size = calculate_size(predicted_biclustering, reference_biclustering,num_rows,num_cols,"union");
        dmax = calculate_dmax(predicted_biclustering,reference_biclustering);
        value = (dmax/union_size);
    end
    
end