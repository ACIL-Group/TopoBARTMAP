%{
This function is a part of implementation of Clustering Error (CE) an external
bicluster evaluation metric introduced by Patrikainen and Meila (2006).

This implimentation is based on Biclustlib, by Victor Alexandre Padilha

Param: predicted, reference: cell array with
{[Bi rows],[Bi columns]; [Bj rows],[Bj columns]} form

%}

function value = calculate_dmax(predicted,reference)
    n_pred_bic = size(predicted,1);
    n_ref_bic = size(reference,1);
    cost_matrix = zeros(n_pred_bic,n_ref_bic);

    for iter_pred=1:1:n_pred_bic
        for iter_ref=1:1:n_ref_bic
            row_intersection = size(intersect(predicted{iter_pred,1},reference{iter_ref,1}),2);
            col_intersection = size(intersect(predicted{iter_pred,2},reference{iter_ref,2}),2);
            cost = row_intersection*col_intersection;
            cost_matrix(iter_pred,iter_ref) = cost;
        end
    end
    
    costofnonassignment = 0.2;
    [assignments, unassignedTracks, unassignedDetections] = assignmunkres(-cost_matrix,costofnonassignment);
    value = 0;

    for itr=1:1:size(assignments,1)
        value = value+cost_matrix(assignments(itr,1),assignments(itr,2));
    end
    
end