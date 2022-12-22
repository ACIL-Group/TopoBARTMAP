%{
This function is a part of implementation of Clustering Error (CE) an external
bicluster evaluation metric introduced by Patrikainen and Meila (2006).

This implimentation is based on Biclustlib, by Victor Alexandre Padilha

Param: biclusters: cell array with
{[Bi rows],[Bi columns]; [Bj rows],[Bj columns]} form

%}

function output = bic2sets(biclusters)
    n_biclusters = size(biclusters,)

end