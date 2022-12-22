%{
This function is written to implement Calinski-Harabasz index for
clustering validity.
%}

function ch_value = CalinskiHarabaszIndex(data,labels)
    if size(labels,1) ~= size(data,2)
        error("Number of labels does not match the number of samples in data");
    end
    clusters = unique(labels);
    
    %% Compute cluster centroid and data centroid (barycentre)
    barycentre = mean(data,2);
    centre_k = zeros(size(data,1) ,size(clusters,1));
    nSamplesPerCluster = 0*clusters;

    %% Compute intra-cluster dispersion
    WGSS_k = zeros(size(clusters));

    for iter=1:1:size(clusters,1)
        samples_k = data(:,find(labels==clusters(iter)));
        centre_k(:,iter) = mean(samples_k,2);
        nSamplesPerCluster(iter) = sum(labels==clusters(iter));
        WGSS_k(iter)= sum(eucDist(samples_k, centre_k(:,iter)).^2);
    end

    %% Compute inter-cluster dispersion
    BGSS = sum(nSamplesPerCluster'.*(eucDist(centre_k,barycentre).^2));

    WGSS = sum(WGSS_k);

    %% Compute Calinski-Harabasz score
    ch_value = (BGSS*(size(labels,1)-size(clusters,1)))/(WGSS*(size(clusters,1)-1));
end
