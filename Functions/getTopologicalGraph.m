%{
This script is written as a part of TopoBARTMAP package to help in producing the Topological
 plot for a specified value of correlation level and number of Eigen values to be considered
 during prototype dimensionality reduction.
%}

function [meta_data] = getTopologicalGraph(TBM,data,correlationLevel,numDim)
	data = normalise(data);
    %generate needed variables
	num_prototypes_A = size(TBM.TAa.FA{1}.P,2);
    num_clusters_A = size(TBM.TAa.tC{1},1);
	num_prototypes_B = size(TBM.TAb.FA{1}.P,2);
	labels_B = TBM.TAb.FA{1}.L;
	labels_A = TBM.TAa.FA{1}.L;
	meta_data = struct();
	meta_data.correlationLevel = correlationLevel;
	meta_data.numReconstructionDim = numDim;
	PCA_ = struct();

	%% Assign the prototypes_B to specific prototypes_A based on greatest correlation
	meta_data.PrototypeAssociation = cell(num_prototypes_A,1);
    assigned = zeros(1,num_prototypes_B);
    %for corrlevel= meta_data.correlationLevel 
    for corrlevel=0.95:-0.05:meta_data.correlationLevel
	    for iterX =1:1:num_prototypes_A
			x_indices = find(labels_A == iterX);
			for iterY= 1:1:num_prototypes_B
				if assigned(1,iterY)~=1
	                y_indices = find(labels_B== iterY);
	            	Bicluster = data(y_indices,x_indices);
	            	num_correlations = 0;
	            	total_correlation = 0;
	            	for i =1:1:size(x_indices,2)-1
	            		mean_x = mean(Bicluster(:,i));
	            		terms_x = Bicluster(:,i) - mean_x;
	            		for j =i+1:1:size(x_indices,2)
	            			mean_y = mean(Bicluster(:,j));
	            		    terms_y = Bicluster(:,j) - mean_y;
	                        denominator = sqrt(sum(terms_x.*terms_x)*sum(terms_y.*terms_y));
	                        if denominator ~= 0
	                            %disp('test_3');
	                            corrl = sum(terms_x.*terms_y)/denominator;
	                        else
	                             %disp('test_4');
	                             corrl = 0;
	                        end
	                        total_correlation = total_correlation + corrl;
	                        num_correlations = num_correlations + 1;
	                    end
	                end
	             	if (total_correlation/num_correlations) > corrlevel
	                    meta_data.PrototypeAssociation{iterX,1} = [meta_data.PrototypeAssociation{iterX,1}, iterY];
	                    assigned(1,iterY) = 1;
	                end
	            end
	        end
	    end
	end
    
    %% t-SNE
    [PCA_.reconstructedData,loss] = tsne(TBM.TAb.FA{1}.P(1:38,:)','NumDimensions',meta_data.numReconstructionDim);
    
    %% Generate the graphs
    % identify nodes that need to included in the graph
    retainedIndices=[];
    for iter=1:1:size(meta_data.PrototypeAssociation,1)
        retainedIndices = [retainedIndices, meta_data.PrototypeAssociation{iter,1}];
    end
    disp(meta_data.PrototypeAssociation)
    retainedIndices = sort(retainedIndices);
    prototypeIndices = 1:1:num_prototypes_B;
    deletedIndices = setdiff(prototypeIndices,retainedIndices); %sorted array
    
    %establish connectivity between retained prototypes
    EdgeMatrixB = TBM.TAb.E{1};
    
    %delete the prototypes that are needed and create new EdgeMatrixB
    meta_data.Graph = struct();
    meta_data.Graph.EdgeMatrix = EdgeMatrixB(retainedIndices,retainedIndices);
    meta_data.Graph.retainedIndices = retainedIndices;
    meta_data.Graph.deletedIndices = deletedIndices;

    meta_data.Graph.G = graph(meta_data.Graph.EdgeMatrix);

    %assign the node locations
    meta_data.Graph.nodeLocations = PCA_.reconstructedData(retainedIndices,:);
    nodeColours = zeros(size(retainedIndices,2),3);
    nodeLabels = zeros(size(retainedIndices,2),1);
    %decide colours for each prototype in A
    colours = [];
    for iter_a = 1:1:num_clusters_A
    	if 1-(2*iter_a/num_clusters_A)>0
    		colour_r = (1-(2*(iter_a-1)/num_clusters_A));
            colour_g = iter_a*2/num_clusters_A;
            colour_b = 0; 
        else
        	colour_r = 0;
            colour_g = 2*(num_clusters_A - iter_a)/num_clusters_A;
            colour_b = (2*iter_a -num_clusters_A)/num_clusters_A;
        end
        colours(iter_a,:) = [colour_r,colour_g,colour_b];
    end
    for iter_a =1:1:num_prototypes_A
    	[~,ia,~] = intersect(retainedIndices,meta_data.PrototypeAssociation{iter_a});
        %identify to which cluster the current prototype belongs
         for iter_c = 1:1:size(TBM.TAa.tC{1},1) 
            if any(ismember(TBM.TAa.tC{1}{iter_c},iter_a))
                for iter_b = 1:1:size(ia,1)
                    nodeColours(ia(iter_b),:) = colours(iter_c,:);
                    nodeLabels(ia(iter_b),1) = retainedIndices(ia(iter_b));
                end
            end
        end
    end	
    meta_data.Graph.nodeColours = nodeColours;
    meta_data.Graph.nodeLabels = nodeLabels;
    meta_data.prototypeBColours = colours;
    
    %Add the plot
    figure()
    if numDim==2
        plot(meta_data.Graph.G,'NodeColor',nodeColours,'NodeLabel',nodeLabels,'XData',meta_data.Graph.nodeLocations(:,1),'YData',meta_data.Graph.nodeLocations(:,2), 'MarkerSize',6);
    end
    if numDim==3
        display('here')
        plot(meta_data.Graph.G,'NodeColor',nodeColours,'NodeLabel',nodeLabels,'XData',meta_data.Graph.nodeLocations(:,1),'YData',meta_data.Graph.nodeLocations(:,2),'ZData',meta_data.Graph.nodeLocations(:,3), 'MarkerSize',6);
    end
    title('Topological Graph for gene clusters')

%% end the function
end
