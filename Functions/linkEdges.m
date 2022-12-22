% ------------------------------------------------------------------------------%
%  Description:								        %
% 								                %
%  This code is an implementation of Topological ARTMap. This method of topo- 	%
% logical clustering is used for identifying cluster from which data originates	%
% . Being an Adaptive Resonance Theory based model, TopoART is fast and can be	%
% implemented as an online algorithm. For this model we use online clustering,	%
% However, here we use the second network of the TopoART as the input to next	%
% module which learns features using unsupervised learning. 			%
% -----------------------------------------------------------------------------	%
% Author: Raghu Yelugam								%
% Date  : April 11,2018								%
% ----------------------------------------------------------------------------	%

function topoClusters = linkEdges(edgeMatrix)

%read the size of the edgeMatrix, and throw an error if its not a square matrix
	if size(edgeMatrix,1) ~= size(edgeMatrix,2)
		error('Error. \t Input matrix must be a square matrix');
	end	%endif
	topoClusters={};							%create an empty cell structure for the clusters record
	tmpSizeEdges=size(edgeMatrix,1);					%Square matrix
	for iter=1:tmpSizeEdges
		%this should render the indices which are greater than 0, ie, that are connected to particular node under consideration.
		connectedNodes = find(edgeMatrix(:,iter))';			%'
		connectedNodes = [iter connectedNodes];
        %disp(connectedNodes);
		if isempty(topoClusters)
			topoClusters{1,1}=connectedNodes;
		else
			mergeClusters=[];				%Record the clusters to be merged
			for iter_a = 1:size(connectedNodes,2);
				temp=connectedNodes(iter_a);
				for iter_b = 1:size(topoClusters,1)
					test=any(topoClusters{iter_b,1}==temp);
					if test
						sizeMerge= size(mergeClusters,2);
						mergeClusters(1,sizeMerge+1) = iter_b;
					end	%endif
				end	%endfor
			end	%endfor
			if isempty(mergeClusters)			%none of the elements of connectedNodes is found in any of the clusters
				temp=size(topoClusters,1);
				topoClusters{temp+1,1} = connectedNodes;	%create a new cluster for these nodes
			else
				mergeClusters = unique(mergeClusters);	%this will reduce repetitions and sort the array
				for iter_c = size(mergeClusters,2):-1:1
					if iter_c~=1
						connectedNodes=[connectedNodes topoClusters{mergeClusters(1,iter_c),1}];				
						topoClusters(mergeClusters(1,iter_c),:) =[];
					elseif iter_c == 1
						topoClusters{mergeClusters(1,iter_c),1}=unique([connectedNodes topoClusters{mergeClusters(1,iter_c),1}]);
					end	%endif
				end	%endfor
			end	%endif
		end	%endif
	end	%endfor
end	%endfunction
