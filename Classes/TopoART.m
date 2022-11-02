%{
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>: Program Description :<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 This is a modified code for the TopoART based on the code Written on April 11. We assume
 the data is not going to change its  order while and after presentation and we record the
 indices as they follow during the training. Accordingly a bucket variable is created  to
  register the indices of the data samples as they are presented and the record the corre-
 sponding label in a different array, 'L'. Please refer to the Data Dictionary for details.

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>: Data Dictionary :<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 :   These are cell data structures  :
 FA : FuzzyART modules
 E : Adjacency matrix to represent edges between the prototypes
 tC : TopoClusters
 c : counter
 bI : Bucket Index
 :  These are expected be arrays arrays   :
 dI : data input
 al_ : Alpha
 be_ : beta
 ph_ : Phi
 ta_ : Tau
 bSB_ : Beta second best
 nFA : number of FuzzyART modules needed, ie the number of heirarchical levels of clusters
 V : vigilance parameter
 cN : Cycle number
 II : Index of the input cI
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>::<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 Authors: Leonardo Enzo Brito Da Silva, Raghu Yelugam
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>::<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                            Last Modified Date: 17 July 2019

                                    : Modification :
  The bI is used to record all the indices instead of the bucketed indices. This makes sure
that every label that is set to -1 while bucketing or not propagated in to the upper layers
is recorded along the sample labels which are registered while training the TopoART modules.

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>::<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%}

classdef TopoART < handle

    properties(SetAccess = public)

        FA = {}
        E = {}
        tC = {}
        c ={}
        V=[]
        dI ={}
        bI = {}
        a_
        b_
        p_
        t_
        bSB_
        nFA
        cN ={}
    end

    methods
        function TA = TopoART(nFA,V,al,be,bSB,ph,ta)       %add other variables that are need to be initialised
            TA.nFA = nFA;
            TA.V(1) = V;
            TA.a_ = al;
            TA.b_ = be;
            TA.p_ = ph;
            TA.bSB_ = bSB;
            TA.t_ = ta;

            for x= 1:1:TA.nFA
                if x>1
                    TA.V(x) = (V(x-1) +1)/2;
                end
                TA.FA{x} = FuzzyART(TA.V(x),TA.a_(x),TA.b_(x));
                TA.c{x} = [];
                TA.tC{x} ={};
                TA.E{x} = [];
                TA.cN{x} = 0;
                TA.bI{x} = [];
                TA.dI{x} = [];
            end
            %if there are anyother thing to initialise
        end

        function Train(TA,cI,II)
            pL = true;
            %disp(sprintf('input = %d',II));                                 %propagate to the next layer
            for x = 1:1:TA.nFA
                TA.cN{x} = TA.cN{x} + 1;
                [A,O] = TA.FA{x}.Search(cI);
                s_ = fuzzyNorm1(cI,TA.FA{x}.P)/fuzzyNorm1(cI);
                tt = size(TA.c{x},2);
                TA.FA{x}.L(size(TA.FA{x}.L,2)+1) = O;
                TA.dI{x}(size(TA.FA{x}.L,2)) = II;
                if tt<O
                    TA.c{x}(1,O) = 1;
                    if size(TA.E{x},2)==0
                        TA.E{x} = 0;
                    else
                        TA.E{x} = [TA.E{x} zeros(size(TA.E{x},1),1);zeros(1,size(TA.E{x},2)) 0];
                    end
                else
                    if O>size(TA.c{x},2)
                        TA.c{x}(1,O) =1;
                    else
                        TA.c{x}(1,O) = TA.c{x}(1,O) +1;
                    end
                    if x < TA.nFA
                        if TA.c{x}(1,O) < TA.p_(x)
                            pL = false;
                        end
                    end
                    if tt > 1
                        ii = O;
                        %tic
                        A(O)= 0;
                        s_(O) = 0;
                        temp = (s_>=TA.V(x));
                        %{
                            s_>=TA.V(x) finds only the prototypes that qualify the vigilance test. Since this includes
                            the first winner, which has the highest activity amongst all, we set the first winner sim-
                            -ilarity and activity to 0. This done, finding the next highest activity prototype finds
                            second winning prototype.
                        %}
                        if any(temp)
                            [~,O] = max(A.*temp);
                            TA.FA{x}.P(:,O) = TA.bSB_(x)*min(TA.FA{x}.P(:,O),cI) + (1-TA.bSB_(x))*TA.FA{x}.P(:,O);
                            TA.E{x}(ii,O) = 1;
                            TA.E{x}(O,ii) = 1;
                        end
                    end
                end
                if (mod(TA.cN{x},TA.t_(x))==0) && (TA.cN{x} ~=0)
                    tC=size(TA.FA{x}.P,2);
                    while ~(tC<1)
                        %disp(size(TA.c{x}))
                        if TA.c{x}(1,tC)>TA.p_(x)
                            tC=tC-1;
                        else
                            tt = find(TA.FA{x}.L==tC);
                            TA.bI{x}(:,size(TA.bI{x},2)+1:size(TA.bI{x},2)+size(tt,2))=tt;
                            TA.FA{x}.L(tt) = -1;
                            tt = find(TA.FA{x}.L>tC);
                            for xx = 1:1:size(tt,2)
                                TA.FA{x}.L(tt(xx))= TA.FA{x}.L(tt(xx)) -1;
                            end
                            TA.FA{x}.P(:,tC) =[];
                            TA.c{x}(:,tC) =[];
                            TA.E{x}(:,tC) =[];
                            TA.E{x}(tC,:) =[];
                            tC=tC-1;
                        end
                    end

                    %TA.linkEdges(TA.E{x});
                end
                if ~pL
                    for y = x+1:1:TA.nFA
                        TA.bI{x}(:,size(TA.bI{x},2)+1)=size(TA.FA{x}.L,2) + 1;
                        TA.dI{y}(1,size(TA.FA{x}.L,2)+1) = II;
                        TA.FA{x}.L(size(TA.FA{x}.L,2)+1) = -1;
                    end
                    break;    %break the for loop
                end
            end
        end

        %this function is for classifying the data samples thrown into bucket to the corresponding cluster
        function O_ = AssignPrototype(TA,cI,x)
            A_ = 1 - fuzzyNorm1(min(TA.FA{x}.P,cI) - TA.FA{x}.P)/fuzzyNorm1(cI);
            [M,I] = max(A_);
            O_ = I;
        end

        % we need one function to record the labels, but we have it as an external function which needs to be integrated to this class
        function O_ = Label(TA,x)
            O_=[];
            for y=1:1:size(TA.FA{x}.L,2)
                cFC = false;                                               %cFC = cluster Not found
                for ii = 1:1:size(TA.tC{x},1)
                    t = TA.tC{x}{ii,1};
                    if ismember(TA.FA{x}.L(y),t)
                        O_(y) = ii;
                        cFC = true;
                        break;
                    end
                end
                if ~ cFC
                    O_(y) = size(TA.tC{x},1) +1;
                end
            end
        end

        function linkEdges(TA,x)
        %read the size of the edgeMatrix, and throw an error if its not a square matrix
            if size(TA.E{x},1) ~= size(TA.E{x},2)
                error('Error. \t Input matrix must be a square matrix');
            end %endif
            topoClusters={};                            %create an empty cell structure for the clusters record
            tmpSizeEdges=size(TA.E{x},1);                    %Square matrix
            for iter=1:tmpSizeEdges
                %this should render the indices which are greater than 0, ie, that are connected to particular node under consideration.
                connectedNodes = find(TA.E{x}(:,iter))';         %'
                connectedNodes = [iter connectedNodes];
                %disp(connectedNodes);
                if isempty(topoClusters)
                    topoClusters{1,1}=connectedNodes;
                else
                    mergeClusters=[];               %Record the clusters to be merged
                    for iter_a = 1:size(connectedNodes,2);
                        temp=connectedNodes(iter_a);
                        for iter_b = 1:size(topoClusters,1)
                            test=any(topoClusters{iter_b,1}==temp);
                            if test
                                sizeMerge= size(mergeClusters,2);
                                mergeClusters(1,sizeMerge+1) = iter_b;
                            end %endif
                        end %endfor
                    end %endfor
                    if isempty(mergeClusters)           %none of the elements of connectedNodes is found in any of the clusters
                        temp=size(topoClusters,1);
                        topoClusters{temp+1,1} = connectedNodes;    %create a new cluster for these nodes
                    else
                        mergeClusters = unique(mergeClusters);  %this will reduce repetitions and sort the array
                        for iter_c = size(mergeClusters,2):-1:1
                            if iter_c~=1
                                connectedNodes=[connectedNodes topoClusters{mergeClusters(1,iter_c),1}];
                                topoClusters(mergeClusters(1,iter_c),:) =[];
                            elseif iter_c == 1
                                topoClusters{mergeClusters(1,iter_c),1}=unique([connectedNodes topoClusters{mergeClusters(1,iter_c),1}]);
                            end %endif
                        end %endfor
                    end %endif
                end %endif
            end %endfor
            TA.tC{x} = topoClusters;
        end %endfunction

        function PrintClusters(TA,x)
            t = size(TA.FA{x}.P,1)/2;
            U = TA.FA{x}.P(1:t,:);
            V = TA.FA{x}.P(t+1:2*t,:);
            vC = 1-V;
            %disp(size(TA.tC{x},1));
            for ii = 1:size(TA.tC{x},1)

                t0 = TA.tC{x}{ii,:};
                t1 = rand(1,3);
                for ll = 1:size(t0,2)
                    disp(ll)
                    tP = genCoordinate(U(:,t0(ll)),vC(:,t0(ll)));
                    %disp(ll);
                    hold on;
                    for it = 1:size(tP,1)-1
                        %disp(it)
                        for ik = it+1:size(tP,1)
                            t2 = [tP(it,:);tP(ik,:)];
                            %disp(t2);
                            if size(t2,2)==2
                                line(t2(:,1),t2(:,2),'Color',t1);
                            elseif size(t2,2)==3
                                line(t2(:,1),t2(:,2),t2(:,3),'Color',t1);
                            end
                        end
                    end
                end
            end
        end
    end
end
