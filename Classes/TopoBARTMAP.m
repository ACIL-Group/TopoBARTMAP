%{
<<<<<<<<<<<<<<<<<<<<<<<<<<<<:: Program Description ::>>>>>>>>>>>>>>>>>>>>>>>>>>>
This code is a MATLAB implementation of TopoBARTMAP. It is base on the MATLAB i-
mplementation of BARTMAP made by Islam Elnabarawy.

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>:: Authors ::<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                  Raghu Yelugam, Leonardo Enzo Brito Da Silva
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>:: Modification ::<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

A cleaner TopoBARTMAP code

Last Modified Date: 22 August 2019

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>::<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%}

classdef TopoBARTMAP < handle

    properties
        TAa
        TAa_Settings

        TAb
        TAb_Settings

        correlationThreshold
        correlationFunction
        matchTrackingStepSize

        crlChkWr_A
        crlChkWr_B
        noiseLabel
    end

    methods

        function TB = TopoBARTMAP(stg)
            TB.TAa = TopoART(stg.A.nFA,stg.A.V,stg.A.al,stg.A.be,stg.A.bSB,stg.A.ph,stg.A.ta);

            TB.TAb = TopoART(stg.B.nFA,stg.B.V,stg.B.al,stg.B.be,stg.B.bSB,stg.B.ph,stg.B.ta);

            TB.correlationThreshold = stg.cTh;
            TB.correlationFunction = stg.cFn;
            TB.matchTrackingStepSize = stg.sS;

            TB.noiseLabel = stg.noiseLabel;

            %these parameters seem redundant
            TB.TAa_Settings = stg.A;
            TB.TAb_Settings = stg.B;
        end

        function Train(TB,DI)
            numFeatures = size(DI,2);               %number of Features
            numObservations = size(DI,1);           %number of Observations
            ttV = normalise(DI);                    %normalise the data with respect to Features
            cDI = [ttV; 1-ttV];                     %Complement code

            %Train TopoART_B
            %disp(sprintf('Training TopoART_B with %d number of samples', numFeatures));
            for xx = 1:1:numFeatures
                %disp(sprintf('presenting the sample - %d',xx));
                TB.TAb.Train(cDI(:,xx),xx);
            end

            %link the unlinked edges
            %disp('Linking ')
            for xx =1:1:TB.TAb_Settings.nFA
                TB.TAb.linkEdges(xx);
            end

            for xx =1:1:TB.TAb_Settings.nFA
                if isempty(TB.TAb.FA{xx}.P)
                    TB.TAb.FA{xx}.P = ones(2*numObservations, 1);
                end
            end

            if TB.noiseLabel
                TB.LabelNoise(cDI,'b');
            end

            %Train the TopoART_A
            ttV = normalise(DI');               %' to Transpose the data to for Clustering on Observations
            cDI = [ttV;1-ttV];                  % and normalise
            tV = [];                            %Record vigilance parameters for each layer of TopoART_A
            for xx = 1:1:TB.TAa_Settings.nFA
                tV(xx) = TB.TAa.V(xx);
            end

            for yy = 1:1:numObservations

%                 disp(sprintf('presenting the sample - %d',yy));
                %set the propagate to Layer variable so that first layer training is done
                pL = true;
                for xx = 1:1:TB.TAa_Settings.nFA
                if pL
                    tPT = TB.TAa.FA{xx}.P;                                      %tPT = temporary ProtoType
                    tLbl =TB.TAa.FA{xx}.L;                                      %tLbl = temporary Label array
                    TB.TAa.FA{xx}.V=TB.TAa.V(xx);                               %This step is to ensure that vigilance is reset
                    TB.TAa.cN{xx} = TB.TAa.cN{xx} +1;                           %Increase the cycle number to account for this cycle
                    %Start the search
                    while true
                        [A,O] = TB.TAa.FA{xx}.Search(cDI(:,yy));                %A - activity vector, O - the highest activity prototype index
                        tt = size(TB.TAa.c{xx},2);                              %record the size of prototype activity counter
                        TB.TAa.FA{xx}.L(1,size(TB.TAa.FA{xx}.L,2)+1) = O;       %Add the new label to the list of Labels
                        if tt<O
                            %disp('found new prototype');
                            %A new prototype is created
                            TB.TAa.c{xx} = [TB.TAa.c{xx} 1];
                            if isempty(TB.TAa.E{xx})                            % If it is the first prototype there are no edges yet
                                TB.TAa.E{xx} =0;
                                break;    %while
                            else
                                TB.TAa.E{xx} = [TB.TAa.E{xx} zeros(size(TB.TAa.E{xx},1),1); zeros(1,size(TB.TAa.E{xx},2)) 0];
                                break;  %while
                            end
                        else
                            cRL =[];
                            vSize = size(TB.TAb.FA{xx}.P,2);
                            cRL = zeros(1,vSize);                               %cRL is correlation and is of size equal to the number of TopoClusters
                            uIX= find(TB.TAa.FA{xx}.L ==O);                     %get the Indices that are of same label
                            uIX = uIX(uIX~=size(TB.TAa.FA{xx}.L,2));            %remove the current observation from the list
                            for jX=1:1:vSize                                    %get all the members that correspond to this particular cluster
                                iIX = find(ismember(TB.TAb.FA{xx}.L,jX));
                                bCL = cDI(iIX,uIX);
                                uD = cDI(iIX,yy);
                                cRL(jX) = TB.BiClusterCorrelation(bCL',uD',TB.correlationFunction);
                            end

                            if any(cRL > TB.correlationThreshold) || (TB.TAa.FA{xx}.V==1)

                                %increase the counter corresponding to the O
                                %tCTC =find(cRL>TB.correlationThreshold);       %temporary Correlated TopoCluster
                                TB.TAa.c{xx}(1,O) = TB.TAa.c{xx}(1,O)+1;
                                %Learning the new pattern is covered by the FuzzyART
                                if xx < TB.TAa.nFA
                                    if TB.TAa.c{xx}(1,O) < TB.TAa.p_(xx)
                                        pL = false;
                                    end
                                end

                                if tt > 1
                                    ii = O;
                                    s_ = fuzzyNorm1(cDI(:,yy),TB.TAa.FA{xx}.P)/fuzzyNorm1(cDI(:,yy));
                                    A(O)= 0;
                                    s_(O)=0;
                                    sGV = [s_>=TB.TAa.V(xx)];   %[s_>= TB.TAa.FA{xx}.V];
                                    while any(sGV)
%                                         disp('Searching for the second best')
                                        [tMax,O] = max(A.*sGV);
                                        cRLSB =[];
                                        vSize = size(TB.TAb.FA{xx}.P,2);
                                        cRLSB = zeros(1,vSize);                             %cRLSB is correlation and is of size equal to the number of TopoClusters
                                        uIX= find(TB.TAa.FA{xx}.L ==O);                     %get the Indices that are of same label
                                        uIX = uIX(uIX~=size(TB.TAa.FA{xx}.L,2));            %remove the current Observation index
                                        for jX=1:1:vSize                                    %get all the members that correspond to this particular cluster
                                            iIX = find(ismember(TB.TAb.FA{xx}.L,jX));
                                            bCL = cDI(iIX,uIX);
                                            uD = cDI(iIX,yy);
                                            cRLSB(jX) = TB.BiClusterCorrelation(bCL',uD',TB.correlationFunction);
                                        end

                                        if any(cRLSB > TB.correlationThreshold)
                                            TB.TAa.FA{xx}.P(:,O) = TB.TAa.bSB_(xx)*min(TB.TAa.FA{xx}.P(:,O),cDI(:,yy)) + (1-TB.TAa.bSB_(xx))*TB.TAa.FA{xx}.P(:,O);
                                            TB.TAa.E{xx}(ii,O) = 1;
                                            TB.TAa.E{xx}(O,ii) = 1;
                                            break;
                                        end
                                        if tMax == 0
                                            break;
                                        end
                                        sGV(O) = 0;
                                        A(O) = 0;
                                    end     %while any(sGV)
                                end % if tt > 1
                                break;  %while true
                            else
                                TB.TAa.FA{xx}.P = tPT;                      %restore the saved prototypes to restart the check
                                TB.TAa.FA{xx}.L = tLbl;                     %restore the saved label array
                                TB.TAa.FA{xx}.V = TB.TAa.FA{xx}.V + TB.matchTrackingStepSize;
                                if TB.TAa.FA{xx}.V>1
                                    TB.TAa.FA{xx}.V=1;
                                end
                            end     %if (~isempty(find(cRL > TB.correlationThreshold,1))) || (TB.TAa.FA{xx}.V==1)
                        end     %if tt<O
                    end %while
                    %check for the deletion
                    if mod(TB.TAa.cN{xx},TB.TAa.t_(xx))==0
                        tC=size(TB.TAa.FA{xx}.P,2);
                        while ~(tC<1)

                            if TB.TAa.c{xx}(1,tC)>TB.TAa.p_(xx)
                                tC=tC-1;                                    %retain the current prototype
                            else
                                tt= find(TB.TAa.FA{xx}.L==tC);              %find members of current prototype
                                %add members of current prototype to the bucket list
                                TB.TAa.bI{xx}(:,size(TB.TAa.bI{xx},2)+1:size(TB.TAa.bI{xx},2)+size(tt,2))=tt;
                                tt = find(TB.TAa.FA{xx}.L>tC);              %find members of other prototypes that are given a label greater than current prototype
                                TB.TAa.FA{xx}.L(TB.TAa.FA{xx}.L==tC)=-1;    %Label the members of current prototype as noise
                                %adjust the labels of other prototypes
                                for xy = 1:1:size(tt,2)
                                    TB.TAa.FA{xx}.L(tt(xy))= TB.TAa.FA{xx}.L(tt(xy)) -1;
                                end
                                % deletion ensues
                                TB.TAa.FA{xx}.P(:,tC) =[];
                                TB.TAa.c{xx}(:,tC) =[];
                                TB.TAa.E{xx}(:,tC) =[];
                                TB.TAa.E{xx}(tC,:) =[];
                                tC=tC-1;
                            end
                        end
                    end
                    %if not propagated to later layers
                    if ~pL
                        %add to respective noise lists
                        for y = xx +1:1:TB.TAa.nFA
                            tt = size(TB.TAa.bI{y},2);
                            TB.TAa.bI{y}(1,tt+1) = yy;
                        end
                        break;
                    end
                end
                end
            end

            %add uncommitted prototype if learning deletes all prototypes
            for xx =1:1:TB.TAb_Settings.nFA
                if isempty(TB.TAa.FA{xx}.P)
                    TB.TAa.FA{xx}.P = ones(2*numFeatures, 1);
                    TB.TAa.E{xx} =0;
                end
            end

            TB.TAa.linkEdges(xx);

            if TB.noiseLabel
                TB.LabelNoise(cDI,'a');
            end

        end

        function LabelNoise(TB,dI,l)
            switch l
                case {'a', 'TAa'}
                    for xx = 1:1:TB.TAa_Settings.nFA
                        tt = TB.TAa.bI{xx};
                        for yy = 1:1:size(tt,2)
                            O = TB.TAa.AssignPrototype(dI(:,tt(yy)),xx);
                            TB.TAa.FA{xx}.L(tt(yy)) = O;
                        end
                    end

                case {'b','TAb'}
                    for xx = 1:1:TB.TAb_Settings.nFA
                        tt = TB.TAb.bI{xx};
                        for yy = 1:1:size(tt,2)
                            O = TB.TAb.AssignPrototype(dI(:,tt(yy)),xx);
                            TB.TAb.FA{xx}.L(tt(yy)) = O;
                        end
                    end
            end
        end

    end %end of methods

    methods(Static)

        %this piece of code might have to be modified and then rerun the tests
       function [Xmean, Ymean, coMask, coCount] = coMean(X, Y)
            %COMEAN Find the means of X and Y only for where X and Y are both nonzero
            coMask = X > 0 & Y > 0;
            Xmean = mean(X(coMask));
            Ymean = mean(Y(coMask));
            coCount = sum(coMask);
        end
        
        function [result] = BiClusterCorrelation(bicluster, user_data, correlation_function)
            %BICLUSTERCORR Compute correlation between a user and a user/item bicluster

            %% preallocate correlation vector
            num_users = size(bicluster, 1);
            bicorr = zeros(1, num_users);

            %% compute the correlation for each pair of users
            for ix = 1:num_users
                switch correlation_function
                    case{'mic','mine','Mine','MIC','mi'}
                        [u1_mean, u2_mean, coMask] = TopoBARTMAP.coMean(user_data, bicluster(ix, :));
                        terms1 = user_data(coMask);
                        terms2 = bicluster(ix,coMask);
                        
                        if size(terms1,2)>1 && size(terms2,2)>1
                            %disp(size(terms1));
                            %disp(size(terms2));
                            minestats = mine(terms1,terms2);
                            bicorr(ix) = minestats.mic;
                        else
                            bicorr(ix) = 0;
                        end
                    case{'pearson','Pearson'}
                        %% calculate the co-mean of this pair of users
                        u1_mean = mean(user_data);
                        u2_mean = mean(bicluster(ix,:));
                        %% compute the terms for all the item values
                        terms1 = user_data - u1_mean;
                        terms2 = bicluster(ix, :) - u2_mean;
                        %% compute the sums to find the user-pair correlation
                        numerator = sum(terms1 .* terms2);
                        root1 = sqrt(sum(terms1 .* terms1));
                        root2 = sqrt(sum(terms2 .* terms2));
                        if root1 == 0 || root2 == 0
                            bicorr(ix) = 0;
                        else
                            bicorr(ix) = numerator / (root1 * root2);
                        end
                    case{'modifiedPearson'}
                        [u1_mean, u2_mean, coMask] = TopoBARTMAP.coMean(user_data, bicluster(ix, :));
                        %% compute the terms for all the item values
                        terms1 = user_data(coMask) - u1_mean;
                        terms2 = bicluster(ix, coMask) - u2_mean;
                        %% compute the sums to find the user-pair correlation
                        numerator = sum(terms1 .* terms2);
                        root1 = sqrt(sum(terms1 .* terms1));
                        root2 = sqrt(sum(terms2 .* terms2));
                        if root1 == 0 || root2 == 0
                            bicorr(ix) = 0;
                        else
                            bicorr(ix) = numerator / (root1 * root2);
                        end
                   
                end
            end

            %% compute the final correlation coefficient for bicluster
            result = mean(bicorr);
        end
    end
end %end of class
