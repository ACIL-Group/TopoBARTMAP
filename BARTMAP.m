%{
<<<<<<<<<<<<<<<<<<<<<<<<<<<<:: Program Description ::>>>>>>>>>>>>>>>>>>>>>>>>>>>
This code is a MATLAB implementation of BARTMAP. It is base on the MATLAB i-
mplementation of BARTMAP made by Islam Elnabarawy.

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>:: Authors ::<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                  Raghu Yelugam, Leonardo Enzo Brito Da Silva
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>:: Modification ::<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

Last Modified Date: 2019 December 01

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>::<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%}

classdef BARTMAP < handle
%please remove the label as handle if needed

    properties
        FAa
        FAa_Settings

        FAb
        FAb_Settings

        correlationThreshold
        correlationFunction
        matchTrackingStepSize

        crlChkWr_A
        crlChkWr_B
        noiseLabel
    end

    methods

        function BM = BARTMAP(stg)
            BM.FAa = FuzzyART(stg.A.V,stg.A.al,stg.A.be);

            BM.FAb = FuzzyART(stg.B.V,stg.B.al,stg.B.be);

            BM.correlationThreshold = stg.cTh;
            BM.correlationFunction = stg.cFn;
            BM.matchTrackingStepSize = stg.sS;

            %these parameters seem redundant
            BM.FAa_Settings = stg.A;
            BM.FAb_Settings = stg.B;
        end

        function Train(BM,DI)
            numFeatures = size(DI,2);               %number of Features
            numObservations = size(DI,1);           %number of Observations
            ttV = normalise(DI);
            cDI = [ttV; 1-ttV];                     %normalise the data

            %Train FuzzyART_B
%             disp('Training FuzzyART_B');
            BM.FAb.Train(cDI);

            %Train the FuzzyART_A
            ttV = normalise(DI');
            cDI = [ttV;1-ttV];                      %Transpose the data and normalise
            tV = BM.FAa.V;                          %This is recording the values of all the vigilance parameters for each layer of TopoART_A

            for yy = 1:1:numObservations

%                 disp(sprintf('presenting the sample - %d',yy));
                tPT = BM.FAa.P;                                 %tPT = temporary ProtoType
                tLbl =BM.FAa.L;                                 %tLbl = temporary Label array
                BM.FAa.V = tV;                                  %This step is to ensure that vigilance is reset
                %Start the search
                while true
                    tt = size(BM.FAa.P,2);
                    [A,O] = BM.FAa.Search(cDI(:,yy));           %A - activity vector, O - the highest activity prototype index
                    BM.FAa.L(1,size(BM.FAa.L,2)+1) = O;  %Add the new label to the list of Labels
                    if tt<O
%                         disp('found new prototype');
                        %A new prototype is created
                        break;    %while
                    else
                        cRL =[];
                        vSize = size(BM.FAb.P,2);
                        cRL = zeros(1,vSize);                   %cRL is correlation and is of size equal to the number of TopoClusters
                        uIX= find(BM.FAa.L ==O);                %get the Indices that are of same label
                        uIX = uIX(uIX~=size(BM.FAa.L,2));       %Remove the current Observations index

                        for jX=1:1:vSize                        %find all the units that correspond to this particular cluster
                            iIX = find(ismember(BM.FAb.L,jX));
                            bCL = cDI(iIX,uIX);
                            uD = cDI(iIX,yy);
                            cRL(jX) = BM.BiClusterCorrelation(bCL',uD',BM.correlationFunction);
                        end

                        if any(cRL > BM.correlationThreshold) || (BM.FAa.V==1)
                            %Learning the new pattern is covered by the FuzzyART
                            break;  %while true
                        else
                            BM.FAa.P = tPT;                     %restore the saved prototypes to restart the check
                            BM.FAa.L = tLbl;                    %restore the saved label array
                            BM.FAa.V = BM.FAa.V + BM.matchTrackingStepSize;
                            if BM.FAa.V>1
                                BM.FAa.V=1;
                            end
                        end     %if (~isempty(find(cRL > BM.correlationThreshold,1))) || (BM.FAa.V==1)
                    end     %if tt<O
                end %while true
            end  %for yy=1:1:numObservations
        end

        % there are other functions to be included. As of now the above code is getting debugged!

        function [r] = eval(TB,cDI,uID,iID,tL)       %define these parameters in the bank
            uC = BM.FAa.FA{tL}.L;
            iC = BM.FAb.FA{tL}.L;
            r=zeros(size(uID));
            for iy = 1:numel(uID)
                uX = uID(ix);
                iX = iID(ix);

                uIX = find(uC==uC(uX));
                if numel(uIX)>1
                    uIX = uIX(uIX ~=uX);
                end

                iIX = find(iC == iC(iX));
                bC = cDI(iIX,uIX);
                uD = cDI(iIX,uX);
                iIX = find(iIX == iX);

                rIX = BM.PredictRating(bC,uD,iIX,BM.correlationFunction);
                if rIX < 0
                    rIX = 0;
                elseif rIX >1
                    rIX = 1;
                end
                rIX = 5*rIX;
                r(iy) =rIX;
            end
        end

    end %end of methods

    methods(Static)

        function [Xmean, Ymean, coMask, coCount] = coMean(X, Y)
            %COMEAN Find the means of X and Y only for where X and Y are both nonzero

            coMask = X > 0 & Y > 0;
            Xmean = mean(X(coMask));
            Ymean = mean(Y(coMask));
            coCount = sum(coMask);
        end

        %this piece of code might have to be modified and then rerun the tests
        function [result,bCr] = BiClusterCorrelation(bicluster, user_data, correlation_function)
            %BICLUSTERCORR Compute correlation between a user and a user/item bicluster

            %% preallocate correlation vector
            num_users = size(bicluster, 1);
            bicorr = zeros(1, num_users);

            %% compute the correlation for each pair of users
            for ix = 1:num_users
                switch correlation_function
                  case {'mic','mine','mi'}
                      minestats = mine(user_data,bicluster(ix,:));
                      bicorr(ix) = minestats.mic;
                  otherwise
                    %% calculate the co-mean of this pair of users
                    [u1_mean, u2_mean, coMask] = TopoBARTMAP.coMean(user_data, bicluster(ix, :));
                    %% compute the terms for all the item values
                    terms1 = user_data(coMask) - u1_mean;
                    terms2 = bicluster(ix, coMask) - u2_mean;
                    %% compute the sums to find the user-pair correlation
                    numerator = sum(terms1 .* terms2);
                    root1 = sqrt(sum(terms1 .* terms1));
                    root2 = sqrt(sum(terms2 .* terms2));
                    if root1 == 0 || root2 == 0
                        r = 0;
                    else
                        r = numerator / (root1 * root2);
                    end

                    %% calculate all the correlation function values
                    r2 = r .^ 2;
                    n = sum(coMask);
                    if n < 2
                        % protect against sqrt of a negative number
                        t = 0;
                    else
                        if r2 == 1
                            % protect against division by zero
                            r2 = r2 - 1e-6;
                        end
                        t = r .* sqrt((n-2)./(1-r2));
                    end

                    %% return the appropriate correlation function value
                    switch correlation_function
                        case { 'r', 'pearson' }
                            bicorr(ix) = r;
                        case { 'abs-r', 'abs-pearson' }
                            bicorr(ix) = abs(r);
                        case { 'r2', 'cod', 'CoD' } % coefficient of determination (r^2)
                            bicorr(ix) = r2;
                        case { 't', 'student' } % t statistic correlation
                            bicorr(ix) = t;
                        case { 'abs-t', 'abs-student' } % absolute t-statistic value
                            bicorr(ix) = abs(t);
                    end
                end
            end

            %% compute the final correlation coefficient for bicluster
            result = mean(bicorr);

        end

        function [r] = PredictRating(bC, uD, iIx,correlationFunction)
            [~,bCr] = TopoBARTMAP.BiClusterCorrelation(bC,uD,correlationFunction);
            ix = find(uD>0);
            ix = ix(ix ~= iIx);
            if ~isempty(ix)
                uMu = mean(uD(ix));
            else
                uMu = 0.5;
            end
            uIx = bC(:,iIx)>0;
            iD = bC(uIx,iIx);
            cD = bCr(uIx);
            if ~isempty(iD) && nnz(cD)>0
                iMu = mean(iD);
                r = sum(cD.*(iD - iMu))/sum(abs(cD));
            else
                r=0;
            end
            r= uMu +r;

        end

    end
end %end of class
