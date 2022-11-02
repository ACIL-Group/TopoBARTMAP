%{
Program Description: 
 This is a modified code for the TopoART based on the code Written on April 11

 Data and Function Dictionary
 P : Prototypes
 L : Labels
 V : vigilance parameter
 cI: Complement coded input
 a_: alpha
 b_: beta

 Authors: Leonardo Enzo Brito Da Silva, Raghu Yelugam
 Date: October 16, 2018
%}

%% Code:
classdef FuzzyART < handle
    
    properties(SetAccess = public)
        P =[]
        L =[]
        dI=[]
        V
        a_
        b_
    end
    
    methods
        function F = FuzzyART(v,a,b)
            F.V = v;
            F.a_ = a;
            F.b_ = b;
        end
        
        function [A,O] = Search(F,cI)
            if isempty(F.P)
                F.P(:,1) = cI;
                O =1;
                A =1;
            else
                A = fuzzyNorm1(cI,F.P)./(F.a_ + fuzzyNorm1(F.P));
                [M,O] = max(A);
                s_ = fuzzyNorm1(cI,F.P(:,O))/fuzzyNorm1(cI);
                if s_ >= F.V
                    F.P(:,O) = F.b_*min(F.P(:,O),cI) + (1-F.b_)*F.P(:,O);
                else
                    ss = true;
                    while ss
                        A(O) = 0;
                        [M,O] = max(A);
                        s_ = fuzzyNorm1(cI,F.P(:,O))/fuzzyNorm1(cI);
                        if s_ >= F.V
                            F.P(:,O) = F.b_*min(F.P(:,O),cI) + (1-F.b_)*F.P(:,O);
                            ss = false;
                        elseif max(A)==0
                            O = size(F.P,2)+1;
                            F.P(:,O)=cI;
                            ss =false;
                        end
                    end
                end
            end
        end
        
        function Train(F,cI) %other args)
            t=size(F.dI,2);
            tt=size(F.L,2);
            for x=1:1:size(cI,2)
               [A,O] = F.Search(cI(:,x));
                F.dI(:,t+x) = cI(:,x);
                F.L(tt+x) = O;
            end
        end
    end
end