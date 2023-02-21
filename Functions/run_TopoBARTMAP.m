function [y, TBM] = run_TopoBARTMAP(x, data, classes)

%% Data info
    [nSamples_B, nSamples_A] = size(data);

%% Parameters for TopoART module A
    Settings = struct();
    Settings.A = struct();
    Settings.A.nFA = 1;                     %number of Fuzzy ART Layers in TopoARTA, >1 if doing hierarchical.
    Settings.A.V = x(1);                    %vigilanceA
    Settings.A.al = 1e-3;                   %alphaA
    Settings.A.be = 1;                      %beta
    Settings.A.bSB = x(2);                  %beta Second Best
    Settings.A.ph = round(x(3));            %Phi for minimum number of members for not deleting a cluster
    Settings.A.ta = round(x(4)*nSamples_A);
    %     Settings.A.ta = max(Settings.A.ph, round(x(4)*nSamples_A)); %Number of cycles before the deletion happens in A

%% Parameters for TopoART module B
    Settings.B = struct();
    Settings.B.nFA = 1;                     %number of Fuzzy ART LAyers in TopoARTB, >1 if doing hierarchical
    Settings.B.V = x(5);                    %vigilanceB
    Settings.B.al = 1e-3;                   %Alpha
    Settings.B.be = 1;                      %beta
    Settings.B.bSB = x(6);                  %beta Second Best
    Settings.B.ph = round(x(7));            %Phi is number of minimum members in the cluster
    Settings.B.ta = round(x(8)*nSamples_B);
    %     Settings.B.ta = max(Settings.B.ph, round(x(8)*nSamples_B)); %Number of cycles before the deletion happens in B

    %% Other parameters
    Settings.cTh = x(9);                %Correlation threshold
    Settings.cFn = 'pearson';           %Correlation function
    Settings.sS = 0.01;                 %This can be used as parameter
    Settings.noiseLabel = 'true';       %should the noise be labelled or not?

%% Train
    TBM = TopoBARTMAP(Settings);
    TBM.Train(data');                   %because of the orientiation of input

%% Compute ARI
    perf = valid_external(classes, TBM.TAa.Label(1)');  % compute external CVI [Rand; AR; Jac; FM]

%% Return fitness value [Rand; AR; Jac; FM]
    y = -perf(4);
    

end
