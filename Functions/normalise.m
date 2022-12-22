% this file is for normalising the data sets presented to ART networks
%this fucntion assumes that the data is in
function data = normalise(data)
    tempMax=max(data,[],2);
    tempMin = min(data,[],2);
    data=(data -tempMin)./(tempMax-tempMin + 1e-16);
end
