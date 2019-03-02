function hD=ExtractResources(idx_data,chEst)
%#codegen
numRx = size(chEst,3);
hD=complex(zeros(length(idx_data),numRx));
for n=1:numRx
        tmp=chEst(:,:,n);
        tmp = tmp(:);
        hD(:,n)=tmp(idx_data(:,n));
end