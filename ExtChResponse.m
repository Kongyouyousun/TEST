function hD=ExtChResponse(chEst, idx_data)
%#codegen
numTx = size(chEst,4);
numRx = size(chEst,3);
hD=complex(zeros(length(idx_data),numRx,numTx));
for n=1:numRx
    for m=1:numTx
        tmp=chEst(:,:,n,m);
        tmp = tmp(:);
        hD(:,n,m)=tmp(idx_data(:,n));
    end
end
