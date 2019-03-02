function G_all=ExtGallResponse(tmp_G_all, idx_data)
% G_all=ExtGallResponse(tmp_G_all,pdschIndices);
% tmp_G_all:2x32x1320x14 idx_data:8360¸öÎ»ÖÃ
numTx = size(tmp_G_all,1);
numRx = size(tmp_G_all,2);
G_all=zeros(numTx,numRx,(numel(idx_data)));

for n=1:numTx
    for m=1:numRx
        tmp=tmp_G_all(n,m,:,:);
        tmp = tmp(:);
        G_all(n,m,:)=tmp(idx_data);
        
    end
end
G_all;