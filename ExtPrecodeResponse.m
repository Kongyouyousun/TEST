function [hD,hD_all]=ExtPrecodeResponse(chEst, idx_data)
%#codegen
numTx = size(chEst,4);
numRx = size(chEst,3);
hD=complex(zeros(numel(idx_data),numRx,numTx));
hD_all=zeros(size(chEst,1)*size(chEst,2),numRx,numTx);  % 18480x2x32
for n=1:numRx
    for m=1:numTx
        tmp=chEst(:,:,n,m);
        tmp = tmp(:);
        hD(:,n,m)=tmp(idx_data);
        hD_all(:,n,m)=tmp;
    end
end
% %需要把UL_precoding_temp中第4,11符号去掉才能对应15840做预编码
% delete=[size(chEst,1)*3+(1:size(chEst,1))  size(chEst,1)*10+(1:size(chEst,1))];  %12*110=1320, 12*8=96
% hD_all(delete,:,:)=[];

%  hD=complex(zeros(numel(idx_data),numTx));
%   for m=1:numTx
%         tmp=chEst(:,:,m);
%         tmp = tmp(:);
%         hD(:,m)=tmp(idx_data);
%     end