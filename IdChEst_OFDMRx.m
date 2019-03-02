function H = IdChEst_OFDMRx(rxG,txG)
numTx = size(txG,3);
numRx = size(rxG,3);
H = complex(zeros(size(rxG,1),size(rxG,2),numRx,numTx)); 

for i = 1:numTx
    for j=1:numRx
        ind = (1:12:length(txG)).';
        y = rxG(ind,1,j)./txG(ind,1,i);
        y = interp1(ind,y,1:length(rxG),'linear','extrap').';
        H(:,:,j,i) = y(:,ones(1,14));
    end
end