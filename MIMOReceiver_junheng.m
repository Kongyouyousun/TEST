function y = MIMOReceiver_junheng(pdschRx, G_all_in, sigma2)
UE_NUM=16;
% noisFac = diag(nVar);%像王东明老师在程序的加噪声的方法.
numData = size(pdschRx, 1);
y = complex(zeros(size(pdschRx)));
Nt=2;
DetSigal=zeros(numData,32);
for n = 1:numData
    G_all=G_all_in(:,:,n);
    MMSE_Filter=zeros(2,2,16);
    rho=zeros(2,16);
    
    for k=1:UE_NUM
        Tmp_inv = inv( G_all((k-1)*Nt+(1:Nt),:)* G_all((k-1)*Nt+(1:Nt),:)' +2*sigma2*eye(Nt));
        MMSE_Filter(:, :, k) =G_all((k-1)*Nt+(1:Nt),(k-1)*Nt+(1:Nt))' * Tmp_inv;
        rho(:,k)=diag(MMSE_Filter(:,:,k)*G_all((k-1)*Nt+(1:Nt),(k-1)*Nt+(1:Nt)));
        tmp_pdschRx=squeeze(pdschRx(n,(k-1)*Nt+(1:Nt))).';
        DetSigal(n,(k-1)*Nt+(1:Nt))=(MMSE_Filter(:, :, k)*tmp_pdschRx./rho(:,k)).';
    end
end
y = DetSigal;
end

  