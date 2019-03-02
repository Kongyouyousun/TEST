function  chEst = ChEst_freq_mmseforCSIrs(EsN0, NumDLRB, sizeFFT,Tap_delay,grid,txgrid,ind)
%%  
%calculate frequency mmse matrices for subcarriers with CSI-RS
% Author:
% Haiyan Ding 2011-11-9
chEst = complex(zeros([size(grid) 32]));
sizeFFT = double(sizeFFT);
for numRx = 1:32
    rxGrid = grid(:,:,numRx);
    for numTx =1:32
        txGrid = txgrid(:,:,numTx);
        Tx = lteExtractResources(ind(:,numTx),txGrid);
        Rx = lteExtractResources(ind(:,numTx),rxGrid);
        hp = Rx./Tx;
        hp = reshape(hp,length(hp)/4,4);
        numRE = NumDLRB*12;
        numCSIRB=size(ind,1)/4;
        idxh = (sizeFFT-NumDLRB*12)/2+1 : (sizeFFT+NumDLRB*12)/2;
        locCSI = ind(1:length(ind)/4,numTx)-size(grid,1)*7;
        idxhp1= locCSI+(sizeFFT-NumDLRB*12)/2;
        m_n = idxhp1*ones(1,length(idxhp1)) - (idxhp1*ones(1,length(idxhp1))).';
        m_n1 = idxh.'*ones(1,length(idxhp1)) - (idxhp1*ones(1,length(idxh))).';
        Rhphp = double(zeros(numCSIRB, numCSIRB));
        Rhhp1 = double(zeros(numRE, numCSIRB));
        %% calculate Rhh with mismatched profile
        L = max(ceil(Tap_delay(end))) + 1;
        alfa = 1;
        l = 0:L-1;
        g = (alfa - 1)/(L-1) .* l + 1;
        
        for i=1:length(g)
            if g(i)<0
                g(i)=0;
            end
        end
        PathGain = g / sqrt(sum(g.^2));
        tap = 0:L-1;
        for k = 1 : length(tap)
            Rhphp = Rhphp+ PathGain(k).^2 * exp(-2j*pi*(m_n )*tap(k)/sizeFFT);
            Rhhp1 = Rhhp1+ PathGain(k).^2 * exp(-2j*pi*(m_n1)*tap(k)/sizeFFT);
        end
        mmseA1 = zeros(numCSIRB, numCSIRB );
        mmseB1 = zeros( numRE, numCSIRB );
        Mat = Rhphp+1/EsN0*eye(numCSIRB);%Ôö
        mmseA1 = Rhphp/Mat;
        mmseB1 = Rhhp1/Mat;
        MMSEBforCSI = mmseB1;
        HpCsi=MMSEBforCSI*hp;
        y=mean(HpCsi,2);
        hD = y(:,ones(1,14));
        chEst(:,:,numRx,numTx)=hD;
    end
end
