lenGrid = size(rxGrid_7500,1)/5;
        pdschRx = [];
        pdschTx = [];
        pdschRxEq = [];
         err_allbws=0;
        for BW = 1:5
            err_allusers=0; 
           grid = rxGrid_7500((1:lenGrid)+(BW-1)*lenGrid,:,:);
%             txgrid = txGrid((1:lenGrid)+(BW-1)*lenGrid,:,:);
%             H = IdChEst_OFDMRx(grid,txgrid);
            % Get LTE PDSCH indices to extract the resource elements
%             pdschIndices = ltePDSCHIndices(gnb, pdsch, pdsch.PRBSet);
            
%             %DMRS  LS channel estimate       
%             hD = ChannelEstimate(grid,txgrid,DMRSIndices);
      
            %CSIRS LS channel estimate
%             hD = ChannelEstimate(grid,txgrid,CSIRSIndices);
            hD_1 = ChannelEstimate(rxGrid_7500((1:lenGrid/2)+(BW-1)*lenGrid,:,:),txGrid_7500((1:lenGrid/2)+(BW-1)*lenGrid,:,:),CSIRSIndices);%1320x12x32x32
            hD_2 = ChannelEstimate(rxGrid_7500((1:lenGrid/2)+(BW-1)*lenGrid+lenGrid/2,:,:),txGrid_7500((1:lenGrid/2)+(BW-1)*lenGrid+lenGrid/2,:,:),CSIRSIndices);%1320x12x32x32
            hD=[hD_1;hD_2];
            
%             EsNO = 10^(SNRdB/10);
%             hD = ChEst_freq_mmseforCSIrs(EsNO, rmcInfo.NDLRB, rmcInfo.Nfft,chinfo.PathSampleDelays,grid,txgrid,CSIRSIndices);
            chEst_1 = ExtChResponse(hD_1, PDSCHIndices)/9; %ori
            chEst_2 = ExtChResponse(hD_2, PDSCHIndices)/9; %ori
            chEst=[chEst_1; chEst_2];
%             chEst = ExtChResponse(hD, PDSCHIndices); 
%             idealchEst = ExtChResponse(H, PDSCHIndices); 
%             idealchEst = idealchEst./(8192/sqrt(6600))*3;
            % Get PDSCH resource elements from the received grid
            for band=1:2
            pdschRx = ExtractResources(PDSCHIndices,grid((band-1)*(lenGrid/2)+(1:lenGrid/2),:,:));
%             pdschTx = ExtractResources(PDSCHIndices,rxGrid_7500);
            pdschRxEq = MIMOReceiver_MMSE(pdschRx, chEst((band-1)*size(chEst_1,1)+(1:size(chEst_1,1)),:,:),nVar);
            
%             for i = 1:32
%                   RMS(i)=sqrt(sum(abs(pdschRxEq(:,i)-pdschTx(:,i)).^2)/length(pdschRxEq));
%             end 
            
            pdschRxEq_all(:,:,BW)=pdschRxEq;%16720x32
              for port=1:32
                  [rxCws,symbols] = ltePDSCHDecode_modify(gnb, pdsch, pdschRxEq_all(:,port,BW));%8360x1
                  codedTrBlock_Eq=rxCws{1,1};
                  codedTrBlock_Eq(codedTrBlock_Eq<0)=0;
                  codedTrBlock_Eq(codedTrBlock_Eq>0)=1;
                  codedTrBlock_Eq(codedTrBlock_Eq == 0) = -1;
                  decState = [];
                  [rxdata,~,decState] = h5gDLSCHDecode_modify(gnb,pdsch,trBlk,codedTrBlock_Eq,decState);%仅仅是一个子带的接收判决后的数据.
%                   txdata=trdataAll(:,port,BW);%拿出所有的子带的发送数据.
                  txdata=trdataAll((band-1)*22152+(1:22152),port,BW);%拿出所有的子带的发送数据.%44304x32
                  err_port=length(find(rxdata~=txdata));
                  rxdataAll(:,port,BW)=rxdata;
                  fprintf('error bits in port %2d 的个数为: %d\n',port,err_port);
                  err_allusers=err_allusers+err_port;%一个子带上的所有用户的错误估计比特.
                  
%                   figure
%                   plot(pdschRxEq(:,port), 'o', 'MarkerEdgeColor', [0 0 0.75], ...
%                       'MarkerFaceColor', [1 0.25 0.25], 'MarkerSize',3);title('均衡后','FontSize',8);
%                   axis([-1,1,-1,1].*2);
          
              end
            end
              err_allbws=err_allbws+err_allusers;
%             figure
%             plot(pdschRx(:,1), 'o', 'MarkerEdgeColor', [0.75 0 0], ...
%                 'MarkerFaceColor', [1 0.25 0.25], 'MarkerSize',3);title('第1层数据均衡前','FontSize',8);
%             axis([-1 1 -1 1]);
%             figure
%             plot(pdschRxEq(:,1), 'o', 'MarkerEdgeColor', [0 0.75 0], ...
%                 'MarkerFaceColor', [1 0.25 0.25], 'MarkerSize',3);title('第1层数据均衡后','FontSize',8);
% %             axis([-1 1 -1 1]);
%             figure
%             plot(pdschRx(:,4), 'o', 'MarkerEdgeColor', [0.75 0 0], ...
%                 'MarkerFaceColor', [1 0.25 0.25], 'MarkerSize',3);title('第4层数据均衡前','FontSize',8);
%             axis([-1 1 -1 1]);
%             figure
%             plot(pdschRxEq(:,4), 'o', 'MarkerEdgeColor', [0 0.75 0], ...
%                 'MarkerFaceColor', [1 0.25 0.25], 'MarkerSize',3);title('第4层数据均衡后','FontSize',8);
%             axis([-1 1 -1 1]);
        end
        fprintf('===============================================================================\n');
        ber_snr(snrIdx,:)=err_allbws./(size(trdataAll,1)*size(trdataAll,2)*5);
        fprintf("SNRIn id is : %2d,对应的误码率为 %d .\n",SNRIn(snrIdx),ber_snr(snrIdx,:));