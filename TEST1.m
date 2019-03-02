  for BW = 1:5
            err_allusers=0; 
            grid = rxGrid((1:lenGrid)+(BW-1)*lenGrid,:,:);
            txgrid = txGrid((1:lenGrid)+(BW-1)*lenGrid,:,:);
%             H = IdChEst_OFDMRx(grid,txgrid);
            % Get LTE PDSCH indices to extract the resource elements
%             pdschIndices = ltePDSCHIndices(gnb, pdsch, pdsch.PRBSet);
            
%             %DMRS  LS channel estimate       
%             hD = ChannelEstimate(grid,txgrid,DMRSIndices);
      
            %CSIRS LS channel estimate
%             hD = ChannelEstimate(grid,txgrid,CSIRSIndices);
            hD = ChannelEstimate_new_15812(grid,txgrid,CSIRSIndices);
 %           
            % ======================== 计算预编码 ====================================
            hD_uplink=permute(hD,[1 2 4 3]);  %转置.Rx 与 Tx互换位置.
            %             hD_uplink=hD;
            [num_n,num_slot,~,~]=size(hD_uplink); %1320 12
            % 		UL_precoding = zeros(2, 2, UE_NUM);			% Precoding matrix for DL/UL transmission of each user. 2x2x16
            DL_precoding = zeros(32, 2, UE_NUM);
            %         UL_precoding_all_temp=zeros(num_n,num_slot,2,32); %1320x14x2x32
            G_hat=zeros(32,32,num_n);  %12x12x1320信道估计矩阵
            W_IS_all = zeros(num_n,num_slot,32, 32);			% Interference suppresion/BD Multi-user precoding,  Unitary matrix. 32x32
            %每个资源块只做1次预编码和干扰抑制
            for i=1:12:num_n %110RB， num_n:1320
                H=squeeze(hD_uplink(i,1,:,:)); %32x32
                % ZF precoding matrix
                G_hat = conj(H)*inv(H.'*conj(H));  %32x32
                %svd分解出W_IS矩阵
                for k = 1:UE_NUM
                    
                    tmp_idx = (k-1)*2 + (1:2);					% for all of the antennas of user k
                    G_k = G_hat(:, tmp_idx); %32x2
                    
                    % Compute the multi-user precoding matrix, Q_k.
                    [Q_k, R_hat_k] = qr(G_k);
                    
                    % Compute the single-user precoding matrix, V_k for DL and U_k for UL.
                    R_k = inv(R_hat_k(1:2, :)); %2x2
                    [U_k, D, V_k] = svd(R_k);
                    
                    % ### Downlink precoding matrix
                    % Compute the downlink compound precoding matrix of user k. W_DL_k is also the precoded RS matrix of user k.
                    % Different user's pilot matrix (W_DL_k) is transmitted in different OFDM symbols.
                    W_DL_k = Q_k(:,1:2)*V_k;						% Be carefull!!!! You should use Q_k(:,1:Nt)*V_k(:,1:RI(k, n_BPU)). Here I want to make all the matrices have the same size for all N_BPU and all users.
                    DL_precoding(:, :, k) = W_DL_k;			% Keep in mind, only 1:RI(k, n_BPU) columns are usefull. 12x4
                    
                    for t=1:num_slot %12
                        %                         W_IS_all(i,t,:, (k-1)*2 + (1:2)) = Q_k(:,1:2);				% Uplink interference suppression matrix of user k. 1320x12x32x32
                        DL_precoding_all_temp(i,t,:,(k-1)*2 + (1:2),BW)= W_DL_k./sqrt(2);  %把所有用户预编码拼成2行 , 1320x12x2x32
                    end
                end  %end for k=1:uenum
            end % end for 1:12:num_n
            for ii=1:110  %1:6
                DL_precoding_all_temp( (ii-1)*12+(1:12),:,:,:,BW)=repmat(DL_precoding_all_temp((ii-1)*12+1,:,:,:,BW),12,1,1,1);
            end
            fprintf('预编码的计算已完成 .\n');
%}            
             % ################### hD*WV=RV 理论值.理论值可以保存吗?可以,但数值必须是在高信噪比下计算得出. ###############################
            % ################## 放在预编码过程的下面,在nn==1时可以算出.#######################################################
            hDt=permute(hD,[3 4 1 2]);
            WVtmp=permute(DL_precoding_all_temp,[3 4 1 2 5]);%32x32x1320x14
            WV=WVtmp(:,:,:,:,BW);
            G_k=complex(zeros(2,32));
            G=complex(zeros(32,32));
            G_all=complex(zeros(32,32,size(grid,1),size(grid,2)));
            for sc=1:size(grid,1)
                for sys=1:size(grid,2)
                    for k=1:16 %uenum
                        hD_k=hDt((k-1)*2+(1:2),:,sc,sys);%2x32
                        for j=1:16
                            WV_j=WV(:,(j-1)*2+(1:2),sc,sys);%32x2
                            G_kj=hD_k*WV_j;%(2x32)*(32x2)=(2x2)
                            G_k(:,(j-1)*2+(1:2))=G_kj;
                        end
                        G((k-1)*2+(1:2),:)=G_k;
                    end
                    G_all(:,:,sc,sys)=G;%32x32x1320x12
                end
            end %over
            
             G_all_li(:,:,:,BW)=ExtGallResponse(G_all,pdschIndices);%32X32X7920X5
             fprintf('RV 计算已完成 .\n');
            
            
            
            % =============================================================            
%             EsNO = 10^(SNRdB/10);
%             hD = ChEst_freq_mmseforCSIrs(EsNO, rmcInfo.NDLRB, rmcInfo.Nfft,chinfo.PathSampleDelays,grid,txgrid,CSIRSIndices);
            chEst = ExtChResponse(hD, PDSCHIndices)/9; %ori
%             chEst = ExtChResponse(hD, PDSCHIndices); 
%             idealchEst = ExtChResponse(H, PDSCHIndices); 
%             idealchEst = idealchEst./(8192/sqrt(6600))*3;
            % Get PDSCH resource elements from the received grid
            pdschRx = ExtractResources(PDSCHIndices,grid);
            pdschTx = ExtractResources(PDSCHIndices,txgrid);
            pdschRxEq = MIMOReceiver_MMSE(pdschRx, chEst,nVar);
            
            for i = 1:32
                  RMS(i)=sqrt(sum(abs(pdschRxEq(:,i)-pdschTx(:,i)).^2)/length(pdschRxEq));
            end 
            
            pdschRxEq_all(:,:,BW)=pdschRxEq;%16720x32
              for port=1:32
                  [rxCws,symbols] = ltePDSCHDecode_modify(gnb, pdsch, pdschRxEq_all(:,port,BW));%8360x1
                  codedTrBlock_Eq=rxCws{1,1};
                  codedTrBlock_Eq(codedTrBlock_Eq<0)=0;
                  codedTrBlock_Eq(codedTrBlock_Eq>0)=1;
                  codedTrBlock_Eq(codedTrBlock_Eq == 0) = -1;
                  decState = [];
                  [rxdata,~,decState] = h5gDLSCHDecode_modify(gnb,pdsch,trBlk,codedTrBlock_Eq,decState);%仅仅是一个子带的接收判决后的数据.
                  txdata=trdataAll(:,port,BW);%拿出所有的子带的发送数据.
                  err_port=length(find(rxdata~=txdata));
                  rxdataAll(:,port,BW)=rxdata;
                  fprintf('error bits in port %2d 的个数为: %d\n',port,err_port);
                  err_allusers=err_allusers+err_port;%一个子带上的所有用户的错误估计比特.
                  
%                   figure
%                   plot(pdschRxEq(:,port), 'o', 'MarkerEdgeColor', [0 0 0.75], ...
%                       'MarkerFaceColor', [1 0.25 0.25], 'MarkerSize',3);title('均衡后','FontSize',8);
%                   axis([-1,1,-1,1].*2);
          
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