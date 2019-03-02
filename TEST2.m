            chEst = ExtChResponse(hD, PDSCHIndices)/9; %ori
%             chEst = ExtChResponse(hD, PDSCHIndices); 
%             idealchEst = ExtChResponse(H, PDSCHIndices); 
%             idealchEst = idealchEst./(8192/sqrt(6600))*3;
            % Get PDSCH resource elements from the received grid
            pdschRx = ExtractResources(PDSCHIndices,grid);
            pdschTx = ExtractResources(PDSCHIndices,txgrid);
%             pdschRxEq = MIMOReceiver_MMSE(pdschRx, chEst,nVar);
            pdschRxEq2 =MIMOReceiver_junheng(pdschRx, G_all_li(:,:,:,BW),sigma2);%
            
%             for i = 1:32
%                   RMS(i)=sqrt(sum(abs(pdschRxEq(:,i)-pdschTx(:,i)).^2)/length(pdschRxEq));
%             end 
            
            pdschRxEq_all(:,:,BW)=pdschRxEq2;%16720x32
              for port=1:32
                  [rxCws,symbols] = ltePDSCHDecode_modify(gnb, pdsch, pdschRxEq_all(:,port,BW));%8360x1
                  codedTrBlock_Eq=rxCws{1,1};
                  codedTrBlock_Eq(codedTrBlock_Eq<0)=0;
                  codedTrBlock_Eq(codedTrBlock_Eq>0)=1;
                  codedTrBlock_Eq(codedTrBlock_Eq == 0) = -1;
                  decState = [];
                  [rxdata,~,decState] = h5gDLSCHDecode_modify(gnb,pdsch,trBlk,codedTrBlock_Eq,decState);%������һ���Ӵ��Ľ����о��������.
                  txdata=trdataAll(:,port,BW);%�ó����е��Ӵ��ķ�������.
                  err_port=length(find(rxdata~=txdata));
                  rxdataAll(:,port,BW)=rxdata;
                  fprintf('error bits in port %2d �ĸ���Ϊ: %d\n',port,err_port);
                  err_allusers=err_allusers+err_port;%һ���Ӵ��ϵ������û��Ĵ�����Ʊ���.
                  
%                   figure
%                   plot(pdschRxEq(:,port), 'o', 'MarkerEdgeColor', [0 0 0.75], ...
%                       'MarkerFaceColor', [1 0.25 0.25], 'MarkerSize',3);title('�����','FontSize',8);
%                   axis([-1,1,-1,1].*2);
          
              end
              err_allbws=err_allbws+err_allusers;
%             figure
%             plot(pdschRx(:,1), 'o', 'MarkerEdgeColor', [0.75 0 0], ...
%                 'MarkerFaceColor', [1 0.25 0.25], 'MarkerSize',3);title('��1�����ݾ���ǰ','FontSize',8);
%             axis([-1 1 -1 1]);
%             figure
%             plot(pdschRxEq(:,1), 'o', 'MarkerEdgeColor', [0 0.75 0], ...
%                 'MarkerFaceColor', [1 0.25 0.25], 'MarkerSize',3);title('��1�����ݾ����','FontSize',8);
% %             axis([-1 1 -1 1]);
%             figure
%             plot(pdschRx(:,4), 'o', 'MarkerEdgeColor', [0.75 0 0], ...
%                 'MarkerFaceColor', [1 0.25 0.25], 'MarkerSize',3);title('��4�����ݾ���ǰ','FontSize',8);
%             axis([-1 1 -1 1]);
%             figure
%             plot(pdschRxEq(:,4), 'o', 'MarkerEdgeColor', [0 0.75 0], ...
%                 'MarkerFaceColor', [1 0.25 0.25], 'MarkerSize',3);title('��4�����ݾ����','FontSize',8);
%             axis([-1 1 -1 1]);
%         end% endof nn==2