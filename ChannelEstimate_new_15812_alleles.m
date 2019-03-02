function chEst = ChannelEstimate_new_15812_alleles(grid,txgrid,csiind,dmrsind)
chEst = complex(zeros([size(grid) 32]));   % Iniitalize Output
%pdsch 信道估计 
% for numRx = 1:32
%     rxGrid = grid(:,:,numRx);
%     for numTx =1:32
%         txGrid = txgrid(:,:,numTx);
%         hp = rxGrid./txGrid;
%         chEst(:,:,numRx,numTx)=hp;
%     end
% end

% chEst = complex(zeros([size(grid) 32]));
% % DMRS 信道估计, 经过测试,DMRS估计是不准的。
% for numRx = 1:32
%     rxGrid = grid(:,:,numRx);
%     for numTx =1:12
%         txGrid = txgrid(:,:,numTx);
%         Tx = lteExtractResources(dmrsind(:,numTx),txGrid);
%         Rx = lteExtractResources(dmrsind(:,numTx),rxGrid);
%         hp = Rx./Tx;
%         hp = reshape(hp,length(hp)/2,2);
%         P_EST = complex(zeros(size(grid,1),1));
%         %%获取每个天线对应的资源网格的DMRS的位置
%         indices = dmrsind(1:length(dmrsind)/2,numTx)-size(grid,1)*2;
%         h1 = mean([hp(:, 1), hp(:, 2)],2);
%         hp =  h1(1:2:end) +  h1(2:2:end);
%         P_EST(indices(1:2:end)) = hp;
%         
%         y =  interp1(find(P_EST~=0),P_EST(P_EST~=0),(1:length(P_EST)).','linear','extrap');
%         y34 = y(:,ones(1,size(txGrid,2)));
%         chEst(:,:,numRx,numTx)=y34;
%     end
% end

chEst = complex(zeros([size(grid) 32]));
% CSIRS 信道估计 ，测试的误码率是准的
for numRx = 1:32
    rxGrid = grid(:,:,numRx);
    for numTx =1:32
        txGrid = txgrid(:,:,numTx);
        Tx = lteExtractResources(csiind(:,numTx),txGrid);
        Rx = lteExtractResources(csiind(:,numTx),rxGrid);
        Tx = reshape(Tx,length(Tx)/4,4);
        Rx = reshape(Rx,length(Rx)/4,4);
        hp = Rx./Tx;
        
        %         P_EST = complex(zeros(size(grid,1),1));
        
        %%获取每个天线对应的资源网格的CSIRS的位置
%         indices = ind(1:length(ind)/4,numTx)-size(grid,1)*7;
        indices = csiind(1:length(csiind)/4,numTx)-size(grid,1)*0;
%         h1_a = mean([hp(:, 1), hp(:, 2)],2);
%         h1_b = mean([hp(:, 3), hp(:, 4)],2);
%         h1 = mean([h1_a, h1_b],2);
%         hp =  h1(1:2:end) +  h1(2:2:end);
%         P_EST(indices(1:2:end)) = hp/2;
%         y =  interp1(find(P_EST~=0),P_EST(P_EST~=0),(1:length(P_EST)).','linear','extrap');
  
        matrix_Hi_r0 = hp(1:2:end,1)/2+hp(2:2:end,1)/2;       
        xIndex = 1:numel(matrix_Hi_r0);
        zIndex = 1:(1/12):numel(matrix_Hi_r0);
        z= interp1(xIndex',matrix_Hi_r0,zIndex','linear');
        Edges = [indices(1)-1 size(grid,1)-indices(end)+1];
%         Edges = [indices(1)-1 size(grid,1)-indices(end)];
        delta          = z(2)-z(1);
        z_before    = z(1)     -  delta*(Edges(1):-1:1)';
        z_after       = z(end) + delta*(1:Edges(2))';
        y1             = [z_before;z;z_after];
        
        matrix_Hi_r1 = hp(1:2:end,2)/2+hp(2:2:end,2)/2;       
        z= interp1(xIndex',matrix_Hi_r1,zIndex','linear');
        delta          = z(2)-z(1);
        z_before    = z(1)     -  delta*(Edges(1):-1:1)';
        z_after       = z(end) + delta*(1:Edges(2))';
        y5             = [z_before;z;z_after];
        
         matrix_Hi_r2 = hp(1:2:end,3)/2+hp(2:2:end,3)/2;       
        z= interp1(xIndex',matrix_Hi_r2,zIndex','linear');
        delta          = z(2)-z(1);
        z_before    = z(1)     -  delta*(Edges(1):-1:1)';
        z_after       = z(end) + delta*(1:Edges(2))';
        y8             = [z_before;z;z_after];
        
         matrix_Hi_r3 = hp(1:2:end,4)/2+hp(2:2:end,4)/2;       
        z= interp1(xIndex',matrix_Hi_r3,zIndex','linear');
        delta          = z(2)-z(1);
        z_before    = z(1)     -  delta*(Edges(1):-1:1)';
        z_after       = z(end) + delta*(1:Edges(2))';
        y12             = [z_before;z;z_after];
        
        y2=(3/4)*y1+(1/4)*y5;
        y3=(2/4)*y1+(2/4)*y5;
        y4=(1/4)*y1+(3/4)*y5;
        
        y6=(2/3)*y5+(1/3)*y8;
        y7=(1/3)*y5+(2/3)*y8;
        
        y9=(3/4)*y8+(1/4)*y12;
        y10=(2/4)*y8+(2/4)*y12;
        y11=(1/4)*y8+(3/4)*y12;
        
%         y15812 = [y1 y5 y8 y12];
        y=[y1 y2 y3 y4 y5 y6 y7 y8 y9 y10 y11 y12];
%         y=mean([y1  y4],2);%改 0301 22：24 
%         y=y1;
%         chEst(:,[1 5 8 12],numRx,numTx)=y15812;
       chEst(:,:,numRx,numTx)=y;
    end
end
% 
% 
