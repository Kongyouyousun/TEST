%% Simulation Length and SNR Points
% Set the length of the simulation in terms of the number of 10ms frames.
% Set the SNR points to simulate.
clear all ,clc
NFrames = 1;      % Number of 10ms frames
SNRIn = 60;  % SNR range

channelType = 'TDL';     % 'CDL' or 'TDL'
nRxAntennas = 32;         % Number of receive antennas at UE

% Set waveform type and OFDM numerology (SCS and CP type)
simulationParameters = [];                      % clear simulationParameters variable
simulationParameters.WaveformType = 'CP-OFDM';   % 'W-OFDM', 'F-OFDM' or 'CP-OFDM'
simulationParameters.SubcarrierSpacing = 15;    % 15,30,60,120,240,480 (kHz)
simulationParameters.CyclicPrefix = 'Normal';   % 'Normal' or 'Extended'
simulationParameters.UseDCSubcarrier = 'On';    % 'On' or 'Off'
simulationParameters.NTxAnts=1;
% Set number of RB and implicit number of transmit antennas
simulationParameters.NDLRB = 110;   % 20MHz at 15kHz subcarrier spacing
% simulationParameters.CellRefP =32;  % Set to 2 or 4 for 'TxDiversity' or 'CDD' transmission scheme

% DL-SCH/PDSCH parameters
simulationParameters.PDSCH.TargetCodeRate = 0.5;    % Code rate used to calculate transport block sizes
simulationParameters.PDSCH.CodingType = 'Turbo';     % Set to 'LDPC' or 'Turbo' 
simulationParameters.PDSCH.PRBSet = (0:simulationParameters.NDLRB-1)';  % Fullband PDSCH allocation
simulationParameters.PDSCH.TxScheme = 'Port0';      % Set to 'Port0', 'TxDiversity' or 'CDD'
simulationParameters.PDSCH.Modulation = {'16QAM'};  % 'QPSK', '16QAM', '64QAM', '256QAM'
simulationParameters.PDSCH.CSI = 'On';              % Turn on LLR scaling after demodulation, set to 'On' or 'Off'
% simulationParameters.PDSCH.NLayers = 4;             % Applicable to CDD (2 or 4, depending on CellRefP)
ncw = 1;                                            % Number of active codewords (can be 2 for CDD)

gnb = lteRMCDL(simulationParameters,ncw); % LTE waveform generation parameterssimulationParameters
pdsch = gnb.PDSCH;                        % Separate out the PDSCH parameters
[~,~,rmcInfo] = lteRMCDLTool(gnb,[]);

trBlkSizes = pdsch.TrBlkSizes;

gnb.Alpha = 0.0125;
% F-OFDM specific parameters
gnb.FilterLength = 513;
gnb.ToneOffset = 2.5;

nTxAntennas = 32;

if strcmpi(channelType,'CDL')
    channel = nr5gCDLChannel; % CDL channel object
    
    % Use CDL-A model (Indoor hotspot model (R1-161736, R1-162960))
    channel.DelayProfile = 'CDL-A';

    antarrays = ...
        [1   1   1   1   1;   % 1 ants
        1   1   2   1   1;   % 2 ants
        2   1   2   1   1;   % 4 ants
        2   2   2   1   1;   % 8 ants
        2   4   2   1   1;   % 16 ants
        4   4   2   1   1;   % 32 ants
        4   4   2   1   2;   % 64 ants
        4   8   2   1   2;   % 128 ants
        4   8   2   2   2;   % 256 ants
        8   8   2   2   2;   % 512 ants
        8  16   2   2   2];  % 1024 ants
    antselected = 1+fix(log2(nTxAntennas));
    channel.TransmitAntennaArray.Size = antarrays(antselected,:);

    if nRxAntennas == 1

        channel.ReceiveAntennaArray.Size = ones(1,5);
    else

        channel.ReceiveAntennaArray.Size = [fix(nRxAntennas/2),1,2,1,1];
    end
    
elseif strcmpi(channelType,'TDL')
    channel = nr5gTDLChannel; % TDL channel object
    % Set the channel geometry
    channel.DelayProfile = 'TDL-A';
    channel.NumTransmitAntennas = nTxAntennas;
    channel.NumReceiveAntennas = nRxAntennas;
    channel.PathGainsOutputPort = true;  % Provide the path gains as an output
    %     channel.MaximumDopplerShift = 70;
    if strcmpi(simulationParameters.WaveformType,'W-OFDM')
        channel.DelaySpread = 30e-9;
    else
        channel.DelaySpread = 300e-9;
    end
else
    error('channelType parameter must be either CDL or TDL');
end


waveformInfo = h5gOFDMInfo(gnb);

chInfo = info(channel);
maxChDelay = ceil(max(chInfo.PathDelays*channel.SampleRate)) + chInfo.ChannelFilterDelay;

maxThroughput = zeros(length(SNRIn),1);

simThroughput = zeros(length(SNRIn),1);
err_AllBWs=0;
for snrIdx = 1:numel(SNRIn)
    txWave = [];
    CSIRSIndices = [];
    DMRSIndices = [];
    Peak = zeros(32,1);
    RMS = zeros(32,1);
    pdschSymbols = zeros(8360,32,5);
    trdataAll=zeros(25456,32,5);
    %产生八层PDSCH数据
    for layer = 1:32
        for BW = 1:5
            % Set the random number generator settings to default values
            rng('default');
            
            SNRdB = SNRIn(snrIdx);
            fprintf('\nSimulating %s (%dx%d) and %s (SCS=%dkHz) with %s channel at %gdB SNR for %d 10ms frame(s)\n',...
                pdsch.TxScheme, nTxAntennas, nRxAntennas, ...
                gnb.WaveformType, gnb.SubcarrierSpacing,channelType, SNRdB, NFrames);
            
            % Initialize variables used in the simulation and analysis
            blkCRC = [];            % Block CRC for all active PDSCH transmissions
            bitTput = [];           % Number of successfully received bits per transmission
            txedTrBlkSizes = [];    % Number of transmitted info bits per transmission

            NSymbols = NFrames * 10 * waveformInfo.SymbolsPerSubframe;

            gnb.NSymbol = 0;

            npdsch = 0;

            trBlk = trBlkSizes(:, mod(npdsch, size(trBlkSizes,2))+1).';
            

            gnb.NSubframe = npdsch;     % The true 1ms subframe number will be fix(gnb.NSymbol/waveformInfo.symbolsPerSubframe)
            

            [~,pdschInfo] = ltePDSCHIndices(gnb, pdsch, pdsch.PRBSet);
            
            % Perform transport channel coding
            pdschInfo.G = 8360;
            codedTrBlkSize = pdschInfo.G;   % Available PDSCH bit capacity
            trdata = randi([0,1],trBlk,1);%25456
            trdataAll(:,layer,BW)=trdata;
            codedTrBlock = h5gDLSCH(gnb, pdsch, codedTrBlkSize*4, trdata);
            
            % LTE PDSCH complex symbol generation
            pdschSym = ltePDSCH(gnb, pdsch, codedTrBlock);
            pdschSymbols(:,layer,BW) = pdschSym;
        end
    end
     use_data_Allbits=size(trdataAll,1)*size(trdataAll,2);%8360x32=814592
    for BW = 1:5
    pdschSymbols(:,:,BW) = squeeze(pdschSymbols(:,:,BW))*eye(32);
    end
    
    pdschIndices =  gnb.NDLRB*12*3+1:14*gnb.NDLRB*12;    
    %generate CSI-RS on 32 antenna
    CSIRS = CSIRSgenerator(gnb, 1);
    %generate DMRS on 12 antenna
    DMRS = DMRSgenerator(gnb, 1);
    
    %%映射参考信号
    for port = 1:32
        temp = [];
        for BW = 1:1            
            % PDSCH mapping in grid associated with PDSCH transmission period
            pdschGrid = lteDLResourceGrid(gnb);

            if BW == 1
                ind1 = 1:length(CSIRS)/5;
                ind2 = 1:length(DMRS)/5;
            else
                ind1 = ind1 + length(CSIRS)/5;
                ind2 = ind2 + length(DMRS)/5;
            end
            %CSI-RS mapping
            [pdschGrid,csiisindices] = CSIRSmapping(pdschGrid,CSIRS(ind1,:),gnb,port);
            
             %DMRS mapping
            if port>=1 && port<=12              
                [pdschGrid,dmrsindices] = DMRSmapping(pdschGrid,DMRS(ind2,:),port);
            end
            temp = [temp;pdschGrid];
        end
        CSIRSIndices = [CSIRSIndices csiisindices];
        if port>=1 && port<=12
            DMRSIndices = [DMRSIndices dmrsindices];
        end
    end
    pdschIndices = ExpungeFrom(pdschIndices,CSIRSIndices(:,1:4));
    pdschIndices = ExpungeFrom(pdschIndices,DMRSIndices(:,1:6));
    %映射PDSCH数据
    txGrid = zeros(gnb.NDLRB*12*5,14,32);
     for port = 1:32
        for BW = 1:1            
            % PDSCH mapping in grid associated with PDSCH transmission period
            pdschGrid = lteDLResourceGrid(gnb);
            pdschGrid(pdschIndices) = pdschSymbols(:,port,BW);
            
            
            if BW == 1
                ind1 = 1:length(CSIRS)/5;
                ind2 = 1:length(DMRS)/5;
            else
                ind1 = ind1 + length(CSIRS)/5;
                ind2 = ind2 + length(DMRS)/5;
            end
            %CSI-RS mapping
            pdschGrid = CSIRSmapping(pdschGrid,CSIRS(ind1,:),gnb,port);
              
            if port>=1 && port<=12               
                %DMRS mapping
                pdschGrid = DMRSmapping(pdschGrid,DMRS(ind2,:),port);
            end
            txGrid((1:length(pdschGrid))+(BW-1)*length(pdschGrid),:,port) = pdschGrid;
        end
        % OFDM modulation of associated resource elements
        [txWaveform,txinfo] = OFDMModulate(gnb, squeeze(txGrid(:,:,port)));
        
        txWave = [txWave txWaveform];
     end

    %% MIMO Fading channel with AWGN
    % Initialize channel time for each subframe
    % Channel model
    channel2.Seed = 10;                   % Random channel seed
    channel2.NRxAnts = 8;             % receive antennas
    channel2.DelayProfile = 'EPA';      % Delay profile  Extended Pedestrian A model (EPA) Extended Vehicular A model (EVA)  Extended Typical Urban model (ETU)
    channel2.DopplerFreq = 20;       % Doppler frequency
    channel2.MIMOCorrelation = 'Low';   % Multi-antenna correlation
    channel2.NTerms = 16;                 % Oscillators used in fading model
    channel2.ModelType = 'GMEDS';         % Rayleigh fading model type
    channel2.InitPhase = 'Random';        % Random initial phases
    channel2.NormalizePathGains = 'On';   % Normalize delay profile power
    channel2.NormalizeTxAnts = 'On';      % Normalize for transmit antennas
    channel2.InitTime = gnb.NSubframe/1000;
    channel2.SamplingRate = length(txWave)*1e3;

    rxWaveform = [];
    for loop = 1:4
    rxWaveform =[rxWaveform  lteFadingChannel(channel2, txWave(:,(1:8)+(loop-1)*8))];    %set channel parameters and signal go through the channel, the derived signal is rxWaveform
    end

    %% Add AWG noise
    sigPow=10*log10(var(txWave));
    nVar=10.^(0.1.*(sigPow-SNRdB));
    rxWaveform=AWGNChannel(rxWaveform,nVar);                                         %rxWaveform went through the channel and add AWG noise
   
        rxGrid = OFDMDemodulate(gnb, rxWaveform);
        lenGrid = size(rxGrid,1)/5;
        pdschRx = [];
        pdschTx = [];
        pdschRxEq = [];
        ber_sim_time=0;
          err_ALLusers=0;
        for BW = 1:1
            grid = rxGrid((1:lenGrid)+(BW-1)*lenGrid,:,:);
            txgrid = txGrid((1:lenGrid)+(BW-1)*lenGrid,:,:);
            % Get LTE PDSCH indices to extract the resource elements
%             pdschIndices = ltePDSCHIndices(gnb, pdsch, pdsch.PRBSet);
            
%             %DMRS  LS channel estimate       
%             hD = ChannelEstimate(grid,txgrid,DMRSIndices);
      
            %CSIRS LS channel estimate
            hD = ChannelEstimate(grid,txgrid,CSIRSIndices);
            chEst = ExtChResponse(hD, pdschIndices); 
                
            % Get PDSCH resource elements from the received grid
            pdschRx = lteExtractResources(pdschIndices,grid);
            pdschTx = lteExtractResources(pdschIndices,txgrid);
            pdschRxEq = 2*MIMOReceiver_MMSE(pdschRx, chEst, nVar);

         for port=1:32
           
           [rxCws,symbols] = ltePDSCHDecode(gnb, pdsch, pdschRxEq(:,port,BW));%8360x1
           codedTrBlock_Eq=rxCws{1,1};
           codedTrBlock_Eq(codedTrBlock_Eq<0)=0;
           codedTrBlock_Eq(codedTrBlock_Eq>0)=1;
           codedTrBlock_Eq(codedTrBlock_Eq == 0) = -1;
           decState = [];
           [rxdata,~,decState] = h5gDLSCHDecode_modify(gnb,pdsch,trBlk,codedTrBlock_Eq,decState);
           txdata=trdataAll(:,port,BW);
           err_port=length(find(rxdata~=txdata));
           fprintf('error in port %d 的个数为: %d\n',port,err_port);
           err_ALLusers=err_ALLusers+err_port;%一个子带上的所有用户的错误估计比特.
           
       end

       fprintf('BW id : %d, errors 的总个数为: %d \n',BW,err_ALLusers);
       err_AllBWs=err_AllBWs+err_ALLusers;
       fprintf('===============================================================================\n');
        end   
       

     ber_sim_time=ber_sim_time+err_AllBWs;

    ber_snrIdx(snrIdx,:)=ber_sim_time./(use_data_Allbits);
    fprintf("SNRIn id is : %d,对应的误码率为 %d .\n",SNRIn(snrIdx),ber_snrIdx(snrIdx,:));
end
