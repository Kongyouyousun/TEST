%% Simulation Length and SNR Points
% Set the length of the simulation in terms of the number of 10ms frames.
% Set the SNR points to simulate.
clear all;
clc;
NFrames = 1;      % Number of 10ms frames
SNRIn = [100];  % SNR range
tic

% == 加载预信息 ===================================
% load("TEMP\DL_precoding_all_temp_0302_1056.mat");
% load("TEMP\G_all_li_0302_1056.mat");
% load("TEMP\out_0302_1056.mat");
%% gNodeB and PDSCH Configuration
% Set the key parameters of the simulation. These are a combination of NR
% and LTE parameters and specify:
%
% * The bandwidth in resource blocks (12 subcarriers per resource block).
% Note that a larger bandwidth occupancy can be used for W-OFDM and F-OFDM
% waveforms compared to CP-OFDM, for example, 108 RBs vs 100 RBs
% * Subcarrier spacing: 15,30,60,120,240,480 (kHz)
% * Waveform: 'F-OFDM', 'W-OFDM' or 'CP-OFDM'
% * Propagation channel model: 'TDL' or 'CDL'
% * The number of LTE PDSCH transmit antennas (Defined by CellRefP: 1, 2 or 4)
% * LTE transmission scheme: 'Port0', 'TxDiversity' or 'CDD'
% * Target code rate
% * Allocated resource blocks (PRBSet)
% * Modulation scheme: 'QPSK', '16QAM', '64QAM', '256QAM'
% * Transport channel coding: 'LDPC','Turbo'

% Define the propagation channel type and geometry
channelType = 'TDL';     % 'CDL' or 'TDL'
nRxAntennas = 32;         % Number of receive antennas at UE

% Set waveform type and OFDM numerology (SCS and CP type)
simulationParameters = [];                      % clear simulationParameters variable
simulationParameters.WaveformType = 'CP-OFDM';   % 'W-OFDM', 'F-OFDM' or 'CP-OFDM'
simulationParameters.SubcarrierSpacing = 15;    % 15,30,60,120,240,480 (kHz)
simulationParameters.CyclicPrefix = 'Extended';   % 'Normal' or 'Extended'
simulationParameters.UseDCSubcarrier = 'On';    % 'On' or 'Off'
simulationParameters.NTxAnts=1;
% Set number of RB and implicit number of transmit antennas
simulationParameters.NDLRB = 110;   % 20MHz at 15kHz subcarrier spacing
% simulationParameters.CellRefP =32;  % Set to 2 or 4 for 'TxDiversity' or 'CDD' transmission scheme

% DL-SCH/PDSCH parameters
simulationParameters.PDSCH.TargetCodeRate = 0.5;    % Code rate used to calculate transport block sizes
simulationParameters.PDSCH.CodingType = 'LDPC';     % Set to 'LDPC' or 'Turbo'
simulationParameters.PDSCH.PRBSet = (0:simulationParameters.NDLRB-1)';  % Fullband PDSCH allocation
simulationParameters.PDSCH.TxScheme = 'Port0';      % Set to 'Port0', 'TxDiversity' or 'CDD'
simulationParameters.PDSCH.Modulation = {'16QAM'};  % 'QPSK', '16QAM', '64QAM', '256QAM'
simulationParameters.PDSCH.CSI = 'On';              % Turn on LLR scaling after demodulation, set to 'On' or 'Off'
% simulationParameters.PDSCH.NLayers = 4;             % Applicable to CDD (2 or 4, depending on CellRefP)
ncw = 1;                                            % Number of active codewords (can be 2 for CDD)

%%
% Call <matlab:doc('lteRMCDL') lteRMCDL> to generate the default LTE eNodeB
% parameters not specified in |simulationParameters|. These will be
% used later to configure the LTE PDSCH (transport block sizes etc) and
% associated dependencies.

gnb = lteRMCDL(simulationParameters,ncw); % LTE waveform generation parameterssimulationParameters
pdsch = gnb.PDSCH;                        % Separate out the PDSCH parameters
[~,~,rmcInfo] = lteRMCDLTool(gnb,[]);

%%
% The |pdsch| structure contains, amongst other fields, the transport
% block sizes and redundancy version sequence for each transport block/codeword
% within a 10ms frame (10 PDSCH transmissions per 10ms at 15kHz SCS). These
% will be used later in the simulation.

% rvSequence = pdsch.RVSeq;
trBlkSizes = pdsch.TrBlkSizes;

%% Waveform Configuration
% Additional parameters can be set depending on the type of modulation
% waveform used.
%
% For W-OFDM you can specify:
%
% * The window roll-off factor
% * The windowing function samples. This parameter is optional and defaults
% to root raised cosine if unspecified
%
% For F-OFDM you can specify:
%
% * The filter length
% * The tone offset

% Set multi-mode waveform specific parameters
% WOLA specific parameters
gnb.Alpha = 0.0125;
% F-OFDM specific parameters
gnb.FilterLength = 513;
gnb.ToneOffset = 2.5;

%% Propagation Channel Model Configuration
% Create the channel model object. Both CDL and TDL channel models are
% supported.

% Given the use of LTE PDSCH for transmission, use the LTE CellRefP
% parameter to define the number of PDSCH transmit antennas
% nTxAntennas = simulationParameters.CellRefP;
nTxAntennas = 32;   
%%
% The sampling rate for the channel model is set using the value returned
% from <matlab:help('h5gOFDMInfo') h5gOFDMInfo>.

waveformInfo = h5gOFDMInfo(gnb);
% channel.SampleRate = waveformInfo.SamplingRate;
% Initialize variables used in the simulation and analysis
% Array to store the maximum throughput for all SNR points
maxThroughput = zeros(length(SNRIn),1);
% Array to store the simulation throughput for all SNR points
simThroughput = zeros(length(SNRIn),1);
trdataAll=zeros(22152,32,5);%普通cp
UE_NUM=16;

for snrIdx = 1:numel(SNRIn)
    for nn =1:2
    txWave = [];
    CSIRSIndices = [];
    DMRSIndices = [];
    PDSCHIndices = [];
    Peak = zeros(32,1);
    RMS = zeros(32,1);
    pdschSymbolsNum = gnb.NDLRB*12*12-gnb.NDLRB*12*(4+2);
    pdschSymbols = zeros(pdschSymbolsNum,32,5);

    %产生八层PDSCH数据
    for layer = 1:32
        for BW = 1:5
            % Set the random number generator settings to default values
            rng('default');
            
            SNRdB = SNRIn(snrIdx);
%             fprintf('\nSimulating %s (%dx%d) and %s (SCS=%dkHz) with %s channel at %gdB SNR for %d 10ms frame(s)\n',...
%                 pdsch.TxScheme, nTxAntennas, nRxAntennas, ...
%                 gnb.WaveformType, gnb.SubcarrierSpacing,channelType, SNRdB, NFrames);
            
            % Initialize variables used in the simulation and analysis
            blkCRC = [];            % Block CRC for all active PDSCH transmissions
            bitTput = [];           % Number of successfully received bits per transmission
            txedTrBlkSizes = [];    % Number of transmitted info bits per transmission
            
            % Reset the channel so that each SNR point will experience the same channelfunction H = h5gPerfectChannelEstimate(enb,channel,pathGains,toffset)
            % realization
%             reset(channel);
            
            % Total number of OFDM symbols in the simulation period
            NSymbols = NFrames * 10 * waveformInfo.SymbolsPerSubframe;
            
            % OFDM symbol number associated with start of each PDSCH transmission
            gnb.NSymbol = 0;
            % Running counter of the number of PDSCH transmission instances
            npdsch = 0;
            
            %     while  gnb.NSymbol < NSymbols
            
            % Extract the current PDSCH transport block size(s)
            % Use a row vector to help with statistics gathering later
            trBlk = trBlkSizes(:, mod(npdsch, size(trBlkSizes,2))+1).';
             
            % Update effective subframe number to be used to generate the
            % LTE PDSCH indices and scrambling etc.
            % For SCS > 15kHz (15kHz * 2^n) this will compress the LTE frame structure
            % into periods of 1/n ms
            gnb.NSubframe = npdsch;     % The true 1ms subframe number will be fix(gnb.NSymbol/waveformInfo.symbolsPerSubframe)
            
            % LTE PDSCH resource element indices
            [~,pdschInfo] = ltePDSCHIndices(gnb, pdsch, pdsch.PRBSet);
            
            % Perform transport channel coding
            pdschInfo.G = pdschSymbolsNum;
            codedTrBlkSize = pdschInfo.G;   % Available PDSCH bit capacity
            trdata = randi([0,1],trBlk,1);
            trdataAll(:,layer,BW)=trdata;
            codedTrBlock = h5gDLSCH(gnb, pdsch, codedTrBlkSize*4, trdata);
            
            % LTE PDSCH complex symbol generation
            pdschSym = ltePDSCH(gnb, pdsch, codedTrBlock);
            pdschSymbols(:,layer,BW) = pdschSym;
        end
    end
    
    % === VBLAST =================
    for BW = 1:5
    pdschSymbols(:,:,BW) = squeeze(pdschSymbols(:,:,BW))*eye(32);
    end
    
    fprintf('产生用户的原始数据,即没有经过预编码.\n');

    %generate CSI-RS on 32 antenna
    CSIRS = CSIRSgenerator(gnb, 1);
    %generate DMRS on 12 antenna
    DMRS = DMRSgenerator(gnb, 1);
    
    %%映射参考信号
    for port = 1:32
        temp = [];
        pdschIndices = 1:12*gnb.NDLRB*12;  
        for BW = 1:5            
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
            [pdschGrid,csiisindices] = CSIRSmapping_new_15812(pdschGrid,CSIRS(ind1,:),gnb,port);
            
             %DMRS mapping
            if port>=1 && port<=12              
                [pdschGrid,dmrsindices] = DMRSmapping(pdschGrid,DMRS(ind2,:),port);%lo = 3;% DMRS放在OFDM3,4的位置上
            end
            temp = [temp;pdschGrid];
        end
        CSIRSIndices = [CSIRSIndices csiisindices];
        DMRSIndices = [DMRSIndices dmrsindices];
        csiposition=[gnb.NDLRB*12*0+1:gnb.NDLRB*12*1 gnb.NDLRB*12*4+1:gnb.NDLRB*12*5 gnb.NDLRB*12*7+1:gnb.NDLRB*12*8 gnb.NDLRB*12*11+1:gnb.NDLRB*12*12];
%         pdschIndices = ExpungeFrom(pdschIndices,gnb.NDLRB*12*7+1:gnb.NDLRB*12*11);
        pdschIndices = ExpungeFrom(pdschIndices,csiposition.');
        pdschIndices = ExpungeFrom(pdschIndices,gnb.NDLRB*12*2+1:gnb.NDLRB*12*4);
        PDSCHIndices = [PDSCHIndices pdschIndices.'];% 计算每个port中的pdsch的位置（序号）
    end
     fprintf('pdsch,csirs已算出位置 .\n');
     
     % == PDSCH加预编码 ======================================
       %% 数据加预编码.                                 
       if(nn==2)
           DL_precoding_data_all=zeros(size(pdschSymbols,1),32,32,5);% 7920 x 32 x 32 x 5
           for BW=1:5 
               tmp=zeros(size(pdschSymbols,1),size(pdschSymbols,2));% 7920 x 32
               dataElement=size(pdschSymbols,1);% 7920
               [DL_precoding_data,~]=ExtPrecodeResponse(DL_precoding_all_temp(:,:,:,:,BW),pdschIndices);
               DL_precoding_data_all(:,:,:,BW)=DL_precoding_data;
               for dataelement=1:dataElement
                   tmp_PrecodedSignal=zeros(1,32);%1x32
                   for uenum = 1:UE_NUM
                       PrecodedSignal= [];
                       dl_precoding_data=squeeze(DL_precoding_data(dataelement,:,((uenum-1)*2+(1:2))));  %32x2
                       PrecodedSignal=dl_precoding_data*squeeze(pdschSymbols(dataelement,(uenum-1)*2+(1:2))).';%(32x2)x(2x1)=(32x1)
                       tmp_PrecodedSignal=tmp_PrecodedSignal+PrecodedSignal.';%1x32
                   end
                   tmp(dataelement,:)=tmp_PrecodedSignal;% 7920 x 32
               end
               pdschSymbols(:,:,BW)=tmp;
           end
           fprintf('数据已加了预编码 .\n');
       end
       pdschSymbols;% 7920 x 32 x 5
    % ====== endof PDSCH加预编码 ==========================

    %映射PDSCH数据
    txGrid = zeros(gnb.NDLRB*12*5,size(pdschGrid,2),32);
     for port = 1:32
        for BW = 1:5            
            % PDSCH mapping in grid associated with PDSCH transmission period
            pdschGrid = lteDLResourceGrid(gnb);
            pdschGrid(PDSCHIndices(:,port)) = pdschSymbols(:,port,BW);
            
            SingleOFDM=var(pdschGrid(:,1:12));
            
            if BW == 1
                ind1 = 1:length(CSIRS)/5;
                ind2 = 1:length(DMRS)/5;
            else
                ind1 = ind1 + length(CSIRS)/5;
                ind2 = ind2 + length(DMRS)/5;
            end
            %CSI-RS mapping
            pdschGrid = CSIRSmapping_new_15812(pdschGrid,CSIRS(ind1,:),gnb,port);
             
            pdschGrid = DMRSmapping(pdschGrid,DMRS(ind2,:),port);

            txGrid((1:length(pdschGrid))+(BW-1)*length(pdschGrid),:,port) = pdschGrid;
        end
        % OFDM modulation of associated resource elements
        [txWaveform,txinfo] = OFDMModulate(gnb, squeeze(txGrid(:,:,port)));
        txWave = [txWave txWaveform];
     end
        fprintf('OFDM已调制好,正准备发送 .\n');
%       hSpec1 = dsp.SpectrumAnalyzer('SampleRate',  length(txWave)*1e3, ...
%                 'SpectrumType', 'Power density', 'PowerUnits', 'dBW', ...
%                 'RBWSource', 'Property',   'RBW', 15000,...
%                 'FrequencySpan', 'Span and center frequency',...
%                 'Span',  110e6, 'CenterFrequency', 0,...
%                 'SpectralAverages', 10, ...
%                 'Title', '发射信号频谱', 'YLimits', [-200 -0],...
%                 'YLabel', 'PSD');
%             step(hSpec1,txWave(:,1));
    
%             % Write the sampling rate and chip rate to the configuration structure to
%             % allow the calculation of ACLR parameters
%             gnb.SamplingRate=txinfo.SamplingRate;
%             gnb.UTRAChipRate = 3.84;              % UTRA chip rate in MCPS
%             % Calculate ACLR measurement parameters
%             aclr = my5g_hACLRParameters(gnb);
%             % Apply required oversampling
%             resampled = resample(txWaveform(:,1),aclr.OSR,1);
%             % Calculate E-UTRA ACLR
%             aclr = hACLRMeasurementEUTRA(aclr, resampled);
%             % hACLRResults.m displays the ACLR and plots the adjacent channel powers.
%             myhACLRResults(aclr);

    %         % Pass data through channel model. Append zeros at the end of the
    %         % transmitted waveform to flush channel content. These zeroes take
    %         % into account any delay introduced in the channel. This is a mix
    %         % of multipath delay and implementation delay. This value may
    %         % change depending on the sampling rate, delay profile and delay
    %         % spread
    %         txWaveform = [txWaveform; zeros(maxChDelay, size(txWaveform,2))]; %#ok<AGROW>
%             channel.SampleRate = length(txWave)*1e3;
%             [rxWaveform,pathGains] = channel(txWave);
    
    %% MIMO Fading channel with AWGN
    % Initialize channel time for each subframe
    % Channel model
%     channel2.Seed = 10;                   % Random channel seed
%     channel2.NRxAnts = 8;             % receive antennas
%     channel2.DelayProfile = 'EPA';      % Delay profile
%     channel2.DopplerFreq = 20;       % Doppler frequency
%     channel2.MIMOCorrelation = 'Low';   % Multi-antenna correlation
%     channel2.NTerms = 16;                 % Oscillators used in fading model
%     channel2.ModelType = 'GMEDS';         % Rayleigh fading model type
%     channel2.InitPhase = 'Random';        % Random initial phases
%     channel2.NormalizePathGains = 'On';   % Normalize delay profile power
%     channel2.NormalizeTxAnts = 'On';      % Normalize for transmit antennas
%     channel2.InitTime = gnb.NSubframe/1000;
%     channel2.SamplingRate = length(txWave)*1e3;
    %channel.Seed = idx;    % Set seed per subframe for larger variation
    %  - Channel estimation averaging affected
%     rxWaveform = [];
%     for loop = 1:4
%         [wavetmp,chinfo]=lteFadingChannel(channel2, txWave(:,(1:8)+(loop-1)*8));
%         rxWaveform =[rxWaveform  wavetmp];
%     end
     
% % Multi-user system with single/multiple streams per user
% prm.numUsers = 16;                 % Number of users
% prm.numSTSVec = 2*ones(1,16);        % Number of independent data streams per user
% prm.numSTS = sum(prm.numSTSVec);  % Must be a power of 2
% prm.numTx = prm.numSTS;         % Number of BS transmit antennas (power of 2)
% prm.numRx = 2*ones(1,16);      % Number of receive antennas, per user (any >= numSTSVec)
% 
% % Each user has the same modulation
% % prm.bitsPerSubCarrier = 4;   % 2: QPSK, 4: 16QAM, 6: 64QAM, 8: 256QAM
% % prm.numDataSymbols = 10;     % Number of OFDM data symbols
% 
% % MS positions: assumes BS at origin
% %   Angles specified as [azimuth;elevation] degrees
% %   az in range [-180 180], el in range [-90 90], e.g. [45;0]
% maxRange = 3000;            % all MSs within 1000 meters of BS
% prm.mobileRanges = 100*ones(1,16);
% azimuth=-180:22.5:180;
% prm.mobileAngles = [azimuth(1:16);zeros(1,16)];
% 
% prm.fc = 3.5e9;               % 3 GHz system
% prm.chanSRate = 122.88e6;       % Channel sampling rate, 100 Msps
% prm.ChanType = 'Scattering'; % Channel options: 'Scattering', 'MIMO'
% 
% numSTS = prm.numSTS;
% numTx = prm.numTx;
% numRx = prm.numRx;
% 
% prm.cLight = physconst('LightSpeed');
% prm.lambda = prm.cLight/prm.fc;
% 
% % % Get transmit and receive array information
% % [isTxURA,expFactorTx,isRxURA,expFactorRx] = helperArrayInfo(prm);
% 
% % Transmit antenna array definition
% %   Array locations and angles
% prm.posTx = [0;0;0];       % BS/Transmit array position, [x;y;z], meters
% 
%     % Uniform Rectangular array
%     antenna = phased.ShortDipoleAntennaElement(...
%     'FrequencyRange',[3.4e9,3.6e9],...
%     'AxisDirection','Z');
%     txarray = phased.UCA('NumElements',32,'Radius',maxRange,'Element',antenna);
% %     txarray = phased.PartitionedArray(...
% %         'Array',array,...
% %         'SubarraySelection',eye(numTx),'SubarraySteering','Custom');
% pattern(txarray,prm.fc,[-180:180],0,'Type','powerdb',...
%     'CoordinateSystem','polar',...
%     'Normalize',true);
% 
% prm.posTxElem = getElementPosition(txarray)/prm.lambda;
% 
% spLoss = zeros(prm.numUsers,1);
% prm.posRx = zeros(3,prm.numUsers);
% toRxAng = zeros(2,prm.numUsers);
% for uIdx = 1:prm.numUsers    
%     % Uniform Linear array
%     rxarray = phased.ULA(numRx(uIdx), ...
%         'ElementSpacing',0.5*prm.lambda, ...
%         'Element',phased.IsotropicAntennaElement);
%     prm.posRxElem = getElementPosition(rxarray)/prm.lambda;
%     
%     % Mobile positions
%     [xRx,yRx,zRx] = sph2cart(deg2rad(prm.mobileAngles(1,uIdx)), ...
%                              deg2rad(prm.mobileAngles(2,uIdx)), ...
%                              prm.mobileRanges(uIdx));
%     prm.posRx(:,uIdx) = [xRx;yRx;zRx];
%     [toRxRange(uIdx),toRxAng(:,uIdx)] = rangeangle(prm.posTx,prm.posRx(:,uIdx));
%     spLoss(uIdx) = fspl(toRxRange(uIdx),prm.lambda);
% end
% % Account for channel filter delay
% prm.numPadZeros = 0;
% % Apply a spatially defined channel to the transmit signal
% [rxSig,chanDelay] = MyhelperApplyMUChannel(txWave,prm,spLoss);
% for i=1:length(rxSig)
%     rxWaveform(:,(i-1)*2+(1:2)) = rxSig{i};
% end
 
% for uIdx = 1:prm.numUsers  
%     % Front-end amplifier gain and thermal noise
%     rxPreAmp = phased.ReceiverPreamp( ...
%         'Gain',spLoss(uIdx), ...        % account for path loss
%         'SampleRate',prm.chanSRate);
%     rxSigAmp = rxPreAmp(rxSig{uIdx});
% %     Scale power for occupied sub-carriers
%     txWaveform(:,(uIdx-1)*2+(1:2))=rxSigAmp;
% end

%       COST 207
% channeltype = 1;
% Cost=cost_207(channeltype,1/122.88e6,0);
%     channel = nr5gTDLChannel; % TDL channel object
%     % Set the channel geometry
%     channel.DelayProfile = 'Custom';
%     channel.NumTransmitAntennas = nTxAntennas;
%     channel.NumReceiveAntennas = nRxAntennas;
%     channel.TransmissionDirection  = 'Downlink';
%     channel.SampleRate = 122.88e6;
%     channel.PathGainsOutputPort = true;  % Provide the path gains as an output
%     %     channel.MaximumDopplerShift = 70;
%     channel.PathDelays = Cost.PathDelays ;
%     channel.AveragePathGains =  Cost.AvgPathGaindB;
    % Pass data through channel model
%     maxChDelay = ceil(max(channel.PathDelays)*channel.SampleRate);
%     txWave = [txWave; zeros(maxChDelay, size(txWave,2))];
%     [rxWaveform,chPathG]= channel(txWave);
%        in = zeros(size(txWave));
%         if channeltype == 3
%             in(1:540:end,:)=1;
%         elseif channeltype == 2
%             in(1:220:end,:)=1;
%         else
%             in(1:25:end,:)=1;
%         end
         % 修改
%         rxWaveform = MIMOFadingChan(txWave,Cost,32);
%           rxWaveform=txWave;
%         if channeltype == 3
%             out(1:540,:) = hh(1:540,:);
%         elseif channeltype == 2
%             out(1:220,:) = hh(1:220,:);
%         else
%             out(1:25,:) = hh(1:25,:);
%         end
%         for antNo = 1:32
%             rxWaveform(:,antNo)=filter(Cost,txWave(:,antNo));
%             rxWaveform(:,antNo)=circonv(txWave(:,antNo),out(:,antNo),8192);
%             %          rxWaveform(:,antNo)=conv2(txWaveform(:,antNo),out(:,antNo),'same');
%         end

% ========================= 新信道 ========================================
    if (nn==1)
        channeltype=1;
        Channel_model_apply_MUChannel_213_00;% 2019/02/13
        in = zeros(size(txWave));  % 114816x32
        if channeltype == 3
            in(1:dot:end,:)=1;   %530
        elseif channeltype == 2
            in(1:210:end,:)=1;
        else
            in(1:80:end,:)=1;
        end
        %         prm.channeltype=channeltype;
        %         [hh,chanObj] = MIMOFadingChan(in,Cost,UE_NUM); %hh:30720x32， 点对点MIMO模型
        %         [hh_temp,~]=MyhelperApplyMUChannel_4_0(in,prm,spLoss);  %hh_temp:30720x32x16
        [hh_temp,~]=MyhelperApplyMUChannel_4_0_RA(in,prm,spLoss); % 2019/02/13
        %         [hh_temp,~]=MyhelperApplyMUChannel_3_0(in,prm,spLoss);  %hh_temp:30720x32x16
        %         hh=zeros(size(txWaveform));
        %         %16个用户信号叠加作为最终基站接收信号
        %         for i=1:UE_NUM
        %             hh=hh+hh_temp{i};
        %         end
        hh=[];
        UE_NUM=16;
        for i=1:UE_NUM
            hh=[hh hh_temp{i}];
        end
       fprintf('信道数据已产生 .\n');
       
        if channeltype == 3
            out = hh(1:2120,:); % 2114-2120
        elseif channeltype == 2
            out = hh(1:820,:);% 812-820
        else
            out = hh(1:80,:);% 74-80
        end

    end % endof nn==1
    for antNo = 1:32
        tmp=conv(txWave(:,antNo),out(:,antNo));  %
        rxWaveform(:,antNo)=tmp(1:end-79);  % 
    end
    fprintf('PDSCH已经过信道 .\n');
% =========================================================================
    
%     % Calculate linear noise gain
%     SNR = 10^(SNRdB/20);
%     
%     % Normalize noise power to take account of sampling rate, which is
%     % a function of the IFFT size used in OFDM modulation, and the
%     % number of transmit antennas
%     N0 = 1/(sqrt(2.0*size(txWaveform,2)*double(waveformInfo.Nfft*4))*SNR);
%     
%     % Create additive white Gaussian noise
%     noise = N0*complex(randn(size(txWaveform)),randn(size(txWaveform)));
%     
%     % Add AWGN to the received time domain waveform
%     rxWaveform = txWaveform + noise;
    %% Add AWG noise
%     sigPow=10*log10(var(txWave));
%     nVar=10.^(0.1.*(sigPow-SNRdB));
%     rxWaveform=AWGNChannel(rxWaveform,nVar);
    if nn==1
        sigma2=0.0000000000001;  %信道估计用完美估计
    else
%         sigma2=(1/16)/(10^(SNRdB/10));
%         sigma2=0.0000000000001;
%         % 旧的sigma2
        sigPow=10*log10(var(txWave));
        nVar=10.^(0.1.*(sigPow-SNRdB));
        sigma2=nVar(1)/16; % 2019/02/15 16:59 改
    end

    fprintf('信号过了信道并加了噪声 .\n');
    rxGrid = OFDMDemodulate(gnb, rxWaveform);
       
    for i=1:size(rxGrid,2)    %[1 2 3 7 8 9 10 14] %DMRS
        noise_dmrs=sqrt(sigma2/2)* (randn(size(rxGrid(:,i,:)))+1i*randn(size(rxGrid(:,i,:)))); %1320x32
        rxGrid(:,i,:)=rxGrid(:,i,:)+noise_dmrs;
    end
    %         % Perfect synchronization. Use information provided by the channel
    %         % to find the strongest multipath component
    %         [offset,mag] = h5gPerfectTimingOffset(pathGains,chInfo,waveformInfo.SamplingRate);
    %         rxWaveform = rxWaveform(1+offset:end, :);
    %         terror=offset/txinfo.SamplingRate*1e9;
    %         fprintf('定时对齐误差： %0.3fns\n',terror);
    %
%             %Add frequency offset impairment to received waveform
%             foffset = 50.0;         % Frequency offset in Hertz
%             fprintf('add frequency error: %0.1fHz\n',foffset);
%             t = (0:length(rxWaveform)-1).'/txinfo.SamplingRate;
%             rxWaveform = rxWaveform.* repmat(exp(1i*2*pi*foffset*t),1,nRxAntennas);
%             %Synchronization
%             foffset_est = my5g_lteFrequencyOffset(gnb,rxWaveform,pdschGrid);
%     %         if strcmpi(simulationParameters.WaveformType,'W-OFDM')
%     %             foffset_est = -foffset_est;
%     %         end
%             fprintf('detect Frequency offset error: %0.2fHz\n',foffset_est);
%             rxWaveform = rxWaveform.* repmat(exp(-1i*2*pi*foffset_est*t),1,nRxAntennas);
%             foffset_est = my5g_lteFrequencyOffset(gnb,rxWaveform,pdschGrid);
%             fprintf('detect Frequency offset error: %0.2fHz\n',foffset_est);
    %
    
    % Perform OFDM demodulation on the received data to recreate the
    % resource grid
%     for port=1:32
%         rxGrid = OFDMDemodulate(gnb, rxWaveform);
        lenGrid = size(rxGrid,1)/5;
        pdschRx = [];
        pdschTx = [];
        pdschRxEq = [];
         err_allbws=0;
        for BW = 1:5
            fprintf('本次是第%d个子带,刚通过信道,开始解映射 .\n',BW);
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
%             hD = ChannelEstimate_new_15812(grid,txgrid,CSIRSIndices);

          if (nn==1)
            % 信道估计
            hD = ChannelEstimate_new_15812(grid,txgrid,CSIRSIndices);
            fprintf('CSIRS信道估计已完成 .\n');
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
%              if BW==5
%              continue;
%              end
             end % endof nn==1
             
%              if (nn==2)         
            % =============================================================            
%             EsNO = 10^(SNRdB/10);
%             hD = ChEst_freq_mmseforCSIrs(EsNO, rmcInfo.NDLRB, rmcInfo.Nfft,chinfo.PathSampleDelays,grid,txgrid,CSIRSIndices);
%             chEst = ExtChResponse(hD, PDSCHIndices)/9; %ori
%             chEst = ExtChResponse(hD, PDSCHIndices); 
%             idealchEst = ExtChResponse(H, PDSCHIndices); 
%             idealchEst = idealchEst./(8192/sqrt(6600))*3;
            % Get PDSCH resource elements from the received grid
            pdschRx = ExtractResources(PDSCHIndices,grid);
%             pdschTx = ExtractResources(PDSCHIndices,txgrid);
            if nn==1
            chEst = ExtChResponse(hD, PDSCHIndices)/9; %ori
            pdschRxEq = MIMOReceiver_MMSE(pdschRx, chEst,sigma2);
            end
            
            if nn==2
            pdschRxEq =MIMOReceiver_junheng(pdschRx, G_all_li(:,:,:,BW),sigma2);%
            end
            
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
%         end% endof nn==2
      end % endof BW==5 
      if (nn==1)
        fprintf('=========================== NO WIS  ====================================================\n');
        ber_snr1(snrIdx,:)=err_allbws./(size(trdataAll,1)*size(trdataAll,2)*5);
        fprintf("SNRIn id is : %2d,对应的误码率为 %d .\n\n\n",SNRIn(snrIdx),ber_snr1(snrIdx,:));
      else  
        fprintf('==========================  YES WIS ===================================================\n');
        ber_snr2(snrIdx,:)=err_allbws./(size(trdataAll,1)*size(trdataAll,2)*5);
        fprintf("SNRIn id is : %2d,对应的误码率为 %d .\n\n\n",SNRIn(snrIdx),ber_snr2(snrIdx,:));
      end
    end  % endof for nn==1:2         
end   
toc
