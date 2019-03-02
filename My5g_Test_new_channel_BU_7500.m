%% Simulation Length and SNR Points
% Set the length of the simulation in terms of the number of 10ms frames.
% Set the SNR points to simulate.
clear all;
clc;
NFrames = 1;      % Number of 10ms frames
SNRIn = [15:30];  % SNR range
tic
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
% trdataAll=zeros(22152,32,5);%普通cp

for snrIdx = 1:numel(SNRIn)
    txWave = [];
    CSIRSIndices = [];
    DMRSIndices = [];
    PDSCHIndices = [];
    Peak = zeros(32,1);
    RMS = zeros(32,1);
    pdschSymbolsNum = gnb.NDLRB*12*12-gnb.NDLRB*72;
    pdschSymbols = zeros(pdschSymbolsNum,32,5);
     txGrid_7500=[];
     trdataAll=zeros(2*22152,32,5);% 扩展cp
    for band=1:2
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
            trdataSband(:,layer,BW)=trdata;%25456  先保留,是为了在接收端比较错误的误比特.

%             trdataAll(:,layer,BW)=trdata;
            codedTrBlock = h5gDLSCH(gnb, pdsch, codedTrBlkSize*4, trdata);
            
            % LTE PDSCH complex symbol generation
            pdschSym = ltePDSCH(gnb, pdsch, codedTrBlock);
            pdschSymbols(:,layer,BW) = pdschSym;
        end
    end
    
    for BW = 1:5
    pdschSymbols(:,:,BW) = squeeze(pdschSymbols(:,:,BW))*eye(32);
    end
    trdataAll((band-1)*22152+(1:22152),:,:)=trdataSband;
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
            [pdschGrid,csiisindices] = CSIRSmapping(pdschGrid,CSIRS(ind1,:),gnb,port);
            
             %DMRS mapping
            if port>=1 && port<=12              
                [pdschGrid,dmrsindices] = DMRSmapping(pdschGrid,DMRS(ind2,:),port);
            end
            temp = [temp;pdschGrid];
        end
        CSIRSIndices = [CSIRSIndices csiisindices];
        DMRSIndices = [DMRSIndices dmrsindices];
        pdschIndices = ExpungeFrom(pdschIndices,gnb.NDLRB*12*7+1:gnb.NDLRB*12*11);
        pdschIndices = ExpungeFrom(pdschIndices,gnb.NDLRB*12*3+1:gnb.NDLRB*12*5);
        PDSCHIndices = [PDSCHIndices pdschIndices.'];
    end
 
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
            pdschGrid = CSIRSmapping(pdschGrid,CSIRS(ind1,:),gnb,port);
             
            pdschGrid = DMRSmapping(pdschGrid,DMRS(ind2,:),port);

            txGrid((1:length(pdschGrid))+(BW-1)*length(pdschGrid),:,port) = pdschGrid;
        end
%         % OFDM modulation of associated resource elements
%         [txWaveform,txinfo] = OFDMModulate(gnb, squeeze(txGrid(:,:,port)));
%         txWave = [txWave txWaveform];
     end
       % OFDM modulation of associated resource elements
%         [txWaveform,txinfo] = OFDMModulate(gnb, txGrid);
%         txWave = [txWave txWaveform];
        txGrid_7500(6600*(band-1)+(1:6600),:,:)=txGrid; %7.5k 对应 2640x12x32
    end   % endof band
%        [txWaveform,txinfo] = OFDMModulate(gnb, txGrid);
       [txWaveform_7500,txinfo] = OFDMModulate_7500Hz(gnb, txGrid_7500);
%        [txWaveform,txinfo] = OFDMModulate_7500Hz_my(gnb, txGrid);
%         txWave = [txWave txWaveform];

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
       channeltype=2;
        Channel_model_apply_MUChannel_213_00;% 2019/02/13
        in = zeros(size(txWaveform_7500));  % 114816x32
        if channeltype == 3
            in(1:2120:end,:)=1;   %2114-2120
        elseif channeltype == 2
            in(1:210:end,:)=1;
        else
             in(1:40:end,:)=1; 
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
  
    if channeltype == 3
        out = hh(1:2120,:);
    elseif channeltype == 2
        out = hh(1:210,:);
    else
        out = hh(1:40,:);
    end
    
  
    rxWaveform_7500 = [];
    for antNo = 1:32
        tmp=conv(txWaveform_7500(:,antNo),out(:,antNo));  %114816x32   40x32
        rxWaveform_7500(:,antNo)=tmp(1:end-209);  %39,209,539      114816x32
    end
    
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
    sigPow=10*log10(var(txWaveform_7500));
    nVar=10.^(0.1.*(sigPow-SNRdB));
    rxWaveform_7500=AWGNChannel(rxWaveform_7500,nVar);
    %
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
        rxGrid_7500 = OFDMDemodulate_7500Hz(gnb, rxWaveform_7500);
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
%     end           
end   
toc
