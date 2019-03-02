clear all
clc
NFrames = 1;      % Number of 10ms frames
SNRIn = 60;  % SNR range

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
nRxAntennas = 4;         % Number of receive antennas at UE

% Set waveform type and OFDM numerology (SCS and CP type)
simulationParameters = [];                      % clear simulationParameters variable
simulationParameters.WaveformType = 'CP-OFDM';   % 'W-OFDM', 'F-OFDM' or 'CP-OFDM'
simulationParameters.SubcarrierSpacing = 60;    % 15,30,60,120,240,480 (kHz)
simulationParameters.CyclicPrefix = 'Extended';   % 'Normal' or 'Extended'
simulationParameters.UseDCSubcarrier = 'On';    % 'On' or 'Off'
simulationParameters.NTxAnts=1;
% Set number of RB and implicit number of transmit antennas
simulationParameters.NDLRB = 68;   % 20MHz at 15kHz subcarrier spacing
% simulationParameters.CellRefP =32;  % Set to 2 or 4 for 'TxDiversity' or 'CDD' transmission scheme

% DL-SCH/PDSCH parameters
simulationParameters.PDSCH.TargetCodeRate = 0.5;    % Code rate used to calculate transport block sizes
simulationParameters.PDSCH.CodingType = 'LDPC';     % Set to 'LDPC' or 'Turbo'
simulationParameters.PDSCH.PRBSet = (0:simulationParameters.NDLRB-1)';  % Fullband PDSCH allocation
simulationParameters.PDSCH.TxScheme = 'Port0';      % Set to 'Port0', 'TxDiversity' or 'CDD'
simulationParameters.PDSCH.Modulation = {'16QAM'};  % 'QPSK', '16QAM', '64QAM', '256QAM'
simulationParameters.PDSCH.CSI = 'On';              % Turn on LLR scaling after demodulation, set to 'On' or 'Off'
% simulationParameters.PDSCH.NLayers = 4;             % Applicable to CDD (2 or 4, depending on CellRefP)
ncw = 1;      % Number of active codewords (can be 2 for CDD)

gnb = lteRMCDL(simulationParameters,ncw); % LTE waveform generation parameterssimulationParameters
pdsch = gnb.PDSCH;                        % Separate out the PDSCH parameters
[~,~,rmcInfo] = lteRMCDLTool(gnb,[]);


% rvSequence = pdsch.RVSeq;
trBlkSizes = pdsch.TrBlkSizes;

% Set multi-mode waveform specific parameters
% WOLA specific parameters
gnb.Alpha = 0.0125;
% F-OFDM specific parameters
gnb.FilterLength = 513;
gnb.ToneOffset = 2.5;
waveformInfo = h5gOFDMInfo(gnb);
%% Propagation Channel Model Configuration
% Create the channel model object. Both CDL and TDL channel models are
% supported.

% Given the use of LTE PDSCH for transmission, use the LTE CellRefP
% parameter to define the number of PDSCH transmit antennas
% nTxAntennas = simulationParameters.CellRefP;
nTxAntennas = 1;

% Initialize variables used in the simulation and analysis
% Array to store the maximum throughput for all SNR points
maxThroughput = zeros(length(SNRIn),1); 
% Array to store the simulation throughput for all SNR points
simThroughput = zeros(length(SNRIn),1);

for snrIdx = 1:numel(SNRIn)
    txWave = [];
    CSIRSIndices = [];
    DMRSIndices = [];
    Peak = zeros(32,1);
    RMS = zeros(32,1);
    pdschSymbols = zeros(816*6,32,2);
    %产生八层PDSCH数据
    for layer = 1:32
        for BW = 1:2
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
            pdschInfo.G = 816*6;
            codedTrBlkSize = pdschInfo.G;   % Available PDSCH bit capacity
            trdata = randi([0,1],trBlk,1);
            codedTrBlock = h5gDLSCH(gnb, pdsch, codedTrBlkSize*4, trdata);
            
            % LTE PDSCH complex symbol generation
            pdschSym = ltePDSCH(gnb, pdsch, codedTrBlock);
            pdschSymbols(:,layer,BW) = pdschSym;
        end
    end
    
    for BW = 1:2
    pdschSymbols(:,:,BW) = squeeze(pdschSymbols(:,:,BW))*eye(32);
    end
    
    pdschIndices = 1:12*gnb.NDLRB*12;    
    %generate CSI-RS on 32 antenna
    CSIRS = CSIRSgenerator(gnb, 1);
    %generate DMRS on 12 antenna
    DMRS = DMRSgenerator(gnb, 1);
    
    %%映射参考信号
    for port = 1:32
        temp = [];
        for BW = 1:2            
            % PDSCH mapping in grid associated with PDSCH transmission period
            pdschGrid = lteDLResourceGrid(gnb);

            if BW == 1
                ind1 = 1:length(CSIRS)/2;
                ind2 = 1:length(DMRS)/2;
            else
                ind1 = ind1 + length(CSIRS)/2;
                ind2 = ind2 + length(DMRS)/2;
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
    end
    pdschIndices = ExpungeFrom(pdschIndices,gnb.NDLRB*12*7+1:gnb.NDLRB*12*11);
    pdschIndices = ExpungeFrom(pdschIndices,gnb.NDLRB*12*3+1:gnb.NDLRB*12*5);
    %映射PDSCH数据
    txGrid = zeros(gnb.NDLRB*12*2,12,32);
     for port = 1:32
        for BW = 1:2            
            % PDSCH mapping in grid associated with PDSCH transmission period
            pdschGrid = zeros(gnb.NDLRB*12,12);
            pdschGrid(pdschIndices) = pdschSymbols(:,port,BW);
            
            
            if BW == 1
                ind1 = 1:length(CSIRS)/2;
                ind2 = 1:length(DMRS)/2;
            else
                ind1 = ind1 + length(CSIRS)/2;
                ind2 = ind2 + length(DMRS)/2;
            end
            %CSI-RS mapping
            pdschGrid = CSIRSmapping(pdschGrid,CSIRS(ind1,:),gnb,port);
             
                pdschGrid = DMRSmapping(pdschGrid,DMRS(ind2,:),port);

            txGrid((1:length(pdschGrid))+(BW-1)*length(pdschGrid),:,port) = pdschGrid;
        end
        % OFDM modulation of associated resource elements
        [txWaveform,txinfo] = OFDMModulate_60kHz(gnb, squeeze(txGrid(:,:,port)));
        txWave = [txWave txWaveform];
     end
     hSpec1 = dsp.SpectrumAnalyzer('SampleRate',  122.88e6, ...
         'SpectrumType', 'Power density', 'PowerUnits', 'dBm', ...
         'RBWSource', 'Property',   'RBW', 60000,...
         'FrequencySpan', 'Span and center frequency',...
         'Span',  110e6, 'CenterFrequency', 0,...
         'SpectralAverages', 10, ...
         'Title', '发射信号频谱', 'YLimits', [-200 -0],...
         'YLabel', 'PSD');
     step(hSpec1,txWave(:,1));    
        
        %         % Write the sampling rate and chip rate to the configuration structure to
        %         % allow the calculation of ACLR parameters
        %         gnb.SamplingRate=txinfo.SamplingRate;
        %         gnb.UTRAChipRate = 3.84;              % UTRA chip rate in MCPS
        %         % Calculate ACLR measurement parameters
        %         aclr = my5g_hACLRParameters(gnb);
        %         % Apply required oversampling
        %         resampled = resample(txWaveform(:,1),aclr.OSR,1);
        %         % Calculate E-UTRA ACLR
        %         aclr = hACLRMeasurementEUTRA(aclr, resampled);
        %         % hACLRResults.m displays the ACLR and plots the adjacent channel powers.
        %         myhACLRResults(aclr);
        
        % %         % LTE case: filter CP-OFDM signal to meet eUTRA ACLR requirements
        % %         if (strcmpi(gnb.WaveformType,'CP-OFDM'))
        % %             txWaveform = cpofdmFilter(gnb,txinfo,txWaveform);
        % %         end
  
        %         % Pass data through channel model. Append zeros at the end of the
        %         % transmitted waveform to flush channel content. These zeroes take
        %         % into account any delay introduced in the channel. This is a mix
        %         % of multipath delay and implementation delay. This value may
        %         % change depending on the sampling rate, delay profile and delay
        %         % spread
        %         txWaveform = [txWaveform; zeros(maxChDelay, size(txWaveform,2))]; %#ok<AGROW>
        %         [rxWaveform,pathGains] = channel(txWaveform);
        %       COST 207
channeltype = 1;
Cost=cost_207(channeltype,1/122.88e6,0);
for antNo = 1:32
    rxWaveform(:,antNo)=filter(Cost,txWave(:,antNo));
end
        
        % %         %% MIMO Fading channel with AWGN
        % %         % Initialize channel time for each subframe
        % %         channel2.InitTime = gnb.NSubframe/1000;
        % %         %channel.Seed = idx;    % Set seed per subframe for larger variation
        % %         %  - Channel estimation averaging affected
        % %         % Pass data through channel model
        % %         rxWaveform = lteFadingChannel(channel2, txWaveform);
        %% Add AWG noise
    sigPow=10*log10(var(txWave));
    nVar=10.^(0.1.*(sigPow-SNRdB));
    rxWaveform=AWGNChannel(rxWaveform,nVar);
    
    % Add frequency offset impairment to received waveform
     foffset = 1600.0;         % Frequency offset in Hertz
     fprintf('add frequency error: %0.1fHz\n',foffset);
     t = (0:length(txWave)-1).'/122.88e6;
     rxWaveform = rxWaveform.* repmat(exp(1i*2*pi*foffset*t),1,32);
     %detect frequency error
     foffset_est = mylteFrequencyOffset_60000Hz(rxWaveform);
     fprintf('detect Frequency offset error: %0.2fHz\n',foffset_est);
     rxWaveform= mylteFrequencyCorrect(rxWaveform,foffset_est);
     %detect frequency error
     foffset_est = mylteFrequencyOffset_60000Hz(rxWaveform);
      fprintf('Frequency offset error after FrequencyCorrect: %0.2fHz\n',foffset_est);
        
                % Perform OFDM demodulation on the received data to recreate the
                % resource grid
                rxGrid = OFDMDemodulate_60kHz(gnb, rxWaveform);              
        lenGrid = size(rxGrid,1)/2;
        pdschRx = [];
        pdschTx = [];
        pdschRxEq = [];
        for BW = 1:2
            grid = rxGrid((1:lenGrid)+(BW-1)*lenGrid,:,:);
            txgrid = txGrid((1:lenGrid)+(BW-1)*lenGrid,:,:);
            % Get LTE PDSCH indices to extract the resource elements
%             pdschIndices = ltePDSCHIndices(gnb, pdsch, pdsch.PRBSet);
            
%             %DMRS  LS channel estimate       
%             hD = ChannelEstimate(grid,txgrid,DMRSIndices);
      
            %CSIRS LS channel estimate
            hD = ChannelEstimate(grid,txgrid,CSIRSIndices);
%             EsNO = 10^(SNRdB/10);
%             hD = ChEst_freq_mmseforCSIrs(EsNO, rmcInfo.NDLRB, rmcInfo.Nfft,chinfo.PathSampleDelays,grid,txgrid,CSIRSIndices);
            chEst = ExtChResponse(hD, pdschIndices); 
                
            % Get PDSCH resource elements from the received grid
            pdschRx = lteExtractResources(pdschIndices,grid);
            pdschTx = lteExtractResources(pdschIndices,txgrid);
            pdschRxEq = MIMOReceiver_MMSE(pdschRx, chEst,nVar);
            
            figure
            plot(pdschRx(:,1), 'o', 'MarkerEdgeColor', [0.75 0 0], ...
                'MarkerFaceColor', [1 0.25 0.25], 'MarkerSize',3);title('第1层数据均衡前','FontSize',8);
%             axis([-1 1 -1 1]);
            figure
            plot(pdschRxEq(:,1), 'o', 'MarkerEdgeColor', [0 0.75 0], ...
                'MarkerFaceColor', [1 0.25 0.25], 'MarkerSize',3);title('第1层数据均衡后','FontSize',8);
%             axis([-1 1 -1 1]);
            figure
            plot(pdschRx(:,4), 'o', 'MarkerEdgeColor', [0.75 0 0], ...
                'MarkerFaceColor', [1 0.25 0.25], 'MarkerSize',3);title('第4层数据均衡前','FontSize',8);
%             axis([-1 1 -1 1]);
            figure
            plot(pdschRxEq(:,4), 'o', 'MarkerEdgeColor', [0 0.75 0], ...
                'MarkerFaceColor', [1 0.25 0.25], 'MarkerSize',3);title('第4层数据均衡后','FontSize',8);
%             axis([-1 1 -1 1]);
        %
        %         % Perfect channel estimation, use the value of the path gains
        %         % provided by the channel
        %         % Get estimates for the transmission period of LTE PDSCH instance
        %         % (in terms of symbols, this is equivalent to an LTE subframe)
        %         gnb.TotSubframes = 1;
        %         estChannelGrid = my5g_PerfectChannelEstimate(gnb,channel,pathGains,offset);
        % %         estChannelGrid = IdChEst_OFDMRx(rxGrid,gnb,channel);
        %
        %         % Noise estimation directly from the noise realization
        %         noiseGrid = h5gOFDMDemodulate(gnb, noise(1+offset:end ,:));
        %         noiseEst = var(noiseGrid(:));
        %
        %         % Get LTE PDSCH indices to extract the resource elements
        %         pdschIndices = ltePDSCHIndices(gnb, pdsch, pdsch.PRBSet);
        %
        %         % Get PDSCH resource elements from the received grid
        %         [pdschRx, pdschHest] = lteExtractResources(pdschIndices, rxGrid, estChannelGrid);
        %
        %         % Decode LTE PDSCH physical channel
        %         % This uses a modified version of ltePDSCHDecode which creates
        %         % an LLR output which is compatible with the LDPC decoder used in
        %         % the NR DL-SCH decoding
        %         [dlschBits, rxSymbols] = h5gPDSCHDecode(gnb, pdsch, pdschRx, pdschHest, noiseEst);
        %
        %         %EVM
        %         eqgrid = lteEqualizeZF(rxGrid,estChannelGrid);
        %         pdschRxEq = lteExtractResources(pdschIndices,eqgrid*(10^(-gnb.PDSCH.Rho/20)));
        %         evm = lteEVM(pdschRxEq,pdschSymbols);
        %
        %         %constellation
        %         figure
        %         plot(pdschRx(:,1), 'o', 'MarkerEdgeColor', [0.75 0 0], ...
        %             'MarkerFaceColor', [1 0.25 0.25], 'MarkerSize',3);title('均衡前','FontSize',8);
        %         figure
        %         plot(pdschRxEq(:,1), 'o', 'MarkerEdgeColor', [0 0.75 0], ...
        %             'MarkerFaceColor', [1 0.25 0.25], 'MarkerSize',3);title('均衡后','FontSize',8);
        %
        %         hNRWaveformsPlotDoubleSidedSpectrum(f,P,gnb.NDLRB,gnb.SubcarrierSpacing,'Spectrum');
        %
        %         % Decode the DL-SCH transport channel
        %         decbits = h5gDLSCHDecode(gnb, pdsch, trBlk, dlschBits);
        %         [number,ratio] = biterr(decbits{1},trdata);
        %
        %         % Update starting symbol number of next PDSCH transmission
        %         gnb.NSymbol = gnb.NSymbol + size(pdschGrid,2);
        %         % Update count of overall number of PDSCH transmissions
        %         npdsch = npdsch + 1;
            end
end


%% Results
% xlabel('SNR (dB)'); ylabel('Throughput (%)'); grid on;
% title(sprintf('%s / %s (%dx%d) / %s / NDLRB=%d / SCS=%dkHz',...
%                pdsch.CodingType, pdsch.TxScheme, nTxAntennas, nRxAntennas, ...
%               gnb.WaveformType, gnb.NDLRB, gnb.SubcarrierSpacing));

% Estimate power spectrum of the provided waveform. This function estimates
% the spectrum of the input waveform. It returns the estimated power
% spectrum and the associated frequency locations.
function [P,f] = estimateSpectrum(waveform,samplingRate)

% Estimate spectrum
spectrumEstimator = dsp.SpectrumEstimator;
spectrumEstimator.SampleRate = samplingRate;
spectrumEstimator.FFTLengthSource = 'Property';
spectrumEstimator.FFTLength = 4096;
spectrumEstimator.SpectralAverages = length(waveform)/spectrumEstimator.FFTLength;
spectrumEstimator.Window = 'Flat Top';
spectrumEstimator.FrequencyRange = 'centered';

P = spectrumEstimator(waveform) * 1e3; % in dBm
f = getFrequencyVector(spectrumEstimator); % in Hz

end

% Plot the double sided spectrum P in the locations f. The inputs NDLRB and
% subcarrierSpacing are used for the figure annotations.
function hNRWaveformsPlotDoubleSidedSpectrum(f,P,NDLRB,subcarrierSpacing,figureTitle)

figure;
grid on;
hold on;
title(figureTitle);

% Plot power spectrum P and zoom on the edge of the allocated band
bw = NDLRB*12*subcarrierSpacing*1e3*1e-6;
hNRWaveformsPlotZoomedSpectrum(f*1e-6,P,bw/2*(-1.25),bw/2*1.25)

legend([num2str(NDLRB) ' RBs. Subcarrier spacing: ' num2str(subcarrierSpacing) ' kHz.'],'Location','best');

% Plot an arrow specifying the used bandwidth
v = axis;
ha = annotation('doublearrow');
ha.Parent = gca;
ha.X = [-1 1]*bw*0.5;
ha.Y = (v(4)-(v(4)-v(3))*0.9)*[1 1];
bwLabel = ['BW: ' num2str(bw,'%0.1f')];
arrowLabel = [num2str(bwLabel) ' MHz'];
t = text(0,(v(4)-(v(4)-v(3))*0.95),arrowLabel);
t.HorizontalAlignment = 'Center';

end

% Plot the spectrum P in the locations f and zoom into the [xmin xmax]
% range.
function hNRWaveformsPlotZoomedSpectrum(f,P,xmin,xmax)

plot(f,10*log10(P));
v = axis;
axis([xmin xmax v(3) v(4)]);
xlabel('Frequency (MHz)')
ylabel('dBm');

end