function [lossSig,chanDelay] = helperApplyMUChannel(sig,prm,spLoss,varargin)
% Apply MIMO channel to input signal
%   Options include:
%       'CDL': NR5G CDL Channel
%       'Scattering': Phased Scattering MIMO channel
%       'MIMO': Comm MIMO channel
%
%   The channel is modeled with a fixed seed so as to keep the same channel
%   realization between sounding and data transmission. In reality, the
%   channel would evolve between the two stages. This evolution is modeled
%   by prepending the preamble signal to the data signal, to prime the
%   channel to a valid state, and then ignoring the preamble portion from
%   the channel output.

%   Copyright 2017 The MathWorks, Inc.

narginchk(3,4);
numUsers = prm.numUsers;
numTx = prm.numTx;
numRx = prm.numRx;
if nargin>3
    % preSig, for data transmission
    preSig = varargin{1}; 
    sigPad = [preSig; zeros(prm.numPadZeros,numTx); ...
              sig; zeros(prm.numPadZeros,numTx)];
else
    % No preSig, for sounding
    preSig = []; 
    sigPad = [sig; zeros(prm.numPadZeros,numTx)];
end

% Create independent channels per user
chan      = cell(numUsers,1);
chanDelay = zeros(numUsers,1);
lossSig   = cell(numUsers,1);

switch prm.ChanType
    case 'CDL'

        % Set the array geometry
        [isTxURA,expFactorTx,isRxURA,expFactorRx] = helperArrayInfo(prm);
        P = 1;      % P: no. of polarizations (fixed to 1)
        Mg = 1;     % Mg: no. of row array panels
        Ng = 1;     % Ng: no. of column array panels

        if isTxURA  % URA
            txM = expFactorTx;  % txM: no. of rows in Tx array
            txN = prm.numSTS;   % txN: no. of columns in Tx array
        else   % ULA
            txM = 1;
            txN = numTx;
        end
        
        % Create independent channels per user
        for uIdx = 1:numUsers

            % NR 5G CDL Channel model (TR 38.901)
            chan{uIdx} = nr5gCDLChannel('DelayProfile','CDL-A',...
                            'CarrierFrequency',prm.fc,....
                            'SampleRate',prm.chanSRate,...
                            'Seed', uIdx); 

            if numUsers>1 % For MU, use custom profile
                % Setup the delay profile (close to CDL-A with 4 clusters only)
                %   With low Doppler shift by default
                %   Make these user-settable
                chan{uIdx}.DelayProfile = 'Custom';
                chan{uIdx}.PathDelays = [0 9].*1e-8;    %  0.1207 0.176].*1e-7;
                chan{uIdx}.AveragePathGains = [0 -20];  % -2.2 -4];
                chan{uIdx}.AnglesAoD = [-178.1 -4.2];   %  4.2 -4.2];
                chan{uIdx}.AnglesAoA = [51.3 -152.7];   % -152.7 -152.7];
                chan{uIdx}.AnglesZoD = [50.2 93.2];     %  93.2 93.2];
                chan{uIdx}.AnglesZoA = [125.4 91.3];    %  91.3 91.3];
                chan{uIdx}.HasLOSCluster = false;
                chan{uIdx}.AngleSpreads = [5 11 3 3];
                chan{uIdx}.XPR = 10;
                chan{uIdx}.NumStrongestClusters = 0;    % No further splitting  
                chan{uIdx}.MaximumDopplerShift = 2;
            end
            
            % Setup the arrays
            chan{uIdx}.TransmitAntennaArray.Size = [txM,txN,P,Mg,Ng];
            chan{uIdx}.TransmitAntennaArray.Element = 'isotropic';

            if isRxURA(uIdx)
                rxM = expFactorRx(uIdx);    % rxM: numRows in Rx array
                rxN = prm.numSTSVec(uIdx);  % rxN: numColumns in Rx array
            else
                rxM = 1;
                rxN = numRx(uIdx);
            end
            chan{uIdx}.ReceiveAntennaArray.Size = [rxM,rxN,P,Mg,Ng];
            
            chanInfo = info(chan{uIdx});
            chanDelay(uIdx) = chanInfo.ChannelFilterDelay;

            % Apply channel
            fadeSig = chan{uIdx}(sigPad);

            % Remove the preamble, if present
            if ~isempty(preSig)
                fadeSig(1:(length(preSig)+prm.numPadZeros),:) = [];
            end
            
            % Apply path loss
            lossSig{uIdx} = fadeSig/sqrt(db2pow(spLoss(uIdx)));
            
        end
        
    case 'Scattering'
        % phased.ScatteringMIMOChannel
        %   No motion => static channel.
        
        % Tx & Rx Arrays
        [isTxURA,expFactorTx,isRxURA,expFactorRx] = helperArrayInfo(prm);
        
%         %   Specify spacing in direct units (meters)
%         if isTxURA % URA
%             txarray = phased.URA([expFactorTx,prm.numSTS], ...
%                 [0.5 0.5]*prm.lambda,'Element', ...
%                 phased.IsotropicAntennaElement('BackBaffled',false));
%         else % ULA
%             txarray = phased.ULA('Element', ...
%                 phased.IsotropicAntennaElement('BackBaffled',false),...
%                 'NumElements',numTx,'ElementSpacing',0.5*prm.lambda);
%         end         
          % Uniform Rectangular array
    antenna = phased.ShortDipoleAntennaElement(...
    'FrequencyRange',[2.9e9,3.5e9],...
    'AxisDirection','Z');
    txarray = phased.UCA('NumElements',32,'Radius',5000,'Element',antenna);
        
        % Create independent channels per user
        for uIdx = 1:numUsers

%             if isRxURA(uIdx) % URA
%                 rxarray = phased.URA([expFactorRx(uIdx),prm.numSTSVec(uIdx)], ...
%                     [0.5 0.5]*prm.lambda,'Element', ...
%                     phased.IsotropicAntennaElement);
%             else % ULA
%                 if numRx(uIdx)>1
%                     rxarray = phased.ULA('Element',phased.IsotropicAntennaElement, ...
%                         'NumElements',numRx(uIdx),'ElementSpacing',0.5*prm.lambda);
%                 else % numRx==1
%                     error(message('comm_demos:helperApplyMUChannel:invScatConf'));
%                     % only a single antenna, but ScatteringMIMOChannel doesnt accept this!
%                     % rxarray = phased.IsotropicAntennaElement;
%                 end                    
%             end
  rxarray = phased.ULA(numRx(uIdx), ...
        'ElementSpacing',0.5*prm.lambda, ...
        'Element',phased.IsotropicAntennaElement);

            Ns = 100;          % Number of scatterers

            % Place scatterers randomly in a circle from the center
            % posCtr = (prm.posTx+prm.posRx(:,uIdx))/2;

            % Place scatterers randomly in a sphere around the Rx
            %   similar to the one-ring model
            posCtr = prm.posRx(:,uIdx);
            radCtr = prm.mobileRanges(uIdx)*0.1;
            scatBound = [posCtr(1)-radCtr posCtr(1)+radCtr; ...
                         posCtr(2)-radCtr posCtr(2)+radCtr; ...
                         posCtr(3)-radCtr posCtr(3)+radCtr];
                       
            % Channel
            chan{uIdx} = phased.ScatteringMIMOChannel(...
                'TransmitArray',txarray,...
                'ReceiveArray',rxarray,...
                'PropagationSpeed',prm.cLight,...
                'CarrierFrequency',prm.fc,...
                'SampleRate',prm.chanSRate, ...
                'SimulateDirectPath',false, ...
                'ChannelResponseOutputPort',true, ...
                'TransmitArrayPosition',prm.posTx,...
                'ReceiveArrayPosition',prm.posRx(:,uIdx),...
                'NumScatterers',Ns, ...
                'ScattererPositionBoundary',scatBound, ...
                'SeedSource','Property', ...
                'Seed',uIdx);

            [fadeSig,~,tau] = chan{uIdx}(sigPad);
            chanDelay(uIdx) = floor(min(tau)*prm.chanSRate);

            % Remove the preamble, if present
            if ~isempty(preSig)
                fadeSig(1:(length(preSig)+prm.numPadZeros),:) = [];
            end

            % Path loss is included in channel
            lossSig{uIdx} = fadeSig;
            
        end
        
    case 'MIMO'

        % Create independent channels per user
        for uIdx = 1:numUsers

            % Using comm.MIMOChannel, with no array information
            chan{uIdx} = comm.MIMOChannel('MaximumDopplerShift',0, ...
                'SpatialCorrelation',false, ...
                'NumTransmitAntennas',numTx, ...
                'NumReceiveAntennas',numRx(uIdx),...
                'RandomStream','mt19937ar with seed', ...
                'Seed',uIdx, ...
                'SampleRate',prm.chanSRate);

            numBytesPerElement = 16;
            maxBytes = 8e9;
            if numTx*numRx*(length(sigPad))*numBytesPerElement > maxBytes
                % If requested sizes are too large, process symbol-wise
                fadeSig = complex(zeros(length(sigPad), numRx));
                symLen = prm.FFTLength+prm.CyclicPrefixLength;
                numSymb = ceil(length(sigPad)/symLen);
                for idx = 1:numSymb
                    sIdx = (idx-1)*symLen+(1:symLen).';
                    fadeSig(sIdx,:) = chan{uIdx}(sigPad(sIdx,:));
                end
            else
                fadeSig = chan{uIdx}(sigPad);
            end

            % Check derived channel parameters
            chanInfo = info(chan{uIdx});
            chanDelay(uIdx) = chanInfo.ChannelFilterDelay;        

            % Remove the preamble, if present
            if ~isempty(preSig)
                fadeSig(1:(length(preSig)+prm.numPadZeros),:) = [];
            end
            
            % Apply path loss
            lossSig{uIdx} = fadeSig/sqrt(db2pow(spLoss(uIdx)));
            
        end
        
end

end

% [EOF]
