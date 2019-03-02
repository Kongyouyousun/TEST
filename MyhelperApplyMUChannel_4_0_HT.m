function [lossSig,chanDelay] = MyhelperApplyMUChannel_4_0_RA(sig,prm,spLoss,varargin)
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
%     sigPad = [sig; zeros(prm.numPadZeros,numRx)];  %原始为[sig; zeros(prm.numPadZeros,numTx)]
    sigPad =[sig; zeros(prm.numPadZeros,numTx)];
end

% Create independent channels per user
chan      = cell(numUsers,1);
chanDelay = zeros(numUsers,1);
lossSig   = cell(numUsers,1);

switch prm.ChanType
    case 'CDL'
        user_AoD=[-180 -157.5 -135 -112.5 -90 -67.5 -45 -22.5 0 22.5 45 67.5 90 112.5 135 157.5];
        user_AoA=[0 22.5 45 67.5 90 112.5 135 157.5 -180 -157.5 -135 -112.5 -90 -67.5 -45 -22.5];
        user_ZoD=90*ones(1,numUsers);
        user_ZoA=90*ones(1,numUsers);
        
        % Set the array geometry
        [isTxURA,expFactorTx,isRxURA,expFactorRx] = MyhelperArrayInfo_2_0(prm);% downlink ,是将helperArrayInfo直接复制过来的.
%         [isTxURA,expFactorTx,isRxURA,expFactorRx] = MyhelperArrayInfo(prm);% uplink
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
                            'NormalizePathGains',true,...
                            'SampleRate',prm.chanSRate,...
                            'Seed', uIdx); 

            if numUsers>1 % For MU, use custom profile
                % Setup the delay profile (close to CDL-A with 4 clusters only)
                %   With low Doppler shift by default
                %   Make these user-settable
                chan{uIdx}.DelayProfile = 'Custom';
                
%                 chan{uIdx}.PathDelays = [0 9].*1e-8;    %  0.1207 0.176].*1e-7;
%                 chan{uIdx}.AveragePathGains = [0 -20];  % -2.2 -4];
%                 chan{uIdx}.AnglesAoD = [-178.1 -4.2]+user_AoD(uIdx);   %  4.2 -4.2];
%                 chan{uIdx}.AnglesAoA = [51.3 -152.7]+user_AoA(uIdx);   % -152.7 -152.7];
%                 chan{uIdx}.AnglesZoD = [user_ZoD(uIdx) user_ZoD(uIdx)];     %  93.2 93.2];
%                 chan{uIdx}.AnglesZoA = [user_ZoA(uIdx) user_ZoA(uIdx)];    %  91.3 91.3];
%                 channeltype=prm.channeltype;
%                 if channeltype==1
%                  % RA模型
%                  chan{uIdx}.PathDelays = [0 0.2 0.4 0.6].*1e-6;  %delay of the RA channel
%                  chan{uIdx}.AveragePathGains = [0 -2 -10 -20];  %power of the RA channel
%                  chan{uIdx}.AnglesAoD = [-178.1 -4.2 4.2 178.1]+user_AoD(uIdx);   %  4.2 -4.2];
%                  chan{uIdx}.AnglesAoA = [51.3 -152.7 152.7 -51.3]+user_AoA(uIdx);   % -152.7 -152.7];
%                  chan{uIdx}.AnglesZoD = [user_ZoD(uIdx) user_ZoD(uIdx) user_ZoD(uIdx) user_ZoD(uIdx)];     %  93.2 93.2];
%                  chan{uIdx}.AnglesZoA = [user_ZoA(uIdx) user_ZoA(uIdx) user_ZoA(uIdx) user_ZoA(uIdx)];    %  91.3 91.3];
%                 elseif channeltype==2
%                 % TU 模型
%                 chan{uIdx}.PathDelays = [0 0.4 1.0 1.6 5 6.6].*1e-6;  %delay of the TU channel 
%                 chan{uIdx}.AveragePathGains = [-3 0 -3 -5 -2 -4];   %power of the TU  channel 
%                 chan{uIdx}.AnglesAoD = [-178.1 -95.5 -4.2 4.2 95.5 178.1]+user_AoD(uIdx);   %  4.2 -4.2];
%                 chan{uIdx}.AnglesAoA = [51.3 -30.7 -152.7 152.7 30.7 -51.3]+user_AoA(uIdx);   % -152.7 -152.7];
%                 chan{uIdx}.AnglesZoD = [user_ZoD(uIdx) user_ZoD(uIdx) user_ZoD(uIdx) user_ZoD(uIdx) user_ZoD(uIdx) user_ZoD(uIdx)];     %  93.2 93.2];
%                 chan{uIdx}.AnglesZoA = [user_ZoA(uIdx) user_ZoA(uIdx) user_ZoA(uIdx) user_ZoA(uIdx) user_ZoA(uIdx) user_ZoA(uIdx)];    %  91.3 91.3];
%                 else 
%                   % HT模型
%                                 % HT   模型             
%                 chan{uIdx}.PathDelays = [0 0.2 0.4 0.6 15 17.2 ].*1e-6;  %delay of the TU channel  15,17.2
%                 chan{uIdx}.AveragePathGains = [0 -2 -4 -7 -6 -12];   %power of the TU  channel  
%                 chan{uIdx}.AnglesAoD = [-178.1 -95.5 -4.2 4.2 95.5 178.1]+user_AoD(uIdx);   %  4.2 -4.2];
%                 chan{uIdx}.AnglesAoA = [51.3 -30.7 -152.7 152.7 30.7 -51.3]+user_AoA(uIdx);   % -152.7 -152.7];
%                 chan{uIdx}.AnglesZoD = [user_ZoD(uIdx) user_ZoD(uIdx) user_ZoD(uIdx) user_ZoD(uIdx) user_ZoD(uIdx) user_ZoD(uIdx)];     %  93.2 93.2];
%                 chan{uIdx}.AnglesZoA = [user_ZoA(uIdx) user_ZoA(uIdx) user_ZoA(uIdx) user_ZoA(uIdx) user_ZoA(uIdx) user_ZoA(uIdx)];    %  91.3 91.3];
%                 end

                %                 % RA模型
%                 chan{uIdx}.PathDelays = [0 0.2 0.4 0.6].*1e-6;  %delay of the RA channel 
%                 chan{uIdx}.AveragePathGains = [0 -2 -10 -20];  %power of the RA channel
%                 chan{uIdx}.AnglesAoD = [-178.1 -4.2 4.2 178.1]+user_AoD(uIdx);   %  4.2 -4.2];
%                 chan{uIdx}.AnglesAoA = [51.3 -152.7 152.7 -51.3]+user_AoA(uIdx);   % -152.7 -152.7];
%                 chan{uIdx}.AnglesZoD = [user_ZoD(uIdx) user_ZoD(uIdx) user_ZoD(uIdx) user_ZoD(uIdx)];     %  93.2 93.2];
%                 chan{uIdx}.AnglesZoA = [user_ZoA(uIdx) user_ZoA(uIdx) user_ZoA(uIdx) user_ZoA(uIdx)];    %  91.3 91.3];
                
%                  TU 模型
                % old:
%                 chan{uIdx}.PathDelays = [0 0.4 1.0 1.6 5 6.6].*1e-6;  %delay of the TU channel 
%                 chan{uIdx}.AveragePathGains = [-3 0 -3 -5 -2 -4];   %power of the TU  channel 
                % new
%                 chan{uIdx}.PathDelays = [0 0.4 1.0 1.6 2 2.5].*1e-6;  %delay of the TU channel
%                 chan{uIdx}.AveragePathGains = [-3 0 -3 -5 -2 -4];   %power of the TU  channel
%                 
%                 chan{uIdx}.AnglesAoD = [-178.1 -95.5 -4.2 4.2 95.5 178.1]+user_AoD(uIdx);   %  4.2 -4.2];
%                 chan{uIdx}.AnglesAoA = [51.3 -30.7 -152.7 152.7 30.7 -51.3]+user_AoA(uIdx);   % -152.7 -152.7];
%                 chan{uIdx}.AnglesZoD = [user_ZoD(uIdx) user_ZoD(uIdx) user_ZoD(uIdx) user_ZoD(uIdx) user_ZoD(uIdx) user_ZoD(uIdx)];     %  93.2 93.2];
%                 chan{uIdx}.AnglesZoA = [user_ZoA(uIdx) user_ZoA(uIdx) user_ZoA(uIdx) user_ZoA(uIdx) user_ZoA(uIdx) user_ZoA(uIdx)];    %  91.3 91.3];

                 % HT   模型  
                 % old :
                chan{uIdx}.PathDelays = [0 0.2 0.4 0.6 15 17.2 ].*1e-6;  %delay of the TU channel  15,17.2
                chan{uIdx}.AveragePathGains = [0 -2 -4 -7 -6 -12];   %power of the TU  channel 
                 % new: 16.5 x 30.72 < 512
%                 chan{uIdx}.PathDelays = [0 0.2 0.4 0.6 0.7 6.6 ].*1e-6;  %delay of the TU channel  15,17.2
%                 chan{uIdx}.AveragePathGains = [0 -2 -4 -7 -6 -12];   %power of the TU  channel 
                 % ITU Channel Model for Vehicular Test Environment, Channel A.             
%                 chan{uIdx}.PathDelays = [0 0.3 0.7 1.09 1.73 2 ].*1e-6;  %delay of the TU channel  15,17.2
%                 chan{uIdx}.AveragePathGains = [0 -1 -9 -10 -15 -20];   %power of the TU  channel 
    
%                 
                chan{uIdx}.AnglesAoD = [-178.1 -95.5 -4.2 4.2 95.5 178.1]+user_AoD(uIdx);   %  4.2 -4.2];
                chan{uIdx}.AnglesAoA = [51.3 -30.7 -152.7 152.7 30.7 -51.3]+user_AoA(uIdx);   % -152.7 -152.7];
                chan{uIdx}.AnglesZoD = [user_ZoD(uIdx) user_ZoD(uIdx) user_ZoD(uIdx) user_ZoD(uIdx) user_ZoD(uIdx) user_ZoD(uIdx)];     %  93.2 93.2];
                chan{uIdx}.AnglesZoA = [user_ZoA(uIdx) user_ZoA(uIdx) user_ZoA(uIdx) user_ZoA(uIdx) user_ZoA(uIdx) user_ZoA(uIdx)];    %  91.3 91.3];

                chan{uIdx}.HasLOSCluster = false;
                chan{uIdx}.AngleSpreads = [5 11 3 3];
                chan{uIdx}.XPR = 10;
                chan{uIdx}.NumStrongestClusters = 0;    % No further splitting  
                chan{uIdx}.MaximumDopplerShift = 10;   % 2019/02/16
%                 chan{uIdx}.MaximumDopplerShift = 0;   % 由10改为0; %2019/02/13
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
%             fadeSig = chan{uIdx}(sigPad(:,(uIdx-1)*2+(1:2)));
            fadeSig = chan{uIdx}(sigPad); % The number of columns of the signal input (32) must be
                                          % equal to the number of transmit antennas defined in TransmitAntennaArray (2).
            % Remove the preamble, if present
            if ~isempty(preSig)
                fadeSig(1:(length(preSig)+prm.numPadZeros),:) = [];
            end
            
            % Apply path loss
%             lossSig{uIdx} = fadeSig/sqrt(db2pow(spLoss(uIdx)));
             lossSig{uIdx} = fadeSig;  %没有路径损耗
           
        end
        
    case 'Scattering'
        % phased.ScatteringMIMOChannel
        %   No motion => static channel.
        % Tx & Rx Arrays
        [isTxURA,expFactorTx,isRxURA,expFactorRx] = MyhelperArrayInfo(prm); % Tx and Rx exchange
        
        %   Specify spacing in direct units (meters) 基站
        if isRxURA % URA
            rxarray = phased.URA([expFactorRx,prm.numSTS], ...
                [0.5 0.5]*prm.lambda,'Element', ...
                phased.IsotropicAntennaElement('BackBaffled',false));
        else % ULA
            rxarray = phased.ULA('Element', ...
                phased.IsotropicAntennaElement('BackBaffled',false),...
                'NumElements',numRx,'ElementSpacing',0.5*prm.lambda);
        end         
        
        % Create independent channels per user 用户
        for uIdx = 1:numUsers

            if isTxURA(uIdx) % URA
                txarray = phased.URA([expFactorTx(uIdx),prm.numSTSVec(uIdx)], ...
                    [0.5 0.5]*prm.lambda,'Element', ...
                    phased.IsotropicAntennaElement);
            else % ULA
                if numTx(uIdx)>1
                    txarray = phased.ULA('Element',phased.IsotropicAntennaElement, ...
                        'NumElements',numTx(uIdx),'ElementSpacing',0.5*prm.lambda);
                else % numRx==1
                    error(message('comm_demos:helperApplyMUChannel:invScatConf'));
                    % only a single antenna, but ScatteringMIMOChannel doesnt accept this!
                    % rxarray = phased.IsotropicAntennaElement;
                end                    
            end

            Ns = 100;          % Number of scatterers

            % Place scatterers randomly in a circle from the center
            % posCtr = (prm.posTx+prm.posRx(:,uIdx))/2;

            % Place scatterers randomly in a sphere around the Rx
            %   similar to the one-ring model
            posCtr = prm.posTx(:,uIdx);
            radCtr = prm.mobileRanges(uIdx)*0.1;
            scatBound = [posCtr(1)-radCtr posCtr(1)+radCtr; ...
                         posCtr(2)-radCtr posCtr(2)+radCtr; ...
                         posCtr(3) posCtr(3)];  %平面模型不需要Z轴
                       
            % Channel
            chan{uIdx} = phased.ScatteringMIMOChannel(...
                'TransmitArray',txarray,...
                'ReceiveArray',rxarray,...
                'PropagationSpeed',prm.cLight,...
                'CarrierFrequency',prm.fc,...
                'SampleRate',prm.chanSRate, ...
                'SimulateDirectPath',false, ...
                'ChannelResponseOutputPort',true, ...
                'ReceiveArrayPosition',prm.posRx,...
                'TransmitArrayPosition',prm.posTx(:,uIdx),...
                'NumScatterers',Ns, ...
                'ScattererPositionBoundary',scatBound, ...
                'SeedSource','Property', ...
                'Seed',uIdx);
            
            [fadeSig,~,tau] = chan{uIdx}(sigPad(:,(uIdx-1)*2+(1:2)));  % 每个用户数据依次输入
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
             lossSig{uIdx} = fadeSig;
%             lossSig{uIdx} = fadeSig/sqrt(db2pow(spLoss(uIdx)));
            
        end
        
end

end

% [EOF]
