%%
% 信道用helperApplyMUChannel，上行链路
%

% Multi-user system with single/multiple streams per user
prm.numUsers = 16;                 % Number of users
prm.numSTSVec = 2*ones(1,16);        % Number of independent data streams per user
prm.numSTS = sum(prm.numSTSVec);  % Must be a power of 2
prm.numTx = prm.numSTS;         % Number of BS receive antennas (power of 2)
prm.numRx = prm.numSTSVec;      % Number of transmit antennas, per user (any >= numSTSVec)

% Each user has the same modulation
% prm.bitsPerSubCarrier = 4;   % 2: QPSK, 4: 16QAM, 6: 64QAM, 8: 256QAM
% prm.numDataSymbols = 10;     % Number of OFDM data symbols

% MS positions: assumes BS at origin
%   Angles specified as [azimuth;elevation] degrees
%   az in range [-180 180], el in range [-90 90], e.g. [45;0]
maxRange = 1000;            % all MSs within 1000 meters of BS
prm.mobileRanges = 100 *ones(1,16);%2019/02/17 14:36
% prm.mobileRanges = 50*ones(1,16);% 2019/02/15
azimuth=-180:22.5:180;
zzimuth=zeros(1,16);
prm.mobileAngles = [azimuth(1:16);zzimuth]; % 2019/02/15
% prm.mobileAngles = [azimuth(2:17);zzimuth];% 做了修改.2019/02/15前
% prm.mobileAngles = [azimuth(1:16)-180;azimuth(1:16)-90];;%做了修改.

prm.fc = 3.5e9;               % 3 GHz system.   2019/02/15前
% prm.fc = 3e9;                  % 3 GHz system.   2019/02/15
% prm.chanSRate = 30.72e6;       % Channel sampling rate, 100 Msps
prm.chanSRate = 122.88e6;       % Channel sampling rate, 100 Msps
prm.ChanType = 'CDL'; % Channel options: 'Scattering', 'MIMO','CDL'

numSTS = prm.numSTS;
numTx = prm.numTx;
numRx = prm.numRx;

prm.cLight = physconst('LightSpeed');
prm.lambda = prm.cLight/prm.fc;

% % Get transmit and receive array information
% [isTxURA,expFactorTx,isRxURA,expFactorRx] = helperArrayInfo(prm);

% Receive antenna array definition
%   Array locations and angles
prm.posTx = [0;0;0];         % BS array position, [x;y;z], meters

    % Uniform Rectangular array
    antenna = phased.ShortDipoleAntennaElement(...
    'FrequencyRange',[2.9e9,3.5e9],...
    'AxisDirection','Z');
    txarray = phased.UCA('NumElements',32,'Radius',maxRange,'Element',antenna);
%     txarray = phased.PartitionedArray(...
%         'Array',array,...
%         'SubarraySelection',eye(numTx),'SubarraySteering','Custom');
% pattern(rxarray,prm.fc,[-180:180],0,'Type','powerdb',...
%     'CoordinateSystem','polar',...
%     'Normalize',true);

prm.posTxElem = getElementPosition(txarray)/prm.lambda;

spLoss = zeros(prm.numUsers,1);
prm.posRx = zeros(3,prm.numUsers);
toRxAng = zeros(2,prm.numUsers);

for uIdx = 1:prm.numUsers    
    % Uniform Linear array
    rxarray = phased.ULA(numRx(uIdx), ...
        'ElementSpacing',0.5*prm.lambda, ...
        'Element',phased.IsotropicAntennaElement);
    prm.posRxElem = getElementPosition(rxarray)/prm.lambda;
    
    % Mobile positions
    [xRx,yRx,zRx] = sph2cart(deg2rad(prm.mobileAngles(1,uIdx)), ...
                             deg2rad(prm.mobileAngles(2,uIdx)), ...
                             prm.mobileRanges(uIdx));
    prm.posRx(:,uIdx) = [xRx;yRx;zRx];
    [toRxRange,toRxAng] = rangeangle(prm.posTx,prm.posRx(:,uIdx));
    spLoss(uIdx) = fspl(toRxRange,prm.lambda);
end

% Account for channel filter delay
prm.numPadZeros = 0;
% Apply a spatially defined channel to the transmit signal
% [rxSig,chanDelay] = MyhelperApplyMUChannel(txWaveform,prm,spLoss);