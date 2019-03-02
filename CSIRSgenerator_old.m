function y = CSIRSgenerator_old(gnb, numTx)
%  LTE Cell-Specific Reference signal generation.
%   Section 6.10.1 of 3GPP TS 36.211 v10.0.0.
%   Generate the whole set per OFDM symbol, for 2 OFDM symbols per slot,
%   for 2 slots per subframe, per antenna port (numTx). 
%   This fcn accounts for the per antenna port sequence generation, while
%   the actual mapping to resource elements is done in the Resource mapper.
%#codegen
% Assumed parameters
nS = gnb.NSubframe;
NcellID = 0;        % One of possible 504 values
Ncp = 1;            % for normal CP, or 0 for Extended CP
NmaxDL_RB = gnb.NDLRB;    % largest downlink bandwidth configuration, in resource blocks
y = complex(zeros(NmaxDL_RB*5, 1,numTx));
l = [8 9 10 11];     % OFDM symbol idx in a slot for common first antenna port
% l = [1 5 8 12];     % OFDM symbol idx in a slot for common first antenna port
% Buffer for sequence per OFDM symbol
seq = zeros(size(y,1)*2, 1); % *2 for complex outputs
hSeqGen = comm.GoldSequence('FirstPolynomial',[1 zeros(1, 27) 1 0 0 1],...
                                'FirstInitialConditions', [zeros(1, 30) 1], ...
                                'SecondPolynomial', [1 zeros(1, 27) 1 1 1 1],...
                                'SecondInitialConditionsSource', 'Input port',... 
                                'Shift', 1600,...
                                'SamplesPerFrame', length(seq));
hInt2Bit = comm.IntegerToBit('BitsPerInteger', 31);

% Generate the common first antenna port sequences
for port = 1:numTx
    for i = 1 % slot wise
        for lIdx = 1:1 % symbol wise
            c_init = mod((2^10)*(14*nS/2+l(lIdx)+1)*(2*NcellID+1) + NcellID,2^31);
            % Convert to binary vector
            iniStates = step(hInt2Bit, c_init);
            % Scrambling sequence - as per Section 7.2, 36.211
            seq = step(hSeqGen, iniStates);
            % Store the common first antenna port sequences
            y(:, lIdx, numTx) = (1/sqrt(2))*complex(1-2.*seq(1:2:end), 1-2.*seq(2:2:end));
        end
    end
end
if numTx ==1
    y = squeeze(y);
end
