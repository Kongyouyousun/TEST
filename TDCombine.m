function y = TDCombine(in, chEst, numTx, numRx)
% LTE transmit diversity combining
%   SFBC and SFBC with FSTD.
persistent hTDDec;
if isempty(hTDDec)
    % OSTBC combiner - always numTx = 2
    hTDDec = comm.OSTBCCombiner('NumTransmitAntennas', 2, ...
        'NumReceiveAntennas', numRx);
end
inLen = size(in, 1);
Index=(2:2:inLen)';
switch numTx
    case 1
        y=in;
    case 2   % For 2TX - SFBC
        in = sqrt(2) * in; % Scale
        y = step(hTDDec, in,chEst);
        % ST to SF transformation.
        % Apply blockwise correction for 2nd symbol combining
        y(Index) = -conj(y(Index));
    case 4   % For 4Tx - SFBC with FSTD
        in = sqrt(2) * in; % Scale
        H = complex(zeros(inLen, 2, numRx));
        idx12 = ([1:4:inLen; 2:4:inLen]); idx12 = idx12(:);
        idx34 = ([3:4:inLen; 4:4:inLen]); idx34 = idx34(:);
        H(idx12, :, :) = chEst(idx12, [1 3], :);
        H(idx34, :, :) = chEst(idx34, [2 4], :);
        y = step(hTDDec, in, H);
        % ST to SF transformation.
        % Apply blockwise correction for 2nd symbol combining
        y(Index) = -conj(y(Index));
end

