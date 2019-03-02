function y = MIMOReceiver_MMSE(in, chEst, nVar)
%#codegen
% MIMO Receiver:
%   Based on received channel estimates, process the data elements
%   to equalize the MIMO channel. Uses the MMSE detector.
% noisFac = numLayers*diag(nVar);
noisFac = diag(nVar);
numData = size(in, 1);
y = complex(zeros(size(in)));
numTx = size(chEst,3);
%% MMSE receiver
for n = 1:numData
    iWn = eye(32);             % Orthonormal matrix
    h = chEst(n,:,:); 
    h = squeeze(h);
%     Q = (h*h' + noisFac)\h';
%     Q = (pinv(h)*h+ noisFac)\pinv(h);
    Q = pinv(h);
    x = in(n, :) * Q.';
    if numTx == 32
        tmp = x * iWn.';
    else
        tmp = x * iWn(:,1:12).';
    end
    y(n, :) = tmp;
end