function H = IdChEst(prmLTEPDSCH, prmMdl, chPathG)
% Ideal channel estimation for LTE subframes
%
%   Given the system parameters and the MIMO channel path Gains, provide
%   the ideal channel estimates for the RE corresponding to the data.
%   Limitation - will work for path delays that are multiple of channel sample
%   time and largest pathDelay < size of FFT
%   Implementation based on FFT of channel impulse response
persistent hFFT; 
if isempty(hFFT) 
   hFFT = dsp.FFT; 
end 
% get parameters
numDataTones = prmLTEPDSCH.NSubcarriers; 
N                        = prmLTEPDSCH.Nfft;
cpLen0               = prmLTEPDSCH.CyclicPrefixLengths(1);
cpLenR               = prmLTEPDSCH.CyclicPrefixLengths(2);
slotLen               = (N*7 + cpLen0 + cpLenR*6);
% Get path delays
pathDelays = prmMdl.PathDelays;
% Delays, in terms of number of channel samples, +1 for indexing
sampIdx = round(pathDelays/(1/prmMdl.SampleRate)) + 1;
[~, numPaths, numTx, numRx] = size(chPathG);

H = complex(zeros(numDataTones, 14, numRx, numTx));
for i= 1:numTx
    for j = 1:numRx
        link_PathG = chPathG(:, :, i, j);
        % Split this per OFDM symbol
        g = complex(zeros(2*7, numPaths));
        for jj = 1:2 % over two slots
            % First OFDM symbol
            g((jj-1)*7+1, :) = mean(link_PathG((jj-1)*slotLen + (1:(N+cpLen0)), :), 1);
%               g((jj-1)*7+1, :) = link_PathG((jj-1)*slotLen + 1, :);
            % Next 6 OFDM symbols
            for k = 1:6
                g((jj-1)*7+k+1, :) = mean(link_PathG((jj-1)*slotLen+cpLen0+k*N+(k-1)*cpLenR + (1:(N+cpLenR)), :), 1);
%                   g((jj-1)*7+k+1, :) = link_PathG((jj-1)*slotLen+cpLen0+k*N+(k-1)*cpLenR + 1, :);
            end
        end
        hImp = complex(zeros(2*7, N));
        hImp(:, sampIdx) = g; % assign pathGains at sample locations
        % FFT processing
        h =  step(hFFT,hImp.'); 

        % Reorder, remove DC, Unpack channel gains
        H(:, :, j, i) = [h(N/2-numDataTones/2+1:N/2, :); h(N/2+1:N/2+numDataTones/2, :)];        
    end
end
