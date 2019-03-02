function reGrid = OFDMDemodulate_60kHz(gnb,waveform,varargin)
    
    if (nargin==3)
      cpFraction = varargin{1};
    else
      cpFraction = 0.55;
    end
    
    % Get dimensionality information derived from the parameters
    info = h5gOFDMInfo(gnb);
    info.NSubcarriers = info.NSubcarriers*2;
    info.Nfft = 2048;
    % Cache the main dims
    nSC = info.NSubcarriers;
    nFFT = info.Nfft;
    cpLengths = 512*ones(1,12);
    symbolsPerSubframe = 48;
    samplesPerSubframe =122880; %info.SamplesPerSubframe;       
           
    % Get dimensional information derived from the resource grid
    nAnts = size(waveform,2);
    
    % Initial symbol number of input data, wrt 1ms subframe
    initialsymbol = 0;  
    if isfield(gnb,'NSymbol')
        initialsymbol = mod(gnb.NSymbol,symbolsPerSubframe);
    end  
    
    % Calculate the number of complete OFDM symbols associated with the grid,
    % accounting for the starting symbol number of the data
    % Number of samples in the input waveform
    nsamples = size(waveform,1);
    % Number of whole subframes in the waveform 
    nwsubframes = fix(nsamples/samplesPerSubframe);
    % Number of remainder samples 
    remsamples = mod(nsamples,samplesPerSubframe);
    % From the remainder samples, subtract the lengths of the individual symbol in a subframe, 
    % starting from the initial symbol associated with the beginning of the input waveform
    sadj = cumsum([remsamples -circshift(cpLengths+nFFT,-initialsymbol)]);    
    % Find the position of the first negative value, to give the number of remainder symbols
    ns = find(sadj(2:end)<0,1);
    % Calculate the total number of OFDM symbols associated with the input waveform 
    nSymbols = nwsubframes*symbolsPerSubframe + ns-1;
     nSymbols = 12;
    
    % Create an empty output resource grid of the required size  
    dims = [nSC nSymbols nAnts];
    reGrid = zeros(dims);
       
    % Waveform dependent initialization
    % Initial sample offset to start of symbol portions of interest
    offset = 0;
    wola = isfield(gnb,'WaveformType') && strcmpi(gnb.WaveformType,'W-OFDM');
    if wola
        
        N1 = info.N1;
        N2 = info.N2;
        TL = nFFT + N1 + N2;                       % Overall length of windowed W-OFDM symbol (not accounting for overlap - it's the window length)       
        o = N1 + N2 - cpLengths(initialsymbol+1);  % Positive if overlap, negative if no overlap but gap between symbols
      
        offset = -o;    
        offset = offset+N2;
     
        if isfield(gnb,'WindowCoeffs') && ~isempty(gnb.WindowCoeffs)
           window = gnb.WindowCoeffs;
           if length(gnb.window) ~= TL
               error('lte:error','The W-OFDM window length must be equal to %d',TL);
           end
        else
           window = sqrt(raised_cosine_window(nFFT,N1+N2));        
        end 
    end
     
    % If filtering was applied to the waveform then match filter before
    % demodulation
    fofdm = isfield(gnb,'WaveformType') && strcmpi(gnb.WaveformType,'F-OFDM');
    if fofdm 
        waveform = fofdm_filter(gnb,nFFT,waveform);
    end

    % Should the DC subcarrier be used for transmission or not
    includezsc = ~(isfield(gnb,'UseDCSubcarrier') && strcmpi(gnb.UseDCSubcarrier,'Off'));
        
    % Calculate position of the first active subcarrier in the FFT output,
    % according to the FFT size and number of actuve subcarriers
    firstActiveSC = (nFFT/2) - nSC/2 + 1;

    % Demodulate all the symbols within the input data
    for symbol = initialsymbol+(0:nSymbols-1) 

        if wola
            % For W-OFDM
            %
            % Get the effective CP of *next* symbol (needed to move
            % to next symbol at the end of the loop)
            cpLength = cpLengths(mod(symbol+1,length(cpLengths))+1);

            % W-OFDM processing steps
            % - Capture the windowed overall range
            % - Fold the bordering CP/CS extension back into the internal ranges
            % - Extract the combined portion for the FFT

            % Extract full extended symbol
            wdata = waveform(mod(offset:offset+TL-1,size(waveform,1))+1,:).*window;            
            % Internal N1 (CP) region
            wdata(end-N2-N1+1:end-N2,:) =  wdata(end-N2-N1+1:end-N2,:) + wdata(1:N1,:);               
            % Internal N2 (CS) region
            wdata(N1+1:N1+N2,:) = wdata(N1+1:N1+N2,:) + wdata(end-N2+1:end,:);

            % Extract central portion 
            symwaveform = wdata(N1+1:N1+nFFT,:);

            % No phase correction required for W-OFDM
            phaseCorrection = 1;              

        else
            % For CP-OFDM/F-OFDM
            %
            % Get cyclic prefix length in samples for the *current* symbol
            cpLength = cpLengths(mod(symbol,length(cpLengths))+1);

            % Position the FFT part of the way through the cyclic
            % prefix. The value of cpFraction should ensure that the 
            % reception is unaffected by the windowing applied in the 
            % lteOFDMModulate function. Default is 55%                
            fftStart = fix(cpLength*cpFraction);                   

            symwaveform = waveform(offset+fftStart+(1:nFFT),:);

            % Create vector of phase corrections, one per FFT sample,
            % to compensate for FFT being performed away from zero phase
            % point on the original subcarriers
            idx = 0:nFFT-1;
            phaseCorrection = repmat(exp(-1i*2*pi*(cpLength-fftStart)/nFFT*idx)',1,size(waveform,2));

        end

        % Extract the appropriate section of input waveform, perform the
        % FFT and apply the phase correction
        fftOutput = fftshift(fft(symwaveform).*phaseCorrection,1)./(nFFT/sqrt(nSC));

        % Extract the active subcarriers for each antenna from the 
        % FFT output, removing the DC subcarrier if unused
        activeSCs = fftOutput(firstActiveSC:firstActiveSC+nSC-(includezsc~=0),:);
        if ~includezsc
            activeSCs(nSC/2+1,:) = [];
        end

        % Assign the active subcarriers into the appropriate column
        % of the received grid, for each antenna
        reGrid(OFDMSymbolIndices(reGrid,symbol-initialsymbol+1)) = activeSCs;

        % Update counter of overall symbol number and position in
        % the input waveform
        offset = offset+nFFT+cpLength;
    end
    
end

%OFDMSymbolIndices Generates indices for a given OFDM symbol of a resource array.
%   IND = OFDMSymbolIndices(GRID,SYMBOL) gives the indices for OFDM symbol number
%   SYMBOL of the resource array GRID, in per-antenna linear indices form 
function ind = OFDMSymbolIndices(reGrid,symbol)
    nSCs = size(reGrid,1);
    nAnts = size(reGrid,3);
    firstSCs = sub2ind(size(reGrid),ones(1,nAnts),repmat(symbol,1,nAnts),1:nAnts);
    ind = repmat(firstSCs,nSCs,1)+repmat((0:nSCs-1).',1,nAnts);
end

% Raised Cosine window creation; creates a window function of length n+w
% with raised cosine transitions on the first and last 'w' samples, where
% 'w' is the "number of time-domain samples (at each symbol edge) over 
% which windowing and overlapping of OFDM symbols is applied" as described
% in the product help
function p = raised_cosine_window(n,w)
    
    p = 0.5*(1-sin(pi*(w+1-2*(1:w).')/(2*w)));
    p = [p; ones(n-w,1); flipud(p)];
    
end

% Apply filtering to input waveform
function waveform = fofdm_filter(gnb,nFFT,waveform)
       
    % Manage the parameter set associated with the filter and resulting
    % coefficients to minimise the need to recreate them
    persistent existingConfig;
    persistent fnum;
    if isempty(existingConfig) 
        existingConfig = [-1 -1 -1];
    end
      
    % Half the filter length
    halfFilt = floor(gnb.FilterLength/2);
     
    % If there is a change to any of the filtering parameters then 
    % redesign the filter  
    currentConfig = [gnb.NDLRB, gnb.FilterLength, gnb.ToneOffset];

    if any(existingConfig ~= currentConfig)
        existingConfig = currentConfig;
       
        % Filter design  
        nSC = gnb.NDLRB*12;
        W = nSC;  % Number of data subcarriers

        n = -halfFilt:halfFilt;  

        % Sinc function
        pb = sinc((W+2*gnb.ToneOffset).*n./nFFT);

        % Sinc truncation window
        w = (0.5*(1+cos(2*pi.*n/(gnb.FilterLength-1)))).^0.6;

        % Normalized lowpass filter coefficients
        fnum = (pb.*w)/sum(pb.*w);
    end
     
    % If no cyclic extension required then just filter
    if isfield(gnb,'CyclicExtension') && strcmpi(gnb.CyclicExtension,'off')    
        % Filter
        waveform = filter(fnum, 1, waveform);
        
        % Remove the filter delay so that the resulting OFDM symbols have
        % the same temporal alignment as other OFDM variants for subsequent
        % processing
        waveform = circshift(waveform,-halfFilt);
    else
        % Create a cyclic, 'loopable' waveform by cyclically extending input,
        % filtering then remove leading transient part
        ewaveform = [waveform(end-gnb.FilterLength:end,:); waveform];

        % Filter then realign
        efwaveform = filter(fnum, 1, ewaveform);
        waveform = efwaveform(end-size(waveform,1)+1:end,:);

        % Remove filter delay so that filtered OFDM symbols have the same
        % temporal alignment as other OFDM variants
        waveform = circshift(waveform,-halfFilt);
    end
   
end
