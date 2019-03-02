function [waveform,info] = OFDMModulate_new_20190216_1555(gnb,grid,varargin)

    if (~isnumeric(grid))
        error('lte:error','The input resource grid must be a numeric array.');
    end

    % Get dimensionality information derived from the parameters
    info = h5gOFDMInfo(gnb);% 30.72MHz
    
%     info.NSubcarriers = info.NSubcarriers*5;
%     info.Nfft = 8192;
%     info.SamplingRate=122880000;
%     info.CyclicPrefixLengths=info.CyclicPrefixLengths*4;
%     info.SymbolLengths=info.SymbolLengths*4;
%     info.SamplesPerSubframe=info.SamplesPerSubframe*4;

    % Cache the main dims
    nSC = info.NSubcarriers*2; %7.5k
%     nFFT = info.Nfft;
    nFFT = info.Nfft*2; % 7.5k
    channeltype=gnb.mychannelType;
%     %     cpLengths = info.CyclicPrefixLengths; % old
%     if channeltype ==1      % new
%         cpLengths = double(info.CyclicPrefixLengths); %160,144....
%     else
%         cpLengths = 512*ones(1,12);  % 1024*ones(1,14)  , 信道2,3都用拓展CP
%         %          cpLengths = double(info.CyclicPrefixLengths);
%     end
%     cpLengths = 512*ones(1,12); % 2019/02/18 16:43
    cpLengths = 1024*ones(1,12); % 2019/02/18 16:43 % 7.5k
    
    symbolsPerSubframe = info.SymbolsPerSubframe;
    SamplesPerSubframe = nFFT*12+sum(cpLengths);%7.5k
    % Get dimensional information derived from the resource grid
    nSymbols = size(grid,2);
    nAnts = size(grid,3);
    
    % Initial symbol number of input wrt subframe
    initialsymbol = 0;  
    if isfield(gnb,'NSymbol')
        initialsymbol = mod(gnb.NSymbol,symbolsPerSubframe);
    end
    
    % Handle optional CP-OFDM windowing argument
    % N is the number of windowed samples
    if (nargin==3)
        N = varargin{1};
        mwltelibrary('validateLTEParameters',setfield(struct(),'Windowing',N),'Windowing');
    else
        N = double(info.Windowing);
    end
    if (~isscalar(N))
        error('lte:error','The Windowing parameter must be an even scalar integer.');
    end
    if (N>(nFFT-info.CyclicPrefixLengths(1)))
        error('lte:error','For the Windowing parameter the value (%d) must be less than or equal to %d (the IFFT size (%d) minus the longest cyclic prefix length (%d))',N,nFFT-info.CyclicPrefixLengths(1),nFFT,info.CyclicPrefixLengths(1));
    end
    if (mod(N,2)~=0)
        error('lte:error','For the Windowing parameter the value (%d) must be even.',N);
    end
    info.Windowing = N;
    
    % Index of first subcarrier in IFFT input
    firstSC = (nFFT/2) - nSC/2 + 1;
    
    % Number of active subcarriers in IFFT
    if (~any(size(grid,1)==[0 nSC info.Nfft]))
        error('lte:error','The input resource grid must contain a whole number of resource blocks i.e. number of rows must be an integer multiple of 12.');
    end
    
    % Handle various empty input cases
    if(size(grid,1)==0)
        nSymbols = 0;
    end
    if (nSymbols==0)
       head = 0;
       N = 0;
    end
    
    % Waveform dependent initialization
    wola = isfield(gnb,'WaveformType') && strcmpi(gnb.WaveformType,'W-OFDM');
    fofdm = isfield(gnb,'WaveformType') && strcmpi(gnb.WaveformType,'F-OFDM');
    % If F-OFDM then zero the amount of OFDM symbol windowing performed
    if fofdm
        info.Windowing = 0;
        N = 0;
    end
    
    % Should the DC subcarrier be used for transmission or not
    includezsc = ~(isfield(gnb,'UseDCSubcarrier') && strcmpi(gnb.UseDCSubcarrier,'Off'));
       
    % Some calculations still need this (those which worked with the number of frames associated with data)
    nSubframes = ceil(nSymbols/symbolsPerSubframe);
    
    if wola 
        
        % Perform W-OFDM specific initialization
        N1 = info.N1;               % Cyclic prefix part
        N2 = info.N2;               % Cyclic suffix part
        TL = nFFT + N1 + N2;        % Overall length of windowed W-OFDM symbol (not accounting for overlap - it's the window length)  
        % Pre-calculate windowing; there is a only a single windows required since it's independent of the CP length
        if isfield(gnb,'WindowCoeffs') && ~isempty(gnb.WindowCoeffs)
            window0 = gnb.WindowCoeffs;
            if length(gnb.window0) ~= TL
                error('lte:error','The W-OFDM window length must be equal to %d',TL);
            end
        else
            window0 = sqrt(raised_cosine_window(nFFT,N1+N2));        
        end  
        window1 = window0;      
        
        % Extension periods (CP and CS) for symbols
        exLengths = [N1*ones(size(cpLengths));N2*ones(size(cpLengths))];
               
        % Inter-symbol distances
        % Current stride length is the distance between the current symbol
        % and the start of the next one
        strideLengths = circshift(cpLengths,-1) + nFFT;
        
        % N1 + N2 is combined symbol extension
        % cpLength is required symbol extension 
        % Overlap is the amount of overlap required to achieve cpLength
        % Initialize starting position of first extended symbol in output waveform
        o = N1 + N2 - cpLengths(1);  % Overlap for first symbol in subframe: Positive if true overlap, negative if no overlap but instead an explicit gap between symbols
                    
        % Where to place the start of the first extended symbol       
        pos = -o;  % Start the first symbol output location such that there is none of the previous symbol present i.e. no overlap part
       
        % General calculation of the midpoint locations between consecutive symbols for any initial pos
        % The mid-point is in the middle of the overlap part 
        jumps = nFFT + fix((cpLengths)/2) + ceil(circshift(cpLengths,-1)/2);
        mp1 = pos + N1 - fix(cpLengths(1)/2) + 1;
        info.Midpoints = cumsum([mp1 repmat(jumps,[1 nSubframes])]);
        
        % Overlap between windowed, extended symbols
        info.WindowOverlap = N1+N2 - cpLengths;        
    
        % Re-adjust for the first actual symbol at the input
        pos = -(N1 + N2 - cpLengths(initialsymbol+1));
        
    else
        
        % Pre-calculate windowing; there are two different windows required:
        % one for the first OFDM symbol of a slot and one for other OFDM
        % symbols, because the first OFDM symbol of a slot has a different
        % cyclic prefix length.    

        % Window size will be nFFT + N + CP length
        window0 = raised_cosine_window(nFFT+cpLengths(1),N);
        window1 = raised_cosine_window(nFFT+cpLengths(2),N);
        
        % Extension periods (prefix, zero suffix) for symbols, accounting 
        % for any required for the windowing
        exLengths = [cpLengths+N; zeros(size(cpLengths))];

        % Inter-symbol distances
        % Current stride length is the distance between the current symbol
        % and the start of the next one
        strideLengths = cpLengths + nFFT;
        
        % Amount of overlap between extended symbols
        o = N;
   
        % Initialize starting position of first extended symbol in output waveform
        % With CP-OFDM, the windowing part will warp around into the end of
        % last symbol in the grid block
        pos = -N;
        
        % Midpoints between the extended symbols
        info.Midpoints = cumsum([1+pos+fix(o/2) repmat(strideLengths,[1 nSubframes])]); 
                   
        % Overlap between windowed, extended symbols           
        info.WindowOverlap = N*ones(size(cpLengths));
    
    end
      
    % Create storage for the returned waveform
    nsamples = sum(cpLengths(mod(initialsymbol + (0:nSymbols-1),symbolsPerSubframe)+1)) + nFFT*nSymbols;
    waveform = zeros(nsamples,nAnts);

    % Modulate each OFDM symbol of the resource grid
    for i = 1:nSymbols

        currentsymbol = initialsymbol+i-1;
        
        % Create IFFT input (map subcarriers)
        if (size(grid,1)==nFFT)
            ifftin = squeeze(grid(:,i,:));
        else
            ifftin = zeros(nFFT,nAnts);
            ifftin(firstSC+(0:nSC/2-1),:) = grid(1:nSC/2,i,:);
            ifftin(firstSC+nSC/2+(includezsc==0)+(0:nSC/2-1),:) = grid(nSC/2+1:end,i,:);
        end

        % Perform IFFT
        iffout = ifft(fftshift(ifftin,1)).*(nFFT/sqrt(nSC));
%             iffout = ifft(fftshift(ifftin,1)).*(8192/sqrt(6600));

        % Add cyclic extension to the symbol 
        exLength = exLengths(:,mod(currentsymbol,length(exLengths))+1);   % Extension lengths (CP/CS) for current symbol
        stride = strideLengths(mod(currentsymbol,length(strideLengths))+1);
        
        % Create extended symbol
        extended = [iffout(end-(exLength(1))+1:end,:); iffout; iffout(1:(exLength(2)),:)];
               
        % Perform windowing, using the appropriate window (first OFDM symbol
        % of half a subframe (0.5ms))
        if (mod(currentsymbol,symbolsPerSubframe/2)==0)
            windowed = extended .* window0;
        else
            windowed = extended .* window1;
        end

        % Perform overlapping and creation of output signal. Note that with
        % windowing the signal "head" gets chopped off the start of the
        % waveform and finally superposed on the end. This means the
        % overall signal can be seamlessly looped when output from an
        % arbitrary waveform generator
        if (i==1)
            % For the first OFDM symbol, chop the windowed "head" (which 
            % will be added to the final waveform end) and then output the
            % rest (rotate around the waveform)
            p = max(-pos,0);
            head = windowed(1:p,:);                 % Part of leading edge of extended symbol which overlaps with previous one (which here will be the last symbol)
            L = size(windowed,1) - size(head,1);    % Length of extended symbol which doesn't overlap with previous symbol
            p = max(pos,0);
            waveform(p+(1:L),:) = windowed(size(head,1)+1:end,:);  % Write non-overlapped (on leading edge) part of waveform to output                      
        else
            % For subsequent OFDM symbols, add the windowed part to the end
            % of the previous OFDM symbol (overlapping them) and then output
            % the rest; 'pos' points to the end of the previous OFDM symbol
            % i.e. the start of the current one, so merge it from N samples
            L = size(windowed,1);
            waveform(pos+(1:L),:) = waveform(pos+(1:L),:) + windowed;            
        end

        % Update 'pos' to point to the start of the next extended OFDM
        % symbol in the output waveform
        pos = pos + stride;         
    end

    % Finally, overlap the "head" with the very end of the signal
    waveform(end-size(head,1)+1:end,:) = waveform(end-size(head,1)+1:end,:) + head;
    
    % If W-OFDM then a final rotation is required to place the effective CP
    % of the first symbol at the beginning of the output waveform
    if wola
        waveform = circshift(waveform,N2,1);
        info.Midpoints = info.Midpoints + N2;
    end
    
    % If F-OFDM then filtering is required
    if fofdm 
        waveform = fofdm_filter(gnb,nFFT,waveform);
    end
    waveform=waveform/3.8;
end

% Raised Cosine window creation; creates a window function of length n+w
% with raised cosine transitions on the first and last 'w' samples, where
% 'w' is the "number of time-domain samples (at each symbol edge) over 
% which windowing and overlapping of OFDM symbols is applied" as described
% in the product help.
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
        nSC = gnb.NDLRB*12*5;
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
    else       
        % Create a cyclic, 'loopable' waveform by cyclically extending the input,
        % filtering and then removing leading transient part
        ewaveform = [waveform(end-gnb.FilterLength:end,:); waveform];

        % Filter then realign
        efwaveform = filter(fnum, 1, ewaveform);
        waveform = efwaveform(end-size(waveform,1)+1:end,:);

        % Remove filter delay so that filtered OFDM symbols have the same
        % temporal alignment as other OFDM variants
        waveform = circshift(waveform,-halfFilt);    
    end
        
end