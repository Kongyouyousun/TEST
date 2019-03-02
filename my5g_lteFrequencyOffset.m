function [foffset, corr] = my5g_lteFrequencyOffset(cfg,waveform,txGrid,varargin)

    if (nargin==4)
        toffset=varargin{1};
    else
        toffset=[];
    end

    % Derive the number of samples per FFT, FFT duration and
    % the number of samples per slot.
    if (isfield(cfg,'NDLRB'))
        if (~isfield(cfg,'CyclicPrefix'))
            cfg.CyclicPrefix='Normal';
            lte.internal.defaultValueWarning('CyclicPrefix','Normal');
        end 
        info = h5gOFDMInfo(cfg);            
        link='Downlink';
    else
        if (~isfield(cfg,'CyclicPrefixUL'))
            cfg.CyclicPrefixUL='Normal';
           	lte.internal.defaultValueWarning('CyclicPrefixUL','Normal');
        end
        info=lteSCFDMAInfo(cfg);                
        cfg.CyclicPrefix=cfg.CyclicPrefixUL;
        link='Uplink';        
        halfsc=repmat(exp(1i*pi/double(info.Nfft)*(0:size(waveform,1)-1))',1,size(waveform,2));  
        waveform=waveform.*halfsc;
    end
    nFFT=double(info.Nfft);
    tFFT=1/(cfg.SubcarrierSpacing*1000);
    samplesPerSlot=info.SamplingRate*1e-3/info.SlotsPerSubframe;
    foffset=0;
    
    if (size(waveform,1)<(samplesPerSlot+nFFT))
        if isempty(waveform)
            foffset = [];
            corr = 0;
            return;
        end
        error('lte:error','The input waveform must span at least 1 slot plus the length of the FFT (%d+%d=%d samples at a sampling rate of %0.2fMs/s)',samplesPerSlot,nFFT,samplesPerSlot+nFFT,info.SamplingRate/1e6);
    end
    
    if (~isfield(cfg,'DuplexMode'))
        cfg.DuplexMode='FDD';
        lte.internal.defaultValueWarning('DuplexMode','FDD');
    end     
    
     % Initial symbol number of input wrt subframe
    initialsymbol = 0;  
    if isfield(cfg,'NSymbol')
        initialsymbol = mod(cfg.NSymbol,info.SymbolsPerSubframe);
    end
    
    % Waveform dependent initialization
    wola = isfield(cfg,'WaveformType') && strcmpi(cfg.WaveformType,'W-OFDM');
    % Derive the length of the cyclic prefixes for 1 slot at
    % the current sampling rate.
    cpLengths=double(info.CyclicPrefixLengths);
    if wola 
        
        % Perform W-OFDM specific initialization
        N1 = info.N1;               % Cyclic prefix part
        N2 = info.N2;               % Cyclic suffix part
        TL = nFFT + N1 + N2;        % Overall length of windowed W-OFDM symbol (not accounting for overlap - it's the window length) 
        
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
        pos = -(N1 + N2 - cpLengths(initialsymbol+1));
%       waveform = circshift(waveform,length(waveform)-N2,1);
        pos = pos+N2;
%          if isfield(cfg,'WindowCoeffs') && ~isempty(cfg.WindowCoeffs)
%            window = cfg.WindowCoeffs;
%            if length(cfg.window) ~= TL
%                error('lte:error','The W-OFDM window length must be equal to %d',TL);
%            end
%         else
%            window = sqrt(raised_cosine_window(nFFT,N1+N2));        
%         end 
    else
        exLengths = [cpLengths; zeros(size(cpLengths))];
        % Inter-symbol distances
        % Current stride length is the distance between the current symbol
        % and the start of the next one
        strideLengths = cpLengths + nFFT;
        pos =0;
    end
        
    % For each antenna:
    for i=1:size(waveform,2)
        p=pos;
        % Form two correlator inputs, the second delayed from
        % the first by nFFT.
        estimates = [];
        nSymbols = size(txGrid,2);
        for add=1:nSymbols
            currentsymbol = initialsymbol+add-1;           
            if wola
                cpLength = cpLengths(mod(currentsymbol+1,length(cpLengths))+1);
                stride = cpLength+nFFT;
                % Extract full extended symbol
                wdata = waveform(mod(p:p+TL-1,size(waveform,1))+1,i);
                % Internal N1 (CP) region
                arm2 = wdata(end-N2-N1+1:end-N2);
                arm1 = wdata(1:N1);
               cpcorrunfilt=arm1.*conj(arm2);
               estimates=[estimates;cpcorrunfilt];
                %Internal N2 (CS) region
                arm1= [arm1; wdata(N1+1:N1+N2)];
                arm2=[arm2; wdata(end-N2+1:end)];
            else
                % Add cyclic extension to the symbol
                exLength = cpLengths(mod(currentsymbol,length(exLengths))+1);   % Extension lengths (CP/CS) for current symbol
                stride = strideLengths(mod(currentsymbol,length(strideLengths))+1);
                arm1=waveform(p+(1:exLength),i);
                arm2=waveform(p + stride-exLength+1:p + stride,i);
            end
            
            % Conjugate multiply the inputs and integrate over the cyclic
            % prefix length.
            cpcorrunfilt=arm1.*conj(arm2);
            estimates=[estimates;cpcorrunfilt];
            p = p + stride;
        end
          
        % Average the estimates, take the angle and compute the
        % corresponding frequency offset.
        res =atan(imag(mean(estimates))/real(mean(estimates)));
        foffset=-res/(2*pi*tFFT);
    end
end
% Raised Cosine window creation; creates a window function of length n+w
% with raised cosine transitions on the first and last 'w' samples, where
% 'w' is the "number of time-domain samples (at each symbol edge) over 
% which windowing and overlapping of OFDM symbols is applied" as described
% in the product help
% function p = raised_cosine_window(n,w)
%     
%     p = 0.5*(1-sin(pi*(w+1-2*(1:w).')/(2*w)));
%     p = [p; ones(n-w,1); flipud(p)];
%     
% end

       

