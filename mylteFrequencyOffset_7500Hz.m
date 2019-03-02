function foffset = mylteFrequencyOffset_7500Hz(waveform,varargin)

    if (nargin==3)
        toffset=varargin{1};
    else
        toffset=[];
    end
 
   % Derive the number of samples per FFT, FFT duration and
    % the number of samples per slot.
%     if (isfield(cfg,'NDLRB'))
%         if (~isfield(cfg,'CyclicPrefix'))
%             cfg.CyclicPrefix='Normal';
%             defaultValueWarning('CyclicPrefix','Normal');
%         end 
        %info=lteOFDMInfo(cfg); 
        info=struct();
        info.SamplingRate=122.88e6;
        info.Nfft=2048*8;
        info.CyclicPrefixLengths=4096*ones(1,12);
        info.SubframeType='Downlink';
        link='Downlink';
%     else
%         if (~isfield(cfg,'CyclicPrefixUL'))
%             cfg.CyclicPrefixUL='Normal';
%            	defaultValueWarning('CyclicPrefixUL','Normal');
%         end
%         %info=lteSCFDMAInfo(cfg);
%         info=struct();
%         info.SamplingRate=prmLTEPDSCH.chanSRate;
%         info.Nfft=prmLTEPDSCH.N;
%         info.CyclicPrefixLengths=[160 144 144 144 144 144 144 160 144 144 144 144 144 144];
%         cfg.CyclicPrefix=cfg.CyclicPrefixUL;
%         info.SubframeType='Uplink';
%         link='Uplink';        
%         halfsc=repmat(exp(1i*pi/double(info.Nfft)*(0:size(waveform,1)-1))',1,size(waveform,2));  
%         waveform=waveform.*halfsc;
%     end
    nFFT=double(info.Nfft);
    tFFT=1/7500;
    samplesPerSlot=info.SamplingRate*1e-3;
    foffset=0;
    
    if (size(waveform,1)<(samplesPerSlot+nFFT))
        if isempty(waveform)
            foffset = [];
            %corr=0;
            return;
        end
        error('lte:error','The input waveform must span at least 1 slot plus the length of the FFT (%d+%d=%d samples at a sampling rate of %0.2fMs/s)',samplesPerSlot,nFFT,samplesPerSlot+nFFT,info.SamplingRate/1e6);
    end
    
%     if (~isfield(cfg,'DuplexMode'))
%         cfg.DuplexMode='FDD';
%         defaultValueWarning('DuplexMode','FDD');
%     end      
    
    % Derive the length of the cyclic prefixes for 1 slot at
    % the current sampling rate.
    cpLengths=double(info.CyclicPrefixLengths);
    cpLengths=cpLengths(1:length(cpLengths)/2);
    cpLength=cpLengths(2);
    
    % Compute the number of samples through a slot where
    % each of the OFDM symbols start.
    symbolStarts=cumsum(cpLengths+nFFT);
    symbolStarts0=[0 symbolStarts(1:end-1)];    
        
    % For each antenna:
    bestCorr=-1;    
    for i=1:size(waveform,2)
        
        % Form two correlator inputs, the second delayed from
        % the first by nFFT.
        arm1=waveform(1:end-nFFT,i);
        arm2=waveform(1+nFFT:end,i);

        % Conjugate multiply the inputs and integrate over the cyclic
        % prefix length.
        cpcorrunfilt=arm1.*conj(arm2);
        cpcorr=conv(cpcorrunfilt,ones(cpLength,1));
        cpcorr=cpcorr(cpLength:end);  
        
        % Average the estimates by combining all available slots.
        cpcorravg=cpcorr(1:(fix(length(cpcorr)/samplesPerSlot)*samplesPerSlot)); 
        nSlots=length(cpcorravg)/(samplesPerSlot);
        cpcorravg=sum(reshape(cpcorravg,samplesPerSlot,length(cpcorravg)/samplesPerSlot),2);  
        
        % Take the absolute value of the averaged output.
        cpcorrmag=abs(cpcorravg);   
       % corr(:,i)=cpcorrmag; %#ok<AGROW>
        
        % If the peak of the correlation on this antenna is higher
        % than for any antenna considered so far:
        if (max(cpcorrmag)>bestCorr)
            
            % update the peak value.
            bestCorr=max(cpcorrmag);
            
            % Extract the timing of the peak correlation relative to 
            % the start of the nearest OFDM symbol in the slot.
            idx=1:1:length(cpcorrmag);
            if (isempty(toffset))                
                cycshift=idx(cpcorrmag==max(cpcorrmag));   
                cycshift=cycshift(1);
            else
                cycshift=toffset;
            end   
            tail=samplesPerSlot-symbolStarts0(end);
            cycshift=mod(cycshift+tail,samplesPerSlot)-tail;
            candidates=-symbolStarts0+cycshift;
            candidateidxs=1:length(symbolStarts0);
            pos=candidateidxs(abs(candidates)==min(abs(candidates)));
            cycshift=min(candidates(pos));
            
            % Form a vector of the locations of all OFDM symbols in the
            % original waveform.
            freqIdxs=[];
            for l=1:nSlots                           
                %info=lteDuplexingInfo(cfg);
                if (strcmp(info.SubframeType,link)==1)
                    freqIdxs=[freqIdxs samplesPerSlot*(l-1)+symbolStarts]; %#ok<AGROW>
                else
                    if (~isempty(freqIdxs))
                        freqIdxs(end)=[];
                    end
                end
%                 if (mod(l,2)==0 && strcmp(cfg.DuplexMode,'TDD')==1)
%                     cfg.NSubframe=mod(cfg.NSubframe+1,10);
%                 end
            end
             if (~isempty(freqIdxs)) 
                freqIdxs(end)=[];

                % Form a vector of all samples of the correlator input
                % within 1/2 the cyclic prefix length.
                estimates=[];
                for add=fix(cpLength/2):fix(cpLength)-1
                    estimates=[estimates; cpcorrunfilt(cycshift+freqIdxs+add)]; %#ok<AGROW>
                end           

                % Average the estimates, take the angle and compute the
                % corresponding frequency offset.
                foffset=-angle(mean(estimates))/(2*pi*tFFT);  
             end
        end
           
    end
        
end

function defaultValueWarning(field,value)
    s=warning('query','backtrace');
    warning off backtrace;        
    warning('lte:defaultValue','Using default value for parameter field %s (%s)',field,value);
    warning(s); 
end
