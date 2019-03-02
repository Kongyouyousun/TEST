%ltePDSCHDecode Physical downlink shared channel decoding
%   [CWS,SYMBOLS] = ltePDSCHDecode(...) returns a cell array CWS of soft
%   bit vectors (codewords) and cell array SYMBOLS of received
%   constellation symbol vectors resulting from performing the inverse of
%   Physical Downlink Shared Channel (PDSCH) processing (TS 36.211 6.4, see
%   <a href="matlab:help('ltePDSCH')">ltePDSCH</a> for details). CWS is optionally scaled by channel state
%   information (CSI) calculated during the equalization process.
%
%   [CWS,SYMBOLS] = ltePDSCHDecode(ENB,CHS,SYM) performs the inverse of
%   Physical Downlink Shared Channel (PDSCH) processing on the matrix of
%   complex modulated PDSCH symbols SYM, cell-wide settings structure ENB
%   and channel-specific configuration structure CHS. The channel inverse
%   processing includes the deprecoding, layer demapping and codeword
%   separation, soft demodulation and descrambling. The deprecoding is
%   performed using matrix pseudo inversion of the precoding matrices.
%   
%   ENB must be a structure including the fields:
%   NCellID         - Physical layer cell identity
%   NSubframe       - Subframe number
%   CellRefP        - Number of cell-specific reference signal antenna  
%                     ports (1,2,4)
%   DuplexMode      - Optional. Duplex mode ('FDD'(default),'TDD')
%   Only required for 'TDD' duplex mode:
%      TDDConfig    - Optional. Uplink/Downlink Configuration (0...6) 
%                     (default 0)
%      SSC          - Optional. Special Subframe Configuration (0...9) 
%                     (default 0)
%   Only required for 'SpatialMux' and 'MultiUser' schemes (see CHS below):
%      NDLRB        - Number of downlink resource blocks  
%      CFI          - Control Format Indicator value (1,2,3)
%      CyclicPrefix - Optional. Cyclic prefix length 
%                     ('Normal'(default),'Extended')              
%
%   CHS must be a structure including the fields:
%   Modulation - A string array or cell array specifying the modulation
%                format for one or two codewords
%                ('QPSK','16QAM','64QAM','256QAM')
%   RNTI       - Radio Network Temporary Identifier (16-bit)
%   TxScheme   - Transmission scheme, one of:
%                'Port0'       - Single-antenna port, Port 0
%                'TxDiversity' - Transmit diversity scheme
%                'CDD'         - Large delay CDD scheme
%                'SpatialMux'  - Closed-loop spatial multiplexing scheme
%                'MultiUser'   - Multi-user MIMO scheme
%                'Port5'       - Single-antenna port, Port 5
%                'Port7-8'     - Single-antenna port, port 7 (when 
%                                NLayers=1); Dual layer transmission, port 
%                                7 and 8 (when NLayers=2)
%                'Port8'       - Single-antenna port, Port 8
%                'Port7-14'    - Up to 8 layer transmission, ports 7-14
%   NLayers    - Number of transmission layers (1...8)
%   CSI        - Optional. Determines if soft bits should be weighted by  
%                CSI ('Off'(default),'On')
%   Only required for 'SpatialMux' and 'MultiUser' schemes:
%      PMISet  - A vector of Precoder Matrix Indications. The vector may 
%                contain either a single value (corresponding to single PMI
%                mode) or multiple values (corresponding to multiple or
%                subband PMI mode). Values are in the range 0...15, with
%                the exact range depending on CellRefP, NLayers and
%                TxScheme. (See <a href="matlab:help('ltePMIInfo')"
%                >ltePMIInfo</a>)
%      PRBSet  - A 1- or 2-column matrix, containing the 0-based Physical 
%                Resource Block (PRB) indices corresponding to the slot
%                wise resource allocations for this PDSCH.
%   Only required for 'Port5', 'Port7-8', 'Port8' and 'Port7-14' schemes:
%      W       - NLayers-by-NTxAnts precoding matrix for the wideband
%                UE-specific beamforming of the PDSCH symbols.
%                   
%   For PRBSet, if a column vector is provided, the resource allocation is
%   the same in both slots of the subframe; the 2-column matrix can be used
%   to specify differing PRBs for each slot in a subframe. Note that the
%   PRB indices are 0-based.
%
%   SYM must be a matrix of NRE-by-NRxAnts complex modulated PDSCH symbols.
%   NRE is the number of QAM symbols per antenna assigned to the PDSCH and
%   NRxAnts is the number of receive antennas.
%
%   [CWS,SYMBOLS] = ltePDSCHDecode(ENB,CHS,SYM,HEST,NOISEEST) performs the
%   decoding of the complex modulated PDSCH symbols SYM using cell-wide
%   settings ENB, channel-specific configuration CHS, channel estimate HEST
%   and the noise estimate NOISEEST. For the TxDiversity transmission
%   scheme the deprecoding is performed using an OSFBC (Orthogonal Space
%   Frequency Block Code) decoder; for the SpatialMux, CDD and MultiUser
%   transmission schemes, the deprecoding is performed using a MIMO MMSE
%   equalizer, equalizing between transmitted and received layers. For the
%   Port0, Port5, Port7-8, Port8 and Port7-14 transmission schemes the
%   reception is performed using MMSE equalization; the input channel
%   estimate HEST is assumed to be with reference to the transmission
%   layers (using the UE-specific reference signals) so the MMSE
%   equalization will produce MMSE equalized layers.
%   
%   HEST is a 3-dimensional NRE-by-NRxAnts-by-CellRefP array for the Port0,
%   TxDiversity, SpatialMux, CDD and MultiUser transmission schemes, where
%   NRE is the number of QAM symbols per antenna assigned to the PDSCH,
%   NRxAnts is the number of receive antennas, and CellRefP is the number
%   of cell-specific reference signal antennas, given by ENB.CellRefP. For
%   the Port5, Port7-8, Port8 and Port7-14 transmission schemes HEST is of
%   size NRE-by-NRxAnts-by-NLayers, where NLayers is the number of
%   transmission layers given by CHS.NLayers.
%
%   NOISEEST is an estimate of the noise power spectral density per RE on
%   the received subframe; such an estimate is provided by the 
%   <a href="matlab: 
%   help('lteDLChannelEstimate')">lteDLChannelEstimate</a> function.
%
%   [CWS,SYMBOLS] = ltePDSCHDecode(ENB,CHS,RXGRID,HEST,NOISEEST) accepts
%   the full received resource grid RXGRID for one subframe, in place of
%   the SYM input; the decoder will internally extract the PDSCH REs to
%   obtain the complex modulated PDSCH symbols. RXGRID is a 3-dimensional
%   M-by-N-by-NRxAnts array of resource elements, where M and N are the
%   number of subcarriers and symbols for one subframe for cell-wide
%   settings ENB and NRxAnts is the number of receive antennas. In this
%   case, HEST is a 4-dimensional M-by-N-by-NRxAnts-by-CellRefP array where
%   M and N are the number of subcarriers and symbols for one subframe for
%   cell-wide settings ENB, NRxAnts is the number of receive antennas, and
%   CellRefP is the number of cell-specific reference signal antenna ports,
%   given by ENB.CellRefP. HEST will be processed to extract the channel
%   estimates relevant to the PDSCH i.e. those in the time/frequency
%   locations corresponding to the PDSCH REs in RXGRID.
%
%   Example:
%   % The complex PDSCH modulated symbols for RMC R.0 are generated and
%   % decoded for cell-wide setting enb.
%   
%   enb = lteRMCDL('R.0'); 
%   codewordBits = randi([0,1],enb.PDSCH.CodedTrBlkSizes(1),1);
%   pdschSym = ltePDSCH(enb,enb.PDSCH,codewordBits);
%   [rxCws,symbols] = ltePDSCHDecode(enb,enb.PDSCH,pdschSym);
%
%   See also ltePDSCH, ltePDSCHIndices, ltePDSCHPRBS, lteDLSCHDecode.

%   Copyright 2009-2018 The MathWorks, Inc.

function [cws,symCombining] = ltePDSCHDecode(varargin)

    % Obtain cell-wide settings enb and PDSCH configuration pdsch.
    enb=varargin{1};
    pdsch=varargin{2};
    
    % If empty input data then empty output
    if isempty(varargin{3})
            % Get the number of codewords
            [Ncodewords, ~] = getNcodewords(pdsch);
            % Create empty cells for the corresponding size of Ncodewords
            cws = cell(1,Ncodewords);
            symCombining = cell(1,Ncodewords);
        return;
    end
        
    % Remove NCodewords field if present; this is a parameter for
    % lteLayerDemap but not a parameter for this function.
    if (isfield(pdsch,'NCodewords'))
        pdsch = rmfield(pdsch,'NCodewords');
    end
    
    [Ncodewords, pdsch.Modulation] = getNcodewords(pdsch);

    % 3-argument signature.
    if (nargin==3)

        % Obtain received PDSCH symbols.
        pdschSymbols=varargin{3};
        
        appendNullSymbols=0;
        if (strcmpi(pdsch.TxScheme,'TxDiversity')==1 && pdsch.NLayers==4 && mod(size(pdschSymbols,1),4)~=0)
            appendNullSymbols=1;
            pdschSymbols=[pdschSymbols; zeros(2,4)];
        end
            
%         % Deprecoding (pseudo-inverse based).        
%         if (any(strcmpi(pdsch.TxScheme,{'Port5' 'Port7-8' 'Port8' 'Port7-14'})))        
%             if(~isfield(pdsch,'W'))
%                  error('lte:error','The function call (ltePDSCHDecode) resulted in an error: Could not find a structure field called W.');
%             end
%             if (~isempty(pdsch.W))
%                rxDeprecoded=pdschSymbols*pinv(pdsch.W);   
%             else
%                 % empty W, no deprecoding
%                 rxDeprecoded = pdschSymbols;
%             end
%         else
%             rxDeprecoded = lteDLDeprecode(enb,pdsch,pdschSymbols);
%         end

        % Layer demapping.
        symCombining = lteLayerDemap(pdsch,pdschSymbols);    
        
        if (appendNullSymbols)
            for i=1:Ncodewords
                temp=symCombining{i};
                symCombining{i}=temp(1:end-2);
            end
        end
        
        csi{1} = ones(size(symCombining{1}));
        if(Ncodewords==2)
            csi{2} = ones(size(symCombining{2}));
        end
    end

    % 5-argument signature.
    if (nargin==5)

        rxSubframe=varargin{3};
        chSubframe=varargin{4};         
        noiseEst=varargin{5};
        
        enb = mwltelibrary('validateLTEParameters',enb,'NDLRB');
        if (length(size(chSubframe))==4 || length(size(rxSubframe))==3 ||...
                (size(chSubframe,2)>8 && size(chSubframe,1)==(enb.NDLRB*12)) || (size(rxSubframe,2)>8 && size(rxSubframe,1)==(enb.NDLRB*12)))
            % Extract PDSCH from received subframe and extract corresponding channel estimates
            ind = ltePDSCHIndices(enb,pdsch,pdsch.PRBSet);
            MN=size(rxSubframe,1)*size(rxSubframe,2); % product of time and frequency dimensions
            NRE=size(ind,1); % number of channel REs in one antenna plane
            R=size(chSubframe,3); % number of receive antennas
            T=size(chSubframe,4); % number of transmit antennas
            rxInd = repmat(ind(:,1),[1 R])+repmat(uint32(0:R-1)*MN,[NRE 1]);
            chInd = repmat(rxInd,[1 1 T])+reshape(repmat(uint32(0:T-1)*MN*R,[NRE*R 1]),[NRE R T]);
            pdschSymbols=rxSubframe(rxInd);
            hest=chSubframe(chInd);            
        else            
            % assume PDSCH symbols are already extracted in the input
            pdschSymbols=rxSubframe;
            hest=chSubframe;
        end

        % perform deprecoding and layer demapping 
        if (strcmpi(pdsch.TxScheme,'TxDiversity')==1)

            appendNullSymbols=0;
            if (strcmpi(pdsch.TxScheme,'TxDiversity')==1 && pdsch.NLayers==4 && mod(size(pdschSymbols,1),4)~=0)
                appendNullSymbols=1;
                pdschSymbols=[pdschSymbols; zeros(2,size(hest,2))];
		        hest=[hest; zeros(2,size(hest,2),size(hest,3))];
            end

            % OSFBC decoding
            [symCombining,csi] = lteTransmitDiversityDecode(pdschSymbols,hest);
            symCombining={symCombining};
            csi={csi};

            if (appendNullSymbols)
                for i=1:Ncodewords
                    temp=symCombining{i};
                    symCombining{i}=temp(1:end-2);
                    csitemp=csi{i};
                    csi{i}=csitemp(1:end-2);
                end
            end

        else
            if (any(strcmpi(pdsch.TxScheme,{'CDD', 'SpatialMux', 'MultiUser'})))
                % MIMO equalization and deprecoding (MMSE based).
                [rxDeprecoded,csi] = lteEqualizeMIMO(enb,pdsch,pdschSymbols,hest,noiseEst);
            else           
                % MMSE for single antenna or UE-specific transmission
                % schemes (Port0, Port5, Port7-8, Port8, Port7-14).
                [pdschRx,csi] = lteEqualizeMMSE(pdschSymbols,hest,noiseEst);               
                rxDeprecoded = lteDLDeprecode(enb,pdsch,pdschRx);
            end

            % Layer demapping.
            symCombining = lteLayerDemap(pdsch,rxDeprecoded);
            csi = lteLayerDemap(pdsch,csi);
        end               
    end           

    % Provide default value for CSI field if absent
    if (~isfield(pdsch,'CSI'))
        pdsch.CSI='Off';
    end
    
    % Symbol demodulation and descrambling of codewords.
    if(Ncodewords == 1)                

        % Soft/hard demodulation of received symbols.
        demodpdschSymb = lteSymbolDemodulate(symCombining{1},pdsch.Modulation{1},'Soft'); 

        % Calculate CSI
        if ~isempty(symCombining{1})
            Qm=length(demodpdschSymb)/length(symCombining{1});
            csi{1}=repmat(csi{1}.',Qm,1);
            csi{1}=reshape(csi{1},numel(csi{1}),1);
        end
        
        % Scrambling sequence generation for descrambling.
        scramblingSeq = ltePDSCHPRBS(enb,pdsch.RNTI,0,length(demodpdschSymb),'signed');

        % Descrambling of received bits.
        cws = {demodpdschSymb.*scramblingSeq};       

        % Scaling LLRs by CSI if enabled
        if (strcmpi(pdsch.CSI,'On')==1)    
            cws{1}=cws{1}.*csi{1};
        end
        
    elseif(Ncodewords == 2)                

        % Soft/hard demodulation of received symbols.
        demodpdschSymbCW0 = lteSymbolDemodulate(symCombining{1},pdsch.Modulation{1},'Soft'); 
        demodpdschSymbCW1 = lteSymbolDemodulate(symCombining{2},pdsch.Modulation{2},'Soft'); 

        % Calculate CSI
        if ~isempty(symCombining{1})
            Qm1=length(demodpdschSymbCW0)/length(symCombining{1});
            csi{1}=repmat(csi{1}.',Qm1,1);
            csi{1}=reshape(csi{1},numel(csi{1}),1);
            Qm2=length(demodpdschSymbCW1)/length(symCombining{2});
            csi{2}=repmat(csi{2}.',Qm2,1);
            csi{2}=reshape(csi{2},numel(csi{2}),1);
        end
        
        % Scrambling sequence generation for descrambling.
        scramblingSeqCW0 = ltePDSCHPRBS(enb,pdsch.RNTI,0,length(demodpdschSymbCW0),'signed');
        scramblingSeqCW1 = ltePDSCHPRBS(enb,pdsch.RNTI,1,length(demodpdschSymbCW1),'signed');

        % Descrambling of received bits.
        CW0 = demodpdschSymbCW0.*scramblingSeqCW0;
        CW1 = demodpdschSymbCW1.*scramblingSeqCW1;

        % Scaling LLRs by CSI if enabled
        if (strcmpi(pdsch.CSI,'On')==1)    
            CW0 = CW0 .* csi{1};
            CW1 = CW1 .* csi{2};
        end
        
        % form final cell array of codewords
        cws = {CW0 CW1};

    end

end



% Calculate the number of codewords. The modulation scheme specified at the
% input in pdsch.Modulation is also returned as a cell array for the right
% number of codewords given the transmission scheme and the number of
% layers.
function [Ncodewords, Modulation] = getNcodewords(pdsch)

    Modulation = pdsch.Modulation;
    % Make modulation scheme a cell array
    if isstring(Modulation)
        Modulation = arrayfun(@char,Modulation,'UniformOutput',false);
    elseif ~iscell(Modulation)
        Modulation = {Modulation};
    end
    % Establish number of received codewords.
    if (strcmpi(pdsch.TxScheme,'TxDiversity')==1 || strcmpi(pdsch.TxScheme,'Port0')==1 || strcmpi(pdsch.TxScheme,'Port5')==1 || strcmpi(pdsch.TxScheme,'Port8')==1 || pdsch.NLayers==1)
        Ncodewords=1;
        Modulation=Modulation(1);
    else
        Ncodewords=size(Modulation,2);
    end

end
