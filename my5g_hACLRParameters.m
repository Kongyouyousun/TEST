function [aclr, nRC, R_C, BWUTRA] = my5g_hACLRParameters(cfg)

    % Calculate the channel and transmission bandwidths. The number of
    % resource blocks (NRB) is established from cfg.NDLRB/cfg.NULRB.
    if (isfield(cfg, 'NDLRB'))
        NRB = cfg.NDLRB;
    else
        NRB = cfg.NULRB;
    end

    % Channel bandwidth of the input signal, calculated from NRB
    aclr.Bandwidth = hNRBToBandwidth(NRB)*1e6;    

    % Transmission bandwidth configuration, calculated from NRB,
    % subcarriers per RB and subcarrier bandwidth
    aclr.BandwidthConfig = double(NRB)*saLteNscRB()*15e3;

    % Find the overall signal bandwidth to support 2nd adjacent E-UTRA
    % signal
    requiredBWEUTRA = ((2*aclr.Bandwidth) + (aclr.Bandwidth/2))*2;

    % Find the overall signal bandwidth to support 2nd adjacent UTRA
    % signal. Establish default UTRA chip rate if absent
    if (~isfield(cfg,'UTRAChipRate'))
        cfg.UTRAChipRate = 3.84;        
        s = warning('query','backtrace');
        warning off backtrace;
        warning('lte:defaultValue',['Using default value for parameter' ...
                'field UTRAChipRate (%0.2f)'],cfg.UTRAChipRate);
        warning(s);                 
    end

    % Establish validity of UTRA chip rate parameters and establish
    % corresponding UTRA chip rates and UTRA channel bandwidths
    nRC = length(cfg.UTRAChipRate);
    R_C = zeros(1, nRC);
    BWUTRA = zeros(1, nRC);
    for k = 1:nRC
        [R_C(k),BWUTRA(k)] = hUTRAParameters(cfg.UTRAChipRate(k));
    end            

    % The overall signal bandwidth is calculated using the widest
    % configured UTRA channel bandwidth
    maxBWUTRA = max(BWUTRA);
    requiredBWUTRA = ((maxBWUTRA*3/2 + (aclr.Bandwidth/2)) + ...
        (maxBWUTRA/2))*2;

    % The ACLR bandwidth is the larger of the required E-UTRA and UTRA
    % bandwidths
    aclr.BandwidthACLR = max(requiredBWEUTRA,requiredBWUTRA);
    
    aclr.Bandwidth=aclr.Bandwidth*cfg.SubcarrierSpacing/15;
    aclr.BandwidthConfig = aclr.BandwidthConfig*cfg.SubcarrierSpacing/15;
    aclr.BandwidthACLR = aclr.BandwidthACLR*cfg.SubcarrierSpacing/15;

    % Calculate the measurement oversampling ratio and sampling rate,
    % allowing for at most 85% bandwidth occupancy
    aclr.OSR = ceil((aclr.BandwidthACLR/0.85)/double(cfg.SamplingRate));        
    aclr.SamplingRate = double(cfg.SamplingRate)*aclr.OSR;   

end