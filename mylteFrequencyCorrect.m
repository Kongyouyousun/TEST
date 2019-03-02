function [out] = mylteFrequencyCorrect(in,foffset)
    info=struct();
    info.SamplingRate=122.88e6;
 
    % Create output storage.
    out=complex(zeros(size(in)));
    
    % Create vector of time samples.
    t=((0:size(in,1)-1)/info.SamplingRate).';   
    
    % For each antenna, apply the frequency offset correction.
    for i=1:size(in,2)
        out(:,i) = in(:,i).*exp(-1i*2*pi*foffset*t);
    end
    
end