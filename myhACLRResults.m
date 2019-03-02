function myhACLRResults(aclr)

    disp(aclr);

    % Plot E-UTRA ACLR
    values = [-aclr.EUTRAdB(1:end/2) 0 -aclr.EUTRAdB(end/2+1:end)];
    tick = 1:numel(values);
    ticklabel = tick-ceil(numel(tick)/2);
    labelvec = tick;
    labelvec(ceil(end/2)) = []; % Do not plot label for 0dB ACLR on channel
    
    figure;
    bar(values, 'BaseValue', -90, 'FaceColor', 'yellow');
    set(gca, 'XTick', tick, 'XTickLabel', ticklabel, 'YLim', [-90 0]);
    for i = labelvec
        text(i, values(i), sprintf('%0.2f dB',values(i)), ...
            'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Top');
    end
    title('E-UTRA Adjacent Channel Leakage Ratio');
    xlabel('Adjacent Channel Offset');
    ylabel('Adjacent Channel Leakage Ratio (dB)');
    
end