function idx = util_chanidx(EEG, labels)
% Return indices of requested channel labels present in EEG (case-insensitive).
% Warns about any label not found so ROI typos surface immediately.
have = {EEG.chanlocs.labels};
idx = [];
for i = 1:numel(labels)
    j = find(strcmpi(have, labels{i}), 1);
    if isempty(j)
        warning('channel "%s" not found in montage', labels{i});
    else
        idx(end+1) = j; %#ok<AGROW>
    end
end
end
