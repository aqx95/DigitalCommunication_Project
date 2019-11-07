% ----- Apply Replication -----
function signal = replicate(orig_signal,value)
    signal = repelem(orig_signal,value);
end