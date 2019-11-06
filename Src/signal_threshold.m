% ----- Apply threshold -----
function y = signal_threshold(s0, threshold)
        for cols = 1:size(s0,2)
            if s0(cols) >= threshold
                s0(cols) = 1;
            else
                s0(cols) = 0;
            end
        end
    y = s0;
end