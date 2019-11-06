% -----Compute Bit error-----
function error_rate = compute_error(input,received)
    error = 0;
    for k = 1: size(input,2)
        if input(k) ~= received(k)
            error = error + 1;      
        end
    end
    error_rate = error / size(input,2);       %calculate bit error rate
end