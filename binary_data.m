% ----- returns N bits of binary data -----
function y = binary_data(N)
    
    y = randi([0 1],1,N); % uniformly generated binary baseband data
    
end