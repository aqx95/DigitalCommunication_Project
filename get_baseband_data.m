% ----- returns N bits of baseband data -----
function y = get_baseband_data(N)
    
    y = randi([0 1],1,N); % uniformly generated binary baseband data
    
end