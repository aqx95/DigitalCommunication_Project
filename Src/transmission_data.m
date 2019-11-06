% ----- Returns transmssion data (consisting of 1s and -1s) ------
function y = transmission_data(binary_data) %s0 is the baseband data
    y = 2.*binary_data - 1;  
end