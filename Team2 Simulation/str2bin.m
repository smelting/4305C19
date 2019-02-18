function binary = str2bin(bin_string)
 binary = zeros(1, length(bin_string));
    for x = 1:length(binary)
        binary(x)= str2double(bin_string(x));
    end
end

