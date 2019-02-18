function binary = hex2bin(ID)

    binary_string = dec2bin(hex2dec(strcat('A', ID))); 
    binary = zeros(1, length(binary_string));
    
    for x = 1:length(binary)
        binary(x)= str2double(binary_string(x));
    end
    binary = binary(5:end);
    
end

