function encoded_data = manchester_encode(binary_data)
encoded_data = zeros(1,2*length(binary_data));
    for i = 1:length(binary_data)
        if(binary_data(i) == 1)
            encoded_data(i*2-1:i*2) = [0 1];
        else
            encoded_data(i*2-1:i*2) = [1 0];
        end
    end
end

