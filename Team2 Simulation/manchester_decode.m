function decoded_data = manchester_decode(binary_data)
temp_data = zeros(1, length(binary_data));
%Manchester Decode
    for i = 1:2:length(binary_data)-1
        if(binary_data(i:i+1) == [0 1])
            temp_data(i) = 1;
        else
            temp_data(i) = 0;
        end
    end
    decoded_data = downsample(temp_data, 2);
end

