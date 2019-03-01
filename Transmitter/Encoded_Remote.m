function [transmit] = Encoded_Remote(data_in,upsample,off_time)
%ENCODED_REMOTE Summary of this function goes here
%   Detailed explanation goes here

off_time = zeros(off_time,1);
encoded_data = [];

%%encode the data based on baseband measurements from remote
for i = 1:1:length(data_in)
    data_in(i)
    if data_in(i) == 0
       encoded_data = [encoded_data; 1; 0; 0; 0;]; 
    elseif data_in(i) == 1
       encoded_data = [encoded_data; 1; 1; 1; 0;];
    end
end


%%Upsample binary stream
upsampled = [];
for i = 1:1:length(encoded_data)
    for x = 1:1:upsample
        if encoded_data(i) == 1
           upsampled = [upsampled; 1]; 
        end
        if encoded_data(i) == 0
            upsampled = [upsampled; 0];
        end
    end
end

%%add some off time to the transmitted signal
upsampled = [upsampled; off_time;]; 

%convert the encoded data into complex form so the pluto can transmit it
transmit = complex(upsampled,zeros(length(upsampled),1));
end

