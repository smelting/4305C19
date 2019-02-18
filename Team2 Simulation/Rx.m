%% Rx
received_signal = tx_out;
% %% Traverse the entire received signal looking for TPMS signal
% 
%concatenate to 1-D array
[x, y] = size(received_signal);
rx = zeros(1, x*y);
samples = 1:x*y;
for a = 0:x-1
    rx(a*y+1:(a+1)*y) = received_signal(a+1,:);
end
symbol_sync_data = symbolSync(rx');
sync_data = carrSync(symbol_sync_data);

%% demodulate
demod_data = fskdemod(sync_data',M,freq_sep,nsamp,Fs);
%% find ID


start_ind = 1;
end_ind = 64;

ID_rx = demod_data(start_ind:end_ind);
%pressure, temp, ID, flags

pressure_rx = demod_data(start_ind - 16:start_ind-8);
temp_rx = demod_data(start_ind - 8:start_ind);
flags_rx = demod_data(end_ind:end_ind+8);
CRC_rx = demod_data(end_ind+8:end_ind+16);
encoded_rx = horzcat(pressure_rx, temp_rx, ID_rx, flags_rx);
%% decode
temp_data = zeros(1, length(encoded_rx));

%Manchester Decode
for i = 1:2:length(encoded_rx)-1
    if(encoded_rx(i:i+1) == [0 1])
        temp_data(i) = 1;
    else
        temp_data(i) = 0;
    end
end

decoded_data = downsample(temp_data, 2);

%figure out how to traverse signal, see if we don't need this
%assume we know ID, and yes need PLL



