received_signal = tx_out;
% %% Traverse the entire received signal looking for TPMS signal
% 
% %concatenate to 1-D array
% [x, y] = size(received_signal);
% rx = zeros(1, x*y);
% samples = 1:x*y;
% 
% for a = 0:x-1
%     rx(a*y+1:(a+1)*y) = received_signal(a+1,:);
% end
% %check BER, 
% power_threshold = 50*bandpower(rx(:,1));
% lower_ind = 1;
% threshold = 20*power_threshold;
% %lower_ind:lower_ind+10
% a = length(rx);
% b = size(bandpower(rx(lower_ind:lower_ind+10)));
% c = bandpower(rx(lower_ind:lower_ind+10));
% dataRx = rx(lower_ind:lower_ind+10);
% 
% while (bandpower(rx(lower_ind:lower_ind+10)) < threshold) && (lower_ind < length(rx))
%     lower_ind = lower_ind + 1; 
%     
%     if lower_ind >= 24991
%         disp('No signal found');    
%         break; 
%     end
% end
% 
% %locating the end point
% upper_ind = length(rx);
% while bandpower(rx(upper_ind-10:upper_ind)) < threshold && upper_ind > 1
%     upper_ind = upper_ind - 1; 
%     if upper_ind <= 10 disp('No signal found');
%         break;
%     end
% end
% 
% %holds the wave form of just the packet (no noise)
% packet_waveform = rx(lower_ind:upper_ind);
% %demodulate the packet
% figure
% plot(1:length(packet_waveform),packet_waveform)
% packet = demodulator( packet_waveform, Fs );
% %decode the packet
% %decode_packet(packet);
% %figure(2)
% %plot(1:length(packet_waveform), real(packet_waveform));
% 