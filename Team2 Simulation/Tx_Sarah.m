%% Tx
% time = 0:0.0002647:2;
time = 0:2.649006623e-4:2;
t = vertcat(time,time,time);
frameSize = 2^10;
numFrames = 200; %200
numSamples = numFrames*frameSize; % Samples to simulate
cdPost = comm.ConstellationDiagram('SymbolsToDisplaySource','Property',...
    'SymbolsToDisplay',frameSize/2,...
    'Name','Baseband with Timing Offset');
% keep random data, 255 is the limit
pressure_dec = 170;
temp_dec = 50;
ID_hex = ('FF103ABC');

%ID (hex values 32 bits in total)
ID = hex2bin(ID_hex);

%tire pressure (8 bits)
pressure = str2bin(dec2bin(pressure_dec, 8));

%tire temperature (8 bits)
temp = str2bin(dec2bin(temp_dec, 8));

%flags (8 bits)
flags = randi([0 1],1,8);

%concatenate data 
concat_data = horzcat(pressure, temp, ID, flags);

%add CRC-8
CRC = comm.CRCGenerator('z^8 + z^7 + z^6 + z^4 + z^2 + 1');
CRC_data = CRC(concat_data');
raw_data = CRC_data';

%manchester encode
encoded_data = manchester_encode(raw_data);

%preamble
preamble = randi([0 1],1,23);

%data to be modulated
full_data = horzcat(preamble, encoded_data);

%ASK Modulation
Fc = 4.33e8;
Fs_ask = 2*Fc;
%ask_data = askMod(full_data,2*Fc);
snr = 25;

% stem(full_data);
% hold on
% plot(ask_data)
% plot(abs(fft(ask_data)));
% plot(ask_data);

%% ASK Modulation

x= full_data;                                    % Binary Information -> modulatedData
bp=(Fs_ask/2)*log2(1+snr); %shannon theorem                                                  % bit period 

disp(' Binary information at Trans mitter :');
disp(x);
%XX representation of transmitting binary information as digital signal XXX
bit=[]; 
for n=1:1:length(x)
    if x(n)==1;
       se=ones(1,100);
    else x(n)==0;
       se=zeros(1,100);
    end
     bit=[bit se];
end
t1=bp/100:bp/100:100*length(x)*(bp/100);
subplot(3,1,1);
plot(t1,bit,'lineWidth',2.5);grid on;
axis([ 0 bp*length(x) -.5 1.5]);
ylabel('amplitude(volt)');
xlabel(' time(sec)');
title('transmitting information as digital signal');
%XXXXXXXXXXXXXXXXXXXXXXX Binary-ASK modulation XXXXXXXXXXXXXXXXXXXXXXXXXXX%
A1=10;                      % Amplitude of carrier signal for information 1
A2=5;                       % Amplitude of carrier signal for information 0
br=1/bp;                                                         % bit rate
f=br*10;                                                 % carrier frequency 
t2=bp/99:bp/99:bp;                 
ss=length(t2);
m=[];
for (i=1:1:length(x))
    if (x(i)==1)
        y=A1*cos(2*pi*f*t2);
    else
        y=A2*cos(2*pi*f*t2);
    end
    m=[m y];
end
t3=bp/99:bp/99:bp*length(x);
subplot(3,1,2);
plot(t3,m);
xlabel('time(sec)');
ylabel('amplitude(volt)');
title('waveform for binary ASK modulation coresponding binary information');

% cdPost(m);
% plot(ask,'b','linewidth',1.5)
% title('ASK Modulation');
% grid on
plot(abs(fft(m, 1024)));
data_out_ask = vertcat(m,m,m);

%FSK data
M = 2;
freq_sep = 73730; %peaks were that far apart maybe ask if we can still use
nsamp = 50; %check, number of samples per symbol
Fs = 100e6/128; %check
fsk_data = fskmod(full_data,M,freq_sep,nsamp,Fs);
plot(abs(fft(fsk_data)))
data_out = vertcat(fsk_data,fsk_data,fsk_data);

%% channel 
%add noise
phase_offset = comm.PhaseFrequencyOffset('PhaseOffset', 50, 'SampleRate', Fs,'FrequencyOffset', 143857490);
timingOffset = 50*0.01; % Samples
phase_data = phase_offset(data_out_ask);
tx_out = awgn(phase_data, snr);
% f = 10;
% phase_offset = exp(i*2*pi*f*t + 3*pi/4);
% noise = 25*randn(3, 7550);
% phase_data = data_out .* phase_offset;
% tx_out = phase_data + noise;

%% Fix
carrSync = comm.CarrierSynchronizer('SamplesPerSymbol',50);
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

%symbol_sync_data = symbolSync(rx');
sync_data = carrSync(rx');
% step(cdPost, sync_data);
% plot(abs(fft(sync_data)));
cdPost(sync_data)
%% traverse signal? possible might be before PLL need to ask
 y = find_mod_scheme(sync_data);
 demod_data = ASK_demodulator(data_out_ask);
%% demodulate

% if y == 1
%     demod_data = ~fskdemod(sync_data,M,freq_sep,nsamp,Fs)';
% else
%     demod_data = ASK_demodulator(sync_data);
% end
% need to fix modulations, channel stuff, PLLs, demodulations
%% find ID
ID_given = ('FF103ABC');
ID_bin = hex2bin(ID_given);
ID_encoded = manchester_encode(ID_bin);
start_ind = find_ID(demod_data, ID_encoded);

%if ID is too close to the beginning or if finds ID in multiple places
%chooses 1
y =1;

if(length(start_ind) > 2)
    start_ind = start_ind(2);
end

end_ind = start_ind + 63;

ID_rx = demod_data(start_ind:end_ind);

%pressure, temp, ID, flags
temp_rx = demod_data(start_ind - 16:start_ind-1);

pressure_rx = demod_data(start_ind - 32:start_ind-17);

flags_rx = demod_data(end_ind+1:end_ind+16);

CRC_rx = demod_data(end_ind+17:end_ind+32);

%% Manchester decode

%convert back to decimal
pressure_value = bi2de(manchester_decode(pressure_rx),'left-msb')
temp_value = bi2de(manchester_decode(temp_rx),'left-msb')

%do you want a timing offset and phase offset, or a phase offset and a freq
%offset

%how to prove that it works?

%what should be included in the vid?