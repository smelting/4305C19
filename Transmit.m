tx = sdrtx('Pluto','Gain',-1, 'CenterFrequency', 433e6);

cd = comm.ConstellationDiagram; % Constellation Diagram setup 

%% General system details
Fc = 4.33e8;
Fs = 32;
t = 1/Fs;
frameSize = 1096;
numFrames = 200;
numSamples = numFrames*frameSize; 
modulationOrder = 2;
filterSymbolSpan = 4;
samplesPerFrame = numSamples/numFrames;
decimation = 8;
%% Preamble Bits
preamble = [1 1 1 0 0 0 0 0 0];
%% Temperature Bits
temp = [0 1 0 1 1 0 1 0];
 
rand_temp = randi([30 99],1);
temp1 = dec2bin(rand_temp);
%% Pressure  Bits
pres = [1 1 1 1 0 0 0 0];
 
rand_pres = randi([150 330],1);
pres1 = dec2bin(rand_pres);
%% Flag Bits
flag = [1 1 1 1 1 0 0 0];
%% ID Bits
ID = [1 0 0 0 1 0 0 0 0 1 1 1 1 1 0 0 0 1 1 0 1 0 0 1 1 1 1 1 1 0 0 1];
 
symbols = ['A':'F' '0':'9'];
stLength = 8;
 
nums = randi(numel(symbols),[1 stLength]);
ID_hex = symbols (nums);
ID_bin = Hex_to_Bin(ID_hex);
%% CRC
crc =[0 1 0 1 1 1 0 0];
%% Encode data 
data = horzcat(ID, temp, pres, flag,crc);
data_encoded = man_encode(data);
data_encoded_preamble = horzcat(preamble, data_encoded);
%% Modulate
M = 4;        % Modulation order
freqsep = 8;  % Frequency separation (Hz)
sps = 8;    % Number of samples per symbol
     % Sample rate (Hz)
modulatedData = transpose(fskmod(data_encoded_preamble,M,freqsep,sps,Fs));
%% TX Filter
b = .1;
 
rrcFilter = rcosdesign(b, filterSymbolSpan, sps, 'sqrt');
txSig = upfirdn(modulatedData, rrcFilter, sps);
add = [0];
txSig = vertcat(txSig,add);

tx.transmitRepeat(txSig);
cd(txSig); 

%% Manchester Encoding from MQP
function [encoded_signal] = man_encode( signal )
%man_encode
% this function takes in a packet and manchester encodes it
%encoded signal is twice the length
encoded_signal = zeros(1, 2*length(signal));
for x = 1:2:length(encoded_signal)
 if signal((x+1)/2) == 1
 encoded_signal(x) = 1;
 else
 encoded_signal(x) = 0;
 end
 encoded_signal(x+1) = xor(encoded_signal(x), 1);
end
end
%% Hex to Binary Function from MQP 
function [ bin ] = Hex_to_Bin( ID_hex )
%Hex_to_Bin
% this function takes in the ID in a hex string and converts it to a
% binary array
binStr = dec2bin(hex2dec(strcat('A', ID_hex)));
bin = zeros(1, length(binStr));
for x = 1:length(bin)
 bin(x)= str2double(binStr(x));
end
 bin = bin(5:end);
end