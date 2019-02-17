% Setup Receiver 
rx=sdrrx('Pluto','OutputDataType','double','SamplesPerFrame',2^15,...
            'CenterFrequency', 433e6);
 
       
 % Scope Setup 
samplesPerStep = rx.SamplesPerFrame/rx.BasebandSampleRate;
steps = 3;
ts = dsp.TimeScope('SampleRate', rx.BasebandSampleRate,...
'TimeSpan', samplesPerStep*steps,...
'BufferLength', rx.SamplesPerFrame*steps);  % Time Scope Setup
 
cd = comm.ConstellationDiagram; % Constellation Diagram setup 
rxSig= rx();
rxSig = rxSig(1:8793);
 
scatterplot(rxSig); 
title('received')
 
%% General system details
sampleRateHz = 1e6; % Sample rate
samplesPerSymbol = 4;
frameSize = 2^10;
numFrames = 200;
numSamples = numFrames*frameSize; % Samples to simulate
modulationOrder = 2;
filterSymbolSpan = 4;
samplesPerFrame = numSamples/numFrames;
t = 1/sampleRateHz;
decimation = 4;
M = 4;        % Modulation order
freqsep = 8;  % Frequency separation (Hz)
sps = 8;    % Number of samples per symbol
     % Sample rate (Hz)
Fs = 32;
% % Setup PLL
% S = comm.SymbolSynchronizer('Modulation', 'OQPSK',...
%                             'TimingErrorDetector', "Zero-Crossing (decision-directed)",...
%                             'SamplesPerSymbol', samplesPerSymbol,...
%                             'DampingFactor', 1);
 
%% Impairments
snr = 15;
timingOffset = samplesPerSymbol*0.01; % Samples
%% Add delay
varDelay = dsp.VariableFractionalDelay;
 
%% Add RX Filters
b = 1;
 
% Rx Filter
rrcFilter = rcosdesign(b, filterSymbolSpan, sps, 'sqrt');
rxFilt = upfirdn(rxSig, rrcFilter, 1, sps);
rxFilt = rxFilt(filterSymbolSpan+1:end-filterSymbolSpan);
 
    
 
 
scatterplot(rxFilt); 
title('filtered');
 
%% Demodulate
demodulatedData = fskdemod(rxFilt, M, freqsep,sps,Fs);
 
%% Cutoff Preamble
demodulated1 = transpose(demodulatedData(10:137));
%demodulated2 = transpose(floor(correctedData(10:137)));
%demodulated1 = data_encoded; 
%% Decode
decodedData = man_decode(demodulated1);
%decodedData = man_decode(demodulated2);
%% Convert Data to readable form
preamble_recovered = demodulatedData(1:9);
 
ID_decoded = decodedData(1:length(ID));
ID_recovered = binaryVectorToHex(ID_decoded);
 
temp_decoded = decodedData(33:40);
temp_recovered = bi2de(temp_decoded);
 
pres_decoded = decodedData(41:48);
pres_recovered = bi2de(pres_decoded,'left-msb');
 
flag_decoded = decodedData(49:56);
flag_recovered = binaryVectorToHex(flag_decoded);
 
crc_decoded =  decodedData(57:end);
crc_recovered = binaryVectorToHex(crc_decoded);
 
 
% Manchester Decoding from MQP
function [ decoded ] = man_decode( encoded )
%man_decode
%This function takes in a manchester encoded signal and decodes it
decoded = zeros(1, length(encoded)/2);
%takes every other value starting at 1
for x = 1:2:length(encoded)
 %making sure two succesive bits are not the same
 if encoded(x) ~= encoded(x+1)
 decoded((x+1)/2) = encoded(x);
 else
 decoded = -1;
 break;
 end
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