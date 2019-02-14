clear all;close all;
%% General system details
sampleRateHz = 1e5; % Sample rate
samplesPerSymbol = 32;
frameSize = 68*samplesPerSymbol;
numFrames = 2;
numSamples = numFrames*frameSize; % Samples to simulate
modulationOrder = 2;
filterSymbolSpan = 8;
N = 2;
M = 2; %Modulation Order
Freq_Step = 1000;

frequencyOffsetHz = sampleRateHz*0.01; % Offset in hertz
phaseOffset = 0; % Radians

%% Impairments
snr = 25;
timingOffset = samplesPerSymbol*0.01; % Samples

preamble = [ 1; 1; 0; 1;];
pressure = [1; 1; 1; 1; 1; 1; 1; 1;];
temp = [1; 1; 1; 1; 0; 0; 0; 0;];
tpmsID = [ 1; 1; 1; 1; 1; 1; 1; 1; 0; 0; 0; 0; 0; 0; 0; 0;...
     1; 1; 1; 1; 1; 1; 1; 1; 0; 0; 0; 0; 0; 0; 0; 0;];
flags = [ 1; 1; 0; 0; 1; 1; 0; 0;];
CRC = [ 0; 0; 1; 0; 1; 0; 1; 1; ];

%% Generate symbols
data = [preamble; pressure; temp; tpmsID; flags; CRC;];
modulatedData = fskmod(data,M,Freq_Step,samplesPerSymbol,sampleRateHz);

%% Add TX/RX Filters
TxFlt = comm.RaisedCosineTransmitFilter(...
    'OutputSamplesPerSymbol', samplesPerSymbol,...
    'FilterSpanInSymbols', filterSymbolSpan);

RxFlt = comm.RaisedCosineReceiveFilter(...
    'InputSamplesPerSymbol', samplesPerSymbol,...
    'FilterSpanInSymbols', filterSymbolSpan,...
    'DecimationFactor', samplesPerSymbol/2);% Set to filterUpsample/2 when introducing timing estimation

%% Add noise source
chan = comm.AWGNChannel( ...
    'NoiseMethod',  'Signal to noise ratio (SNR)', ...
    'SNR',          snr, ...
    'SignalPower',  20, ...
    'RandomStream', 'mt19937ar with seed');

%% Add delay
varDelay = dsp.VariableFractionalDelay;

%% Setup visualization object(s)
sa = dsp.SpectrumAnalyzer('SampleRate',sampleRateHz,'ShowLegend',true);
ds = dsp.TimeScope('SampleRate',sampleRateHz,'ShowLegend',true,'TimeSpan',0.025);


% Precalculate constants
normalizedOffset = 1i.*2*pi*frequencyOffsetHz./sampleRateHz;

offsetData = zeros(size(modulatedData));
offsetDataFreq = zeros(size(modulatedData));
for k=1:frameSize:(numSamples - frameSize)
    timeIndex = (k:k+frameSize-1).';
    % Filter signal
    filteredTXData = step(TxFlt, modulatedData(timeIndex));
    
    
    % Pass through channel
    noisyData = step(chan, filteredTXData);
    % Time delay signal
    offsetData = step(varDelay, noisyData, k/frameSize*timingOffset); % Variable delay
    
    %step(cdPre,noisyData(timeIndex));%step(cdPost,offsetDataFreq(timeIndex));pause(0.1); %#ok<*UNRCH>
    
end

for k=1:frameSize:numSamples*samplesPerSymbol/2
    
    timeIndex = (k:k+frameSize-1).';
    freqShift = exp(normalizedOffset*timeIndex + phaseOffset);
    
    % Offset data and maintain phase between frames
    offsetDataFreq(timeIndex) = noisyData(timeIndex).*freqShift;
    filteredData = step(RxFlt, offsetDataFreq);
end

symsync = comm.SymbolSynchronizer( ... 
'SamplesPerSymbol', 2, ...
'DampingFactor', sqrt(2)/2, ...
'NormalizedLoopBandwidth', 0.01);
%rxSym = symsync(filteredData);
rxSym = downsample(filteredData,2);
sa(rxSym);
ds(rxSym);
pause;
sa.release();
lenRx = length(rxSym);
t = (1:lenRx)/sampleRateHz;
freq = max_frequencies(rxSym, sampleRateHz,2);
offset = -(freq(1) + freq(2))/2;
mod_sig = exp(1i.*2*pi*offset.*t).';
rxSym = rxSym .* mod_sig;
sa(rxSym);
receivedData = fskdemod(rxSym,M,Freq_Step,samplesPerSymbol,sampleRateHz);
pause;
figure(1);
stem(data,'red');
title('Transmitted Data');
figure(2);
stem(receivedData,'blue');
title('Received Data');
[num,BER] = biterr(data,receivedData);
BER