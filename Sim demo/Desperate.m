%% Tx system specs
sampleRateHz = 1e6; % Sample rate 
%100e6/128 apparently

sPS = 4; %samples per symbols
frameSize = 20;
numFrames = 100;
numSamples = numFrames*frameSize; % Samples to simulate
modulationOrder = 2; %only 4 and 2
filterSymbolSpan = 4;
encoderState=0;
encodedData = [];


%% Data setup
data = randi([0 1], numSamples, 1);

%% Visualizaton 
cdPLL = comm.ConstellationDiagram('ReferenceConstellation', [-1 1],...
    'SymbolsToDisplaySource','Property',...
    'SymbolsToDisplay',frameSize/2,...
    'Name','Baseband with Timing Offset after PLL');
%% Encoder

if modulationOrder== 4
    tempData = [];
   for i= 1:2:length(data)
       tempData = [data(i);data(i+1)];
       encodedData = [encodedData, tempData];
   end
end
encodedData = encodedData';
%% Modulation
if modulationOrder == 2
    mod = comm.BPSKModulator();
    modulatedData = mod.step(data);
    demod = comm.BPSKDemodulator();
else
    mod = comm.QPSKModulator('BitInput', true);
    modulatedData = mod.step(data);
    demod = comm.QPSKDemodulator();
end
%% Matched Filter Setup
TxFlt = comm.RaisedCosineTransmitFilter(...
    'OutputSamplesPerSymbol', sPS,...
    'FilterSpanInSymbols', filterSymbolSpan);

RxFlt = comm.RaisedCosineReceiveFilter(...
    'InputSamplesPerSymbol', sPS,...
    'FilterSpanInSymbols', filterSymbolSpan,...
    'DecimationFactor', 1) ;% Set to filterUpsample/2 when introducing timing estimation

%% Frequency Distortion and AWGN setup
noiseStrength = 1/2;
freqShift = 0*10^3; 

%% System loop

filteredData = [];%zeros(length(modulatedData)*2,1);
allfilteredTXData = [];
allnoisyData = [];
alloffsetData = [];
allfilteredData = [];
allfilteredDataRef = [];
downsampledRxData = [];
allDownsampledRxData = [];
allDemodulatedData = [];
modtemp = [];

nFD= [];
filteredData = [];
pllData = [];
pll = comm.SymbolSynchronizer(...
    'SamplesPerSymbol', sPS, ...
    'NormalizedLoopBandwidth', 0.01, ...
    'DampingFactor', 1.0, ...
    'DetectorGain', 10, ...
    'TimingErrorDetector', 'Zero-Crossing (decision-directed)');

for k=1:frameSize:(numSamples)
    
    timeIndex = (k:k+frameSize-1).';
    
    % Filter signal
    for fltfix= 1:sPS
    modtemp = [modulatedData(timeIndex);0;0;0;0]; %4
    end
    filteredTXData = step(TxFlt, modtemp);
    for fltfix= 1:(sPS*filterSymbolSpan)
    filteredTXData =[filteredTXData; 0;]; %10
    end
    allfilteredTXData = [allfilteredTXData; filteredTXData];
    
    % Pass through channel
    noiseChan = noiseStrength*randn([1 length(filteredTXData)]); %%AWGN
    noisyData = filteredTXData + noiseChan';
    allnoisyData = [allnoisyData; noisyData];
    
    % Time delay signal
%     offsetData = step(varDelay, noisyData, timingOffset); % Variable
    offset= exp(1i*2*pi*freqShift*[1:length(filteredTXData)]); %rx_mystery_waveform.*exp(1i*2*pi*.0195*[1:length(rx_mystery_waveform)]);
    for i = 1:length(filteredTXData)
    offsetData(i) = noisyData(i) *offset(i);
    end
    if k==1
    offsetData= offsetData';
    end
    alloffsetData = [alloffsetData, offsetData];

    % Filter signal
    %filteredData(timeIndex) = step(RxFlt, offsetData);
    nFD = step(RxFlt, offsetData);
    filteredData = [filteredData; nFD]; %#ok<AGROW>
%     filteredDataRef = step(RxFlt, noisyData);
%     allfilteredDataRef = [allfilteredDataRef; filteredDataRef];
    pllData = pll(step(RxFlt,offsetData));
    pllDataRef = real(pll(step(RxFlt,offsetData)));

    step(cdPLL,pllData);
%% PLL
    
    %downsample filtered data
    %by Sean Brady
    for i =1:frameSize
        downsampledRxData = [downsampledRxData,filteredData((sPS*filterSymbolSpan)+1+4*(i-1))];
    end

    %Demod
    %by Sean Brady
    
    demodulatedData = demod.step(downsampledRxData');
%    allDemodulatedData = [allDemodulatedData,demodulatedData'];
    tdemodulatedData = demod.step(modulatedData);
    % Visualize Error
%     step(cdPre,noisyData);step(cdPost,nFD);pause(0.01); %#ok<*UNRCH>
end
%% Visuals

%testing
t=1:sPS*frameSize;
plot(t,allfilteredTXData(1:sPS*frameSize)); %filterSymbolSpan*sPS =starting point for rx filter output
hold on;
plot(t, filteredData(1:sPS*frameSize),'-o','MarkerIndices',1:sPS*frameSize);
hold on;
plot(t, noisyData(1:sPS*frameSize));
hold off;
legend('TX','RX','NOISE');

cdPre = comm.ConstellationDiagram('ReferenceConstellation', [-1 1],...
    'Name','Baseband');

%% PLL setup

inputDataFull = filteredData;
InputFrameLen = frameSize;

%TED
MaxOutputFrameLen = ceil(InputFrameLen*11/sPS/10);
alpha = 0.5;

b = @(Mu) [alpha*Mu^2 - alpha*Mu;... %Interpolator
    -alpha*Mu^2 - (1-alpha)*Mu + 1;...
    -alpha*Mu^2 + (1+alpha)*Mu;...
    alpha*Mu^2 - alpha*Mu];

LoopFilterState = 0;
LoopPreviousInput = 0;
Trigger = false;
TriggerHistory = false(1, sPS);
Mu = 0;
M1Counter = 0;
InterpFilterState = complex(zeros(3, 1),zeros(3, 1));
TEDBuffer = complex(zeros(1, sPS), zeros(1, sPS));
maxOutputSize = ceil(InputFrameLen*11/double(sPS)/10);
SymbolHolder = complex(zeros(maxOutputSize, 1), zeros(maxOutputSize, 1));

%Loop Filter
DetectorGain = 2.7;
zeta = 1/sqrt(2);
NormalizedLoopBandwidth = 0.1;

% Calculate
Kp = DetectorGain;
K0 = -1;
theta = NormalizedLoopBandwidth/sPS/(zeta + 0.25/zeta);
d  = (1 + 2*zeta*theta + theta^2) * K0 * Kp;
ProportionalGain = (4*zeta*theta) /d;
IntegratorGain   = (4*theta*theta)/d;

%interpolator
filt = dsp.FIRFilter;
filt.NumeratorSource='Input port';
timeCorrect = zeros(size(inputDataFull));
index = 1;
buffer = [0 0 0];

LoopFilter = dsp.IIRFilter( ...
    'Structure', 'Direct form II transposed', ...
    'Numerator', [1 0], 'Denominator', [1 -1]);
NumTrigger  = 0;


