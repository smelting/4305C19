clear all

%% Tx system specs
sampleRateHz = 1e6; % Sample rate
sPS = 2; %samples per symbols
frameSize = 73;
numFrames =10;
numSamples = numFrames*frameSize; % Samples to simulate
modulationOrder = 2; %only 4 and 2
filterSymbolSpan = 4;
encoderState=0;
encodedData = [];

%% Packet Structure

%9 bit preamble
preamble = [1;1;1;1;1;1;0;0;0];
%32b ID
ID = randi([0 1], 32, 1);
%8b temp
temp = randi([0 1], 8, 1);
%8b pressure
pressure = randi([0 1], 8, 1);
%8b flags= error message
flag = randi([0 1], 8, 1);
%8b CRC = check for data corruption, need to find in real life
CRC = randi([0 1], 8, 1);
%73b frameSize


%% Data setup
data = [];
frameData = [preamble; ID; temp; pressure; flag; CRC];

for i =1:numFrames
    data = [data; frameData];
end

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
    modl = comm.BPSKModulator();
    modulatedData = modl.step(data);
    demod = comm.BPSKDemodulator();
else
    modl = comm.QPSKModulator('BitInput', true);
    modulatedData = modl.step(data);
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
noiseStrength = 0;
freqShift = 2*10^2; 

%% PLL setup

%TED
MaxOutputFrameLen = ceil(frameSize*11/sPS/10);
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
maxOutputSize = ceil(frameSize*11/double(sPS)/10);
SymbolHolder = complex(zeros(maxOutputSize, 1), zeros(maxOutputSize, 1));

%Loop Filter
DetectorGain = 3.07;
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
%timeCorrect = zeros(size(filteredData));
index = 1;
buffer = [0 0 0];

LoopFilter = dsp.IIRFilter( ...
    'Structure', 'Direct form II transposed', ...
    'Numerator', [1 0], 'Denominator', [1 -1]);



%% System loop

filteredData = [];%zeros(length(modulatedData)*2,1);
allfilteredTXData = [];
allnoisyData = [];
alloffsetData = [];
allfilteredData = [];
downsampledRxData = [];
allDownsampledRxData = [];
allDemodulatedData = [];
modtemp = [];

nFD= [];
allintOut= [];
NumTrigger = 0;
tempbER =0;
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
    offset= exp(1i*2*pi*freqShift*[1:length(filteredTXData)]); %rx_mystery_waveform.*exp(1i*2*pi*.0195*[1:length(rx_mystery_waveform)]);
    for i = 1:length(filteredTXData)
    offsetData(i) = noisyData(i) *offset(i);
    end
    if k==1
    offsetData= offsetData';
    end
    alloffsetData = [alloffsetData, offsetData];

    % Filter signal
    nFD = step(RxFlt, offsetData);
    if k== 1
    nFD = nFD(5:length(nFD));
    end
    filteredData = [filteredData; nFD]; %#ok<AGROW>
    
    %PLL
    io=[];
    ir=[];
    for i =1:frameSize+2
%% Using Trigger count as output input count
        NumTrigger = NumTrigger + Trigger;
        %% Interpolator
        % Piecewise parabolic interpolator

        if i == 1
            for ls= 5:2:length(filteredData)
                if real(filteredData(ls)) > 0
                    io = [io,1];
                else
                    io= [io,-1];
                end
            end
        end
        % Update coefficents based on new Mu
        intOut = step(filt,filteredData(sPS*(i-1)+5),fliplr(b(Mu).'));
        if i >2
        allintOut = [allintOut, intOut];
        end
        if real(intOut) > 0
          ir = [ir,1];
        else
          ir= [ir,-1];
        end

        if Trigger % Interpolation output as symbols
            SymbolHolder(NumTrigger) = intOut;
        end
        %% Timing error detector (TED)
        % For zero-crossing, TED wants to make every other sample time
        % aligned with a zero crossing in the eye diagram.  Therefore, for
        % zero-cross to operate we must have 2 SPS.

        % TED calculation occurs on a Trigger
        if Trigger && all(~TriggerHistory(2:end))

            % The above condition allows TED update after a skip. If we want
            % TED update to happen only at regular strobings, need to check
            % TriggerHistory(1) == true in addition to the condition above.

            % Calculate the midsample from interpolator output
            t1 = TEDBuffer(end/2 + 1 - rem(sPS,2));
            t2 = TEDBuffer(end/2 + 1);
            midSample = (t1+t2)/2;

            % Rice Notation
            % x -> real(in)
            % y -> imag(in)
            % a0 -> sign(x)
            % a1 -> sign(y)

            % Data aided method
            e = real(midSample) * (sign(real(TEDBuffer(1))) - sign(real(intOut))) + ...
                imag(midSample) * (sign(imag(TEDBuffer(1))) - sign(imag(intOut)));
        else
            e = 0;
        end

        % Handle bit stuffing/skipping
        switch sum([TriggerHistory(2:end), Trigger])
            case 0
                % No Trigger in awhile
            case 1
                % Update buffer (Normal)
                TEDBuffer = [TEDBuffer(2:end), intOut];
            otherwise % > 1
                % Stuff condition
                TEDBuffer = [TEDBuffer(3:end), 0, intOut];
        end
        %% Loop filter
        % The output v essentially adjust the period of the Trigger

        %loopFiltOut = LoopPreviousInput + LoopFilterState;
        %v = e*ProportionalGain + loopFiltOut;
        %LoopFilterState = loopFiltOut;
        %LoopPreviousInput = e*IntegratorGain;

        v = LoopFilter(e*IntegratorGain);
        v = e*ProportionalGain + v;

        %% Interpolation controller
        % We want to get a Trigger signal every SPS samples, which should
        % happen when we acquire lock
        % Essentially a Trigger should occur at the start of a symbol,
        % therefore based on the start of the symbol we will then utilize
        % the next sample (which should also be within that symbol),
        % interpolate across both with the appropriate delay and produce an
        % output
        %
        % Modulo-1 counter interpolation controller which generates/updates
        % Trigger signal and fractional interpolation interval.
        %
        % Trigger - Trigger edge found and to output data
        % Mu - fractional interpolation interval

        W = v + 1/sPS; % W should be small or == SPS when locked

        TriggerHistory = [TriggerHistory(2:end), Trigger];

        % Determine if we have an underflow, aka we are at the start edge
        % of a sample
        Trigger = (M1Counter < W);

        % With a new underflow we must update interpolator coefficients
        if Trigger
            Mu = M1Counter / W;
        end

        % Update counter
        M1Counter = mod((M1Counter - W), 1);
    end
    ir = ir(3:75);
    kw=0;
    for i = 1:73
        if ir(i)==io(i)
            kw=kw+1;
        end
    end
    tempbER= tempbER + kw/73;
    % Output 1 SPS
    y = SymbolHolder(1:NumTrigger, 1);
    timeCorrect(index:index+length(y)-1) = y;
    index = index + length(y);
    
end

bER= tempbER/numFrames;

demodulatedData = demod.step(allintOut');







