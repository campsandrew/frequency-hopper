%% General simulation parameters
Rsym = 0.2e6;              % Symbol rate in Hertz
ModulationOrder = 4;       % QPSK alphabet size
Interpolation = 2;         % Interpolation factor
Decimation = 1;            % Decimation factor
Tsym = 1 / Rsym;           % Symbol time in sec
Fs = Rsym * Interpolation; % Sample rate

%% Generate Hopping Channels
numChannels = 10;
channelWidth = 1e6;
initFrequency = 2.4e9;
channels = zeros(numChannels, 1);
for i = 1:numChannels
    channels(i) = initFrequency + (i-1) * channelWidth; 
end

%% Tx parameters
RolloffFactor = 0.5;
ScramblerBase = 2;
ScramblerPolynomial = [1 1 1 0 1];
ScramblerInitialConditions = [0 0 0 0];
RaisedCosineFilterSpan = 10;

%% Rx parameters
DesiredPower                  = 2;            % AGC desired output power (in watts)
AveragingLength               = 50;           % AGC averaging length
MaxPowerGain                  = 60;           % AGC maximum output power gain
MaximumFrequencyOffset        = 6e3;
K = 1;
A = 1/sqrt(2);
PhaseRecoveryLoopBandwidth    = 0.01;         % Normalized loop bandwidth for fine frequency compensation
PhaseRecoveryDampingFactor    = 1;            % Damping Factor for fine frequency compensation
TimingRecoveryLoopBandwidth   = 0.01;         % Normalized loop bandwidth for timing recovery
TimingRecoveryDampingFactor   = 1;            % Damping Factor for timing recovery
TimingErrorDetectorGain       = 2.7*2*K*A^2+2.7*2*K*A^2;
PreambleDetectorThreshold     = 0.8;

%% Frame Specifications
BarkerCode = [+1 +1 +1 +1 +1 -1 -1 +1 +1 -1 +1 -1 +1];
ModulatedHeader = sqrt(2)/2 * (-1-1i) * BarkerCode.';
BarkerLength = length(BarkerCode);
HeaderLength = BarkerLength * 2;

%% Message Info
Message = 'ECE 531 Project Message';
MessageLength = length(Message) + 1;
NumberOfMessage = 100;
PayloadLength = NumberOfMessage * MessageLength * 7;
FrameSize = (HeaderLength + PayloadLength) / log2(ModulationOrder);

%% Message generation
msgSet = zeros(NumberOfMessage * MessageLength, 1); 
for msgCnt = 0 : NumberOfMessage - 1
    msgSet(msgCnt * MessageLength + (1 : MessageLength)) = sprintf('%s\n', Message);
end
integerToBit = comm.IntegerToBit(7, 'OutputDataType', 'double');
MessageBits = integerToBit(msgSet);

% For BER calculation masks
BerMask = zeros(NumberOfMessage * length(Message) * 7, 1);
for i = 1 : NumberOfMessage
    BerMask( (i-1) * length(Message) * 7 + ( 1: length(Message) * 7) ) = ...
        (i-1) * MessageLength * 7 + (1: length(Message) * 7);
end

%% Pluto TX
tx = sdrtx(..., 
    'Pluto', ...
    'RadioID',                      'usb:0', ...
    'CenterFrequency',              channels(1), ...
    'BasebandSampleRate',           Fs, ...
    'SamplesPerFrame',              Interpolation * FrameSize, ...
    'Gain',                         0);
BitGenerator = QPSKBitsGenerator( ...
    'NumberOfMessage',              NumberOfMessage, ...
    'MessageLength',                MessageLength, ...
    'MessageBits',                  MessageBits, ...
    'ScramblerBase',                ScramblerBase, ...
    'ScramblerPolynomial',          ScramblerPolynomial, ...
    'ScramblerInitialConditions',   ScramblerInitialConditions);
QPSKModulator = comm.QPSKModulator( ...
    'BitInput',                     true, ...
    'PhaseOffset',                  pi/4, ...
    'OutputDataType',               'double');
TxFilter = comm.RaisedCosineTransmitFilter( ...
    'RolloffFactor',                RolloffFactor, ...
    'FilterSpanInSymbols',          RaisedCosineFilterSpan, ...
    'OutputSamplesPerSymbol',       Interpolation);

%% Pluto RX
rx = sdrrx(..., 
    'Pluto', ...
    'RadioID',                      'usb:0', ...%'sn:104473dc5993001904000f0002c42965db', ...
    'CenterFrequency',              channels(1), ...
    'BasebandSampleRate',           Fs, ...
    'SamplesPerFrame',              Interpolation * FrameSize, ...
    'GainSource',                   'Manual', ...
    'Gain',                         30, ...
    'OutputDataType',               'double');
AGC = comm.AGC(...
    'DesiredOutputPower',           DesiredPower, ...
    'AveragingLength',              AveragingLength, ...
    'MaxPowerGain',                 MaxPowerGain);
RxFilter = comm.RaisedCosineReceiveFilter(...
    'RolloffFactor',                RolloffFactor, ...
    'FilterSpanInSymbols',          RaisedCosineFilterSpan, ...
    'InputSamplesPerSymbol',        Interpolation, ...
    'DecimationFactor',             Decimation);
CoarseFreqEstimator = comm.CoarseFrequencyCompensator(...
    'Modulation',                   'QPSK', ...
    'Algorithm',                    'Correlation-based', ...
    'MaximumFrequencyOffset',       MaximumFrequencyOffset, ....
    'SampleRate',                   Fs / Decimation);
CoarseFreqCompensator = comm.PhaseFrequencyOffset(...
    'PhaseOffset',                  0, ...
    'FrequencyOffsetSource',        'Input port', ...
    'SampleRate',                   Fs / Decimation);
FineFreqCompensator = comm.CarrierSynchronizer(...
    'Modulation',                   'QPSK', ...
    'ModulationPhaseOffset',        'Auto', ...
    'SamplesPerSymbol',             Interpolation / Decimation, ...
    'DampingFactor',                PhaseRecoveryDampingFactor, ...
    'NormalizedLoopBandwidth',      PhaseRecoveryLoopBandwidth);
TimingRec = comm.SymbolSynchronizer(...
    'TimingErrorDetector',          'Gardner (non-data-aided)', ...
    'SamplesPerSymbol',             Interpolation / Decimation, ...
    'DampingFactor',                TimingRecoveryDampingFactor, ...
    'NormalizedLoopBandwidth',      TimingRecoveryLoopBandwidth, ...
    'DetectorGain',                 TimingErrorDetectorGain);
PrbDet = comm.PreambleDetector(ModulatedHeader, ...
    'Input',                        'Symbol', ...
    'Threshold',                    PreambleDetectorThreshold);
FrameSync = FrameSynchronizer(...
    'OutputFrameLength',            FrameSize, ...
    'PreambleLength',               HeaderLength / 2);
QPSKDemodulator = comm.QPSKDemodulator(...
    'PhaseOffset', pi/4, ...
    'BitOutput', true);
Descrambler = comm.Descrambler(...
    ScramblerBase, ...
    ScramblerPolynomial, ...
    ScramblerInitialConditions);
ErrorRateCalc = comm.ErrorRate( ...
    'Samples', 'Custom', ...
    'CustomSamples', BerMask);
BitToInteger = comm.BitToInteger(7, 'OutputDataType', 'int8');
IntegerToBit = comm.IntegerToBit(7, 'OutputDataType', 'double');

%% Simulation Variables
state = 0;
StopTime = 100;
StopSettleTime = 5;
currentTime = 0;
currentSettleTime = StopSettleTime;
FrameTime = Interpolation * FrameSize / Fs;
MeanFreqOff = 0;
Cnt = 0;
BERThreshold = 0.01;

%% Main Simulation Loop
while 1
    
    % State Machine
    if state == 0
        % Generate Data
        [transmittedBin, ~] = BitGenerator();
        modulatedData = QPSKModulator(transmittedBin);           
        transmittedSignal = TxFilter(modulatedData);
        
        % Send Transmission
        tx.transmitRepeat(transmittedSignal);
        state = 1;
    elseif state == 1
        % Get Received Signal through Matched Filter
        RCRxSignal = RxFilter(AGC(rx()));
        
        % Coarse Frequency Estimation and Compensation
        [~, freqOffsetEst] = CoarseFreqEstimator(RCRxSignal);
        freqOffsetEst = (freqOffsetEst + Cnt*MeanFreqOff) / (Cnt + 1);
        Cnt = Cnt + 1;
        MeanFreqOff = freqOffsetEst;
        coarseCompSignal = CoarseFreqCompensator(RCRxSignal, -freqOffsetEst);
        
        % Timing Corrected Signal
        timingRecSignal = TimingRec(coarseCompSignal);

        % Fine Frequency Compensation
        fineCompSignal = FineFreqCompensator(timingRecSignal);

        % Frame Syncronization
        [prbIdx, dtMt] = PrbDet(fineCompSignal);
        [symFrame, isFrameValid] = FrameSync(fineCompSignal, prbIdx, dtMt);
        
        % Get Message and Bit Error Rate
        if isFrameValid
            % Phase ambiguity estimation and demodulation
            phaseEst = round(angle(mean(conj(ModulatedHeader) .* symFrame(1:HeaderLength/2)))*2/pi)/2*pi;
            phShiftedData = symFrame .* exp(-1i*phaseEst);
            demodOut = QPSKDemodulator(phShiftedData);

            % Performs descrambling on only payload part
            deScrData = Descrambler(demodOut(HeaderLength + (1:PayloadLength)));

            % Get bit error rate and message
            charSet = BitToInteger(deScrData);
            message = sprintf('%s', char(charSet));
            BER = ErrorRateCalc(MessageBits, deScrData);            
            disp("BER: " + BER(1));
%             disp(message);

            if BER(1) > BERThreshold && currentSettleTime > StopSettleTime
                state = 2;
            end
        end
    elseif state == 2
        % Method One Random Selection
        %lo = RandomChannel(channels);
        
        % Method Two Analyze All Channels and Pick Best Channel
        lo = BestChannel(channels, tx, rx, Fs);

        tx.CenterFrequency = lo;
        rx.CenterFrequency = lo;
        reset(ErrorRateCalc);
        release(tx);
        release(rx);
        state = 0;
        currentSettleTime = 0;
        disp("Channel Switched: " + lo);
        pause(2000);
    else
        release(tx);
        release(rx);
        break;
    end
    
    % Check if simulation is done
    currentSettleTime = currentSettleTime + FrameTime;
    currentTime = currentTime + FrameTime;
    if currentTime > StopTime
        state = -1;
    end
end

%% Random Channel Selection Method
function lo = RandomChannel(channels)
    lo = channels(randi(length(channels)));
end

%% Best Channel Selection Method
function lo = BestChannel(channels, tx, rx, Fs)
    % Do spectral analysis
    scope = dsp.SpectrumAnalyzer('SampleRate', Fs, 'SpectralAverages', 100);
    values = zeros(1,length(channels));
    frames = 50;
    for i = 1:length(channels)
        channel = channels(i);
        %Make tx channel different?
        tx.CenterFrequency = channel;
        rx.CenterFrequency = channel;
        %Make sure to transmit nothing to avoid noise
        tx.transmitRepeat(zeros(tx.SamplesPerFrame,1));
        %Grab some samples for analysis
        for j = 1:frames
            scope(rx());
        end
        data = getSpectrumData(scope);
        %Center is at 768
        spectrum = data{:,{'Spectrum'}}{1};
        %frequency = data{:,{'FrequencyVector'}}{1};
        %Store spectral value at channel
        values(i) = spectrum(768);
    end
    [val, idx] = min(values);
    lo = channels(idx);
end