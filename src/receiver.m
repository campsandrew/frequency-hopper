%% General simulation parameters
Rsym = 0.2e6;             % Symbol rate in Hertz
ModulationOrder = 4;      % QPSK alphabet size
Interpolation = 2;        % Interpolation factor
Decimation = 1;           % Decimation factor
Tsym = 1/Rsym;  % Symbol time in sec
Fs   = Rsym * Interpolation; % Sample rate

%% Frame Specifications
% [BarkerCode*2 | 'Hello world 000\n' | 'Hello world 001\n' ... | 'Hello world 099\n'];
BarkerCode      = [+1 +1 +1 +1 +1 -1 -1 +1 +1 -1 +1 -1 +1];     % Bipolar Barker Code
BarkerLength    = length(BarkerCode);
HeaderLength    = BarkerLength * 2;                   % Duplicate 2 Barker codes to be as a header
Message         = 'Hello world';
MessageLength   = length(Message) + 5;                % 'Hello world 000\n'...
NumberOfMessage = 1;                                  % Number of messages in a frame
PayloadLength   = NumberOfMessage * MessageLength * 7; % 7 bits per characters
FrameSize       = (HeaderLength + PayloadLength)/ log2(ModulationOrder);                                    % Frame size in symbols
FrameTime       = Tsym*FrameSize;

%% Rx parameters
RolloffFactor     = 0.5;                      % Rolloff Factor of Raised Cosine Filter
ScramblerBase     = 2;
ScramblerPolynomial           = [1 1 1 0 1];
ScramblerInitialConditions    = [0 0 0 0];
RaisedCosineFilterSpan = 10;                  % Filter span of Raised Cosine Tx Rx filters (in symbols)
DesiredPower                  = 2;            % AGC desired output power (in watts)
AveragingLength               = 50;           % AGC averaging length
MaxPowerGain                  = 60;           % AGC maximum output power gain
MaximumFrequencyOffset        = 6e3;
% Look into model for details for details of PLL parameter choice. 
% Refer equation 7.30 of "Digital Communications - A Discrete-Time Approach" by Michael Rice.
K = 1;
A = 1/sqrt(2);
PhaseRecoveryLoopBandwidth    = 0.01;         % Normalized loop bandwidth for fine frequency compensation
PhaseRecoveryDampingFactor    = 1;            % Damping Factor for fine frequency compensation
TimingRecoveryLoopBandwidth   = 0.01;         % Normalized loop bandwidth for timing recovery
TimingRecoveryDampingFactor   = 1;            % Damping Factor for timing recovery
% K_p for Timing Recovery PLL, determined by 2KA^2*2.7 (for binary PAM),
% QPSK could be treated as two individual binary PAM,
% 2.7 is for raised cosine filter with roll-off factor 0.5
TimingErrorDetectorGain       = 2.7*2*K*A^2+2.7*2*K*A^2;
PreambleDetectorThreshold     = 0.8;

%% Message generation and BER calculation parameters
msgSet = zeros(100 * MessageLength, 1);
for msgCnt = 0 : 99
    msgSet(msgCnt * MessageLength + (1 : MessageLength)) = ...
        sprintf('%s %03d\n', Message, msgCnt);
end
integerToBit = comm.IntegerToBit(7, 'OutputDataType', 'double');
MessageBits = integerToBit(msgSet);

% For BER calculation masks
BerMask = zeros(NumberOfMessage * length(Message) * 7, 1);
for i = 1 : NumberOfMessage
    BerMask( (i-1) * length(Message) * 7 + ( 1: length(Message) * 7) ) = ...
        (i-1) * MessageLength * 7 + (1: length(Message) * 7);
end

% Pluto receiver parameters
PlutoCenterFrequency      = 915e6;
PlutoGain                 = 30;
PlutoFrontEndSampleRate   = Fs;
PlutoFrameLength          = Interpolation * FrameSize;

% Experiment parameters
PlutoFrameTime = PlutoFrameLength / PlutoFrontEndSampleRate;
StopTime = 10;

rx  = QPSKReceiver(...
    'ModulationOrder', ModulationOrder, ...
    'SampleRate', Fs, ...
    'DecimationFactor', Decimation, ...
    'FrameSize', FrameSize, ...
    'HeaderLength', HeaderLength, ...
    'NumberOfMessage', NumberOfMessage, ...
    'PayloadLength', PayloadLength, ...
    'DesiredPower', DesiredPower, ...
    'AveragingLength', AveragingLength, ...
    'MaxPowerGain', MaxPowerGain, ...
    'RolloffFactor', RolloffFactor, ...
    'RaisedCosineFilterSpan', RaisedCosineFilterSpan, ...
    'InputSamplesPerSymbol', Interpolation, ...
    'MaximumFrequencyOffset', MaximumFrequencyOffset, ...
    'PostFilterOversampling', Interpolation/Decimation, ...
    'PhaseRecoveryLoopBandwidth', PhaseRecoveryLoopBandwidth, ...
    'PhaseRecoveryDampingFactor', PhaseRecoveryDampingFactor, ...
    'TimingRecoveryDampingFactor', TimingRecoveryDampingFactor, ...
    'TimingRecoveryLoopBandwidth', TimingRecoveryLoopBandwidth, ...
    'TimingErrorDetectorGain', TimingErrorDetectorGain, ...
    'PreambleDetectorThreshold', PreambleDetectorThreshold, ...
    'DescramblerBase', ScramblerBase, ...
    'DescramblerPolynomial', ScramblerPolynomial, ...
    'DescramblerInitialConditions', ScramblerInitialConditions,...
    'BerMask', BerMask,...
    'PrintOption',true);
    
% Create and configure the Pluto System object.
radio = sdrrx('Pluto');
radio.CenterFrequency       = PlutoCenterFrequency;
radio.BasebandSampleRate    = PlutoFrontEndSampleRate;
radio.SamplesPerFrame       = PlutoFrameLength;
radio.GainSource            = 'Manual';
radio.Gain                  = PlutoGain;
radio.OutputDataType        = 'double';
radio.RadioID = 'usb:1';


% Initialize variables
currentTime = 0;
BER = [];
rcvdSignal = complex(zeros(PlutoFrameLength,1));
while currentTime <  StopTime
    rcvdSignal = radio();
    [~, ~, ~, BER] = rx(rcvdSignal);
    currentTime= currentTime + (radio.SamplesPerFrame / radio.BasebandSampleRate);
end

release(rx);
release(radio);

