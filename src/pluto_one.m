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
BarkerLength = length(BarkerCode);
HeaderLength = BarkerLength * 2;

%% Message Info
Message = 'Hello World1';
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

%     'RadioID',                      'sn:104473dc599300131e00150082672a34e9', ...
%% Pluto TX
tx = sdrtx(..., 
    'Pluto', ...
    'RadioID',                      'usb:0', ...%'sn:104473dc5993001904000f0002c42965db', ...
    'CenterFrequency',              915e6, ...
    'BasebandSampleRate',           Fs, ...
    'SamplesPerFrame',              Interpolation * FrameSize, ...
    'Gain',                         0);
hTx = QPSKTransmitter(...
    'UpsamplingFactor',             Interpolation, ...
    'RolloffFactor',                RolloffFactor, ...
    'RaisedCosineFilterSpan',       RaisedCosineFilterSpan, ...
    'MessageBits',                  MessageBits, ...
    'MessageLength',                MessageLength, ...
    'NumberOfMessage',              NumberOfMessage, ...
    'ScramblerBase',                ScramblerBase, ...
    'ScramblerPolynomial',          ScramblerPolynomial, ...
    'ScramblerInitialConditions',   ScramblerInitialConditions);

%% Pluto RX
rx = sdrrx(..., 
    'Pluto', ...
    'RadioID',                      'usb:0', ...%'sn:104473dc5993001904000f0002c42965db', ...
    'CenterFrequency',              915e6, ...
    'BasebandSampleRate',           Fs, ...
    'SamplesPerFrame',              Interpolation * FrameSize, ...
    'GainSource',                   'Manual', ...
    'Gain',                         30, ...
    'OutputDataType',               'double');
hRx  = QPSKReceiver(...
    'ModulationOrder',               ModulationOrder, ...
    'SampleRate',                    Fs, ...
    'DecimationFactor',              Decimation, ...
    'FrameSize',                     FrameSize, ...
    'HeaderLength',                  HeaderLength, ...
    'NumberOfMessage',               NumberOfMessage, ...
    'PayloadLength',                 PayloadLength, ...
    'DesiredPower',                  DesiredPower, ...
    'AveragingLength',               AveragingLength, ...
    'MaxPowerGain',                  MaxPowerGain, ...
    'RolloffFactor',                 RolloffFactor, ...
    'RaisedCosineFilterSpan',        RaisedCosineFilterSpan, ...
    'InputSamplesPerSymbol',         Interpolation, ...
    'MaximumFrequencyOffset',        MaximumFrequencyOffset, ...
    'PostFilterOversampling',        Interpolation/Decimation, ...
    'PhaseRecoveryLoopBandwidth',    PhaseRecoveryLoopBandwidth, ...
    'PhaseRecoveryDampingFactor',    PhaseRecoveryDampingFactor, ...
    'TimingRecoveryDampingFactor',   TimingRecoveryDampingFactor, ...
    'TimingRecoveryLoopBandwidth',   TimingRecoveryLoopBandwidth, ...
    'TimingErrorDetectorGain',       TimingErrorDetectorGain, ...
    'PreambleDetectorThreshold',     PreambleDetectorThreshold, ...
    'DescramblerBase',               ScramblerBase, ...
    'DescramblerPolynomial',         ScramblerPolynomial, ...
    'DescramblerInitialConditions',  ScramblerInitialConditions,...
    'BerMask',                       BerMask, ...
    'PrintOption',                   true);

%% State Variables
state = 0;
TxStopTime = 15;
RxStopTime = 10;
rounds = 2;
currentTime = 0;
currentRound = 0;
FrameTime = Interpolation * FrameSize / Fs;

%% TX/RX State Machine
while 1
    
    if state == 0
        data = step(hTx);
        step(tx, data);
        
        if currentTime >= TxStopTime
            disp("SWITCH TO RX");
            currentTime = 0;
            state = 1;
        end
    elseif state == 1
        rcvdSignal = rx();
        [~, ~, ~, ~, message] = hRx(rcvdSignal);
        disp(message);
        
        if currentTime >= RxStopTime
            rounds = rounds + 1;
            currentTime = 0;
            state = 0;
        end
    else
        disp("CLEANUP");
        release(tx);
        release(rx);
        release(hRx);
        release(hTx);
        break
    end
    
    if currentRound >= rounds
        state = 2;
    end
    
    currentTime = currentTime + FrameTime;
    
%     if currentTime >= TxStopTime
%         tx.Gain = -80;
% %         step(tx, complex(zeros(8426, 1), zeros(8426, 1)));
%         rcvdSignal = rx();
%         [~, ~, ~, ~, message] = hRx(rcvdSignal);
%         disp(message);
%     else
%         tx.Gain = 0;
%         
%     end
%     
%     
%     message = receive(rx, hRx);
%     disp(message);
    
%     if state == 0
%         transmit(tx, hTx);
%         
%         if currentTime >= TxStopTime
%             currentTime = 0;
%             state = 1;
%         end
%     elseif state == 1
%         message = receive();
%         
%         if contains(message, 'Freq')
%             currentTime = 0;
%             currentRound = currentRound + 1;
%             state = 0;
%         end
%         
%         if currentTime >= RxStopTime || contains(message, 'None')
%             currentTime = 0;
%             currentRound = currentRound + 1;
%             state = 0;
%         end
%     elseif state == 2
%         state = 3;
%     else
%         release(hRx);
%         release(hTx);
%         release(rx);
%         release(tx);
%         break
%     end
%     
%     if currentRound >= rounds
%         state = 2;
%     end
    
%     currentTime = currentTime + FrameTime;
end