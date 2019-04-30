%% General simulation parameters
Rsym = 0.2e6;              % Symbol rate in Hertz
ModulationOrder = 4;       % QPSK alphabet size
Interpolation = 2;         % Interpolation factor
Decimation = 1;            % Decimation factor
Tsym = 1 / Rsym;           % Symbol time in sec
Fs = Rsym * Interpolation; % Sample rate
FrameSize = 1024;
FrameTime = Interpolation * FrameSize / Fs;

%% Generate Hopping Channels
numChannels = 5;
channelWidth = 1e6;
initFrequency = 2.4e9;
channels = zeros(numChannels, 1);
for i = 1:numChannels
    channels(i) = initFrequency + (i-1) * channelWidth; 
end

%% Pluto TX
tx = sdrtx(..., 
    'Pluto', ...
    'RadioID',                      'usb:0', ...
    'CenterFrequency',              channels(4), ...
    'BasebandSampleRate',           Fs, ...
    'SamplesPerFrame',              Interpolation * FrameSize, ...
    'Gain',                         0);

%% Simulation Variables
state = 0;
timePerChannel = 10;
maxChanges = 100;
currentTime = 0;
numChanges = 0;
currentChannel = 1;
%% Main Simulation Loop
while 1
    
    % State Machine
    if state == 0
        data = complex(rand(Interpolation * FrameSize, 1), rand(Interpolation * FrameSize, 1));
        tx(data);
        
        if currentTime > timePerChannel
            state = 1;
        end
    elseif state == 1
        currentTime = 0;
        currentChannel = currentChannel + 1;
        lo = channels(mod(currentChannel,length(channels))+1);
        %tx.CenterFrequency = lo;
        numChanges = numChanges + 1;
        state = 0;
        disp("Channel Switched: " + lo);
    else
        release(tx);
        break;
    end
    
    % Check if simulation is done
    currentTime = currentTime + FrameTime;
    if numChanges > maxChanges
        state = -1;
    end
end

%% Random Channel Selection
function lo= RandomChannel(channels)
    lo = channels(randi(length(channels)));
end