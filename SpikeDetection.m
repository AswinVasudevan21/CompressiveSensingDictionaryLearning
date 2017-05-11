%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SpikeDetection.m
% Given raw data (single channel) and spike detection threshold, computes
% spike snippets

% INPUT VARIABLES: 
% 1. Electrode readings fed from FPGA Hub (d)
% 2. Spike detection threshold (T_detection, output of SpikeDetectionThreshold function)
% INPUT SIZES: 
% 1. d=int32(6,000,000x1) (5 minutes of data at 20kHz);
% 2. T_detection=int32(1x1)


% OUTPUT: Detected spike snippets
% OUTPUT Sizes: spikes=int16(ns x 32), where ns is the number of spikes detected

function spikes = SpikeDetection(d, T_detection)

sf = 20000; % Sampling frequency 20kHz

% Total spike length is 0.0016 seconds - we will take 
before = 0.0003 * sf;    % 0.0003 seconds of data before threshold crossing point and
after = 0.0013 * sf - 1; % 0.0013 seconds of data after threshold crossing point.


% Detection Method : Negative Amplitude Thresholding
% Find where signal crosses threshold below 0
foundInd = find(-1*d > T_detection);
maxCovered = 0; % Largest index covered so far in spikes

spikes = zeros(10000,32);
spikeCount = 1;

% Aligning Spikes
for a = 1:length(foundInd)
    % If the threshold crossing points are already included in the previous spikes,
    % continue to the next index (this resolves overlapping spikes in case threshold crossing points are close together).
    if maxCovered >= foundInd(a)
        continue;
    end
    
    fromInd = foundInd(a)-before;
    toInd = foundInd(a)+after;
    
    if fromInd <= 0 || toInd > length(d) || fromInd < maxCovered
        continue;
    end
    
    
    spikes(spikeCount,:) = d(fromInd:toInd)';
    
    spikeCount = spikeCount+1;
    
    maxCovered = toInd;

end

spikes = spikes(1:spikeCount-1,:);

end


