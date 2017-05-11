%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SpikeDetectionThreshold.m
% Given raw data (single channel), computes spike detection threshold

% INPUT VARIABLES: Electrode readings fed from FPGA Hub (d)
% INPUT SIZES: d=int16(6,000,000x1) (5 minutes of data at 20kHz);

% OUTPUT: Threshold
% Threshold = 5 * median(abs(data))/0.6745

function T_detection = SpikeDetectionThresholdForRecon(d)

T_detection = 5 * median(abs(d))/0.2;

end


