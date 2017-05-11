function y_snle = snle(y, goodWires, varargin )
%
% usage: y_snle = snle( y, varargin )
%
% INPUTS:
%   y - input data, m x n where m is the number of wires and n is the
%       number of points
%   goodWires - which wires have good data for which it is worth
%       calculating the smoothed non-linear energy. This is a vector with a
%       true value for each good wire, and a false/zero value for each bad
%       wire
%
% varargs:
%   windowsize - length of the windowing function, in samples; default = 12
%   windowfunction - type of windowing function to use; default = 'hanning'
%
% OUTPUTS:
%   y_snle - the smoothed nonlinear energy of the input data
%
% From Mukhopadyay and Ray, "A New Interpretation of Nonlinear
%   Energy Operator and Its Efficacy is Spike Detection",
%   IEEE Trans Biomed Eng, 1998

windowSize = 10;
windowFunct = 'hanning';
snlePeriod = 3;

for iarg = 1 : 2 : nargin - 2
    switch lower(varargin{iarg})
        case 'windowsize',
            windowSize = varargin{iarg + 1};
        case 'windowfunction',
            windowFunct = varargin{iarg + 1};
        case 'snlePeriod'
            snlePeriod = varargin{iarg + 1};
    end
    
end

numSamps = size(y, 2);

% make sure the window length is odd so there is no phase shift in the
% filtered data
if round(windowSize / 2) == windowSize / 2
    windowSize = windowSize + 1;
end
phaseShift = ceil(windowSize / 2);

% can add more options for the window later; see option in the signal
% processing toolbox
switch lower(windowFunct)
    case 'hanning',
        w = hann(windowSize);
end

y_snle = y;
L = length(y)-(2*snlePeriod);
% T variable added to extend window, 20150201 Matt Gaidica
y_snle(:,1+snlePeriod:end-snlePeriod) = y(:,1+snlePeriod:end-snlePeriod).^2 - y(:,1:L) .* y(:,(end-L+1):end);
% old code below
%y_snle(:, 2 : end-1) = y(:, 2 : end-1).^2 - y(:, 1 : end-2) .* y(:, 3 : end);

for iRow = 1 : size(y, 1)
    if goodWires(iRow)
        temp = conv(y_snle(iRow, :), w);
        y_snle(iRow, :) = temp(phaseShift : phaseShift + numSamps - 1);
    end
end