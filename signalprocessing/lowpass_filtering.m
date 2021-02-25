function output_signal = lowpass_filtering(input_signal, varargin)
% -------------------------------------------------------
%
%    LOWPASS_FILTERING - Lowpass filtering of a signal of any kind
%
%    Ver. 1.0.0
%
%    Created:       Alexander Kramlich        (17.07.2015)
%    Last modified: Alexander Kramlich        (19.07.2015)
%
%    Institute of Biomedical Engineering
%    Universitaet Karlsruhe (TH)
%
%    http://www.ibt.kit.edu
%
%    Copyright 2000-2015 - All rights reserved.
%
% ------------------------------------------------------
%
%   LOWPASS_FILTERING(input_signal)
%   LOWPASS_FILTERING(input_signal, 'parameterName', 'parameterValue')
%
%                           --- Input Arguments ---
%   input_signal:
%       Single signal (e.g. NLEO of an ECG) or several signals in a matrix
%       of dimension MxN, where M is the total number of signals and N the
%       length of one single signal. The given function operates along the
%       second dimension of the matrix.
%
%   Name-Value Pair Arguments:
%       'FilterType'
%           --- 'iir' / 'fir' (default) ---
%           Type of the filter.
%       'HalfPowerFrequency'
%           Cutofffrequency of the low-pass postfilter.
%       	Default value: 24 Hz (Nguyen et al. 2008)
%           Note: HalfPowerFrequency (IIR) -- 3dB,
%                 CutoffFrequency (FIR) -- 6dB.
%       'FilterOrder'
%           Filter order of the filter.
%           Defalt value: 2 (Nguyen et al. 2008)
%
%                           --- Output Arguments ---
%   output_signal:
%       Filtered signal of the dimension of the input signal.
%

%% Validation of input
p = inputParser;
defaultFilterType = 'fir';
defaultHalfPowerFrequency = 24;
defaultFilterOrder = 2;
defaultSamplingRate = 2035; 
expectedModes = {'iir', 'fir'};

addRequired(p, 'input_signal',...
            @(x) validateattributes(x, {'numeric'}, {'nonempty', '2d'}));
  
addParameter(p, 'FilterType',...
                defaultFilterType,...
                @(x) any(validatestring(x, expectedModes)));
addParameter(p, 'HalfPowerFrequency',...
                defaultHalfPowerFrequency,...
                @(x) validateattributes(x, {'numeric'}, {'nonempty', '>', 0}));
addParameter(p, 'FilterOrder',...
                defaultFilterOrder,...
                @(x) validateattributes(x, {'numeric'}, {'nonempty', '>', 0}));
addParameter(p, 'SamplingRate',...
                defaultSamplingRate,...
                @(x) validateattributes(x, {'numeric'}, {'nonempty', '>', 0}));
            
parse(p, input_signal, varargin{:})

%% Function's local variables
% Number of samples in one single signal
N = size(p.Results.input_signal, 2);

% Cutofffrequency of the low-pass filter
f_g = p.Results.HalfPowerFrequency;

% Filter Order
n = p.Results.FilterOrder;

% Sampling rate
f_s = p.Results.SamplingRate;

%% Low-pass filtering
switch p.Results.FilterType
    case 'iir'
        % Infinite impulse response filter
        LP2 = designfilt('lowpassiir',...
                         'FilterOrder', n,...
                         'HalfPowerFrequency', f_g,...
                         'SampleRate', f_s,...
                         'DesignMethod', 'butter');
        
        output_signal = padarray(input_signal, [0 2*N]);        % Add padding
        output_signal = filtfilt(LP2, double(output_signal)')';
        output_signal = output_signal(:, 2*N+1:3*N);            % Remove padding
    case 'fir'
        % Finite impulse response filter
        LP2 = designfilt('lowpassfir',...
                         'FilterOrder', n,...
                         'CutoffFrequency', f_g,...
                         'SampleRate', f_s);
                    
        output_signal = filtfilt(LP2, input_signal')';
    otherwise
        % No filtering
end
