function output_signal = method_wavelet(input_signal, varargin)
% -------------------------------------------------------
%
%    METHOD_WAVELET - Activity detection algorithm based on the wavelet
%                     transform of the signal
%
%    Ver. 1.3.1
%
%    Created:       Alexander Kramlich        (13.08.2015)
%    Last modified: Alexander Kramlich        (21.01.2016)
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
%   Steps:
%   1. Wavelet transform
%   2. Low-pass filter
%   3. Thresholding
%   4. Postprocessing
%
% ------------------------------------------------------
%
%   WAVELET_TRANSFORM(signal)
%   WAVELET_TRANSFORM(signal, 'parameterName', 'parameterValue')
%
%                           --- Input Arguments ---
%   signal:
%       Single signal (e.g. EGM, ECG, etc.) or several signals in a matrix
%       of dimension MxN, where M is the total number of signals and N the
%       length of one single signal.
%   
%   Name-Value Pair Arguments:
%       'WaveletShape'
%          ---
%           'db1' / 'db2' / 'db3' / 'db4' / 'db5' /
%           'coif1' / 'coif2' / 'coif3' / 'coif4' / 'coif5' /
%           'sym2' / 'sym3' / 'sym4' / 'sym5' /
%           'sym6' / 'sym7' / 'sym8' / 'sym9' / 'sym10' /
%           'bior1.1' / 'bior1.3' / 'bior1.5' (default) /
%           'bior2.2' / 'bior2.4' / 'bior2.6' / 'bior2.8' /
%           'bior3.1' / 'bior3.3' / 'bior3.5' / 'bior3.7' / 'bior3.9' /
%           'bior4.4' / 'bior5.5' / 'bior6.8' /
%           ---
%           Shape of the wavelet used for the wavelet transform.
%       'HalfPowerFrequency'
%           Cutofffrequency of the second low-pass filter.
%           Default value: 40 Hz (optimum of the performance evaluation
%           with bior1.5)
%       'ThresholdFactor'
%           Threshold factor k which is multiplied with the standard
%           deviation of the signal to calculate the threshold.
%           T = k * std(signal)
%           Threshold is then used to obtain the stepfunction by comparison
%           of the signal with the threshold.
%           Default value: 0.9 (optimum of the performance evaluation
%           with bior1.5)
%       'Thresholding'
%           --- 'on' (default) / 'off' ---
%           Option to switch off thresholding.
%           Note: Switching off thresholding also switches off
%           postprocessing.
%       'Postprocessing'
%           --- 'on' (default) / 'off' ---
%           Option to switch off postprocessing.
%
%                           --- Output Arguments ---
%   output_signal:
%       Stepfunction, which is a logic array with a length of the original
%       signal. 1 (true) indicates the presence of activity, whereas
%       0 (false) stands for inactivity.
%       The value of the stepfunction is set to 1 (true), if the input
%       signal (e.g. NLEO of an ECG) is above threshold. Otherwise the
%       value of the stepfunction is 0 (false).
%

%% Validation of input
p = inputParser;
defaultWaveletShape = 'bior1.5';
defaultHalfPowerFrequency = 35;
defaultThresholdFactor = 0.4;
defaultThresholding = 'on';
defaultPostprocessing = 'on';
defaultMinimumInaktivityLength = 42;            % Schilling: 42
defaultMinimumAktivityLength = 10;              % Schilling: 10
defaultSamplerate = 2035;

expectedWaveletShapes = {'db1', 'db2', 'db3', 'db4', 'db5',...             % Daubechies
                         'coif1', 'coif2', 'coif3', 'coif4', 'coif5',...   % Coiflets
                         'sym2', 'sym3', 'sym4', 'sym5',...                % Symlets
                         'sym6', 'sym7', 'sym8', 'sym9', 'sym10',...
                         'bior1.1', 'bior1.3', 'bior1.5',...               % Biorthogonals
                         'bior2.2', 'bior2.4', 'bior2.6', 'bior2.8',...
                         'bior3.1', 'bior3.3', 'bior3.5', 'bior3.7', 'bior3.9',...
                         'bior4.4',...
                         'bior5.5',...
                         'bior6.8'};
expectedModes = {'on', 'off'};

addRequired(p, 'input_signal', @(x) validateattributes(x, {'numeric'}, {'nonempty', '2d'}));

addParameter(p, 'WaveletShape',...
                defaultWaveletShape,...
                @(x) any(validatestring(x, expectedWaveletShapes)));
addParameter(p, 'HalfPowerFrequency',...
                defaultHalfPowerFrequency,...
                @(x) validateattributes(x, {'numeric'}, {'nonempty', '>', 0}));
addParameter(p, 'ThresholdFactor',...
                defaultThresholdFactor,...
                @(x) validateattributes(x, {'numeric'}, {'nonempty', '>', 0}));
addParameter(p, 'Thresholding',...
                defaultThresholding,...
                @(x) any(validatestring(x, expectedModes)));
addParameter(p, 'Postprocessing',...
                defaultPostprocessing,...
                @(x) any(validatestring(x, expectedModes)));
addParameter(p, 'samplerate',...
                defaultSamplerate,...
                @(x) validateattributes(x, {'numeric'}, {'nonempty', '>', 0}));
addParameter(p, 'MinimumInaktivityLength',...
                defaultMinimumInaktivityLength,...
                @(x) validateattributes(x, {'numeric'}, {'nonempty', '>', 0}));
addParameter(p, 'MinimumAktivityLength',...
                defaultMinimumAktivityLength,...
                @(x) validateattributes(x, {'numeric'}, {'nonempty', '>', 0}));
            
parse(p, input_signal, varargin{:})

%% Function's local variables
% Wavelet Shape
wavelet_shape = p.Results.WaveletShape;

% Halfpower frequency
f_g = p.Results.HalfPowerFrequency;

% Threshold Factor
k = p.Results.ThresholdFactor;

% Thresholding
thresholding_mode = p.Results.Thresholding;

% Postprocessing
postprocessing_mode = p.Results.Postprocessing;
samplerate = p.Results.samplerate;
MinimumInaktivityLength = p.Results.MinimumInaktivityLength;
MinimumAktivityLength = p.Results.MinimumAktivityLength;

%% Method
% 1. Wavelet transform
output_signal = wavelet_transform(input_signal,...
                                  'WaveletShape', wavelet_shape,...
                                  'Levels', [4, 5, 6],...
                                  'TransformMode', 'phasefree',...
                                  'OutputMode', 'sum_of_squares');
% 2. Low-pass filter
output_signal = lowpass_filtering(output_signal,...
                                  'FilterType', 'iir',...
                                  'HalfPowerFrequency', f_g,...
                                  'FilterOrder', 2);
% 3. Thresholding
switch thresholding_mode
    case 'on'
        output_signal = compare_with_threshold(output_signal, k);
    case 'off'
        postprocessing_mode = 'off';
    otherwise
end
% 4. Postprocessing
switch postprocessing_mode
    case 'on'
        output_signal = postprocessing(output_signal,...
                                      'MinimumInaktivityLength', MinimumInaktivityLength,...
                                      'MinimumAktivityLength', MinimumAktivityLength,...
                                      'SamplingRate', samplerate);
    case 'off'
    otherwise
end