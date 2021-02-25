function stepfunction = compare_with_threshold(input_signal, threshold_factor)
% -------------------------------------------------------
%
%    COMPARE_WITH_THRESHOLD - Function for the comparison of a signal of
%                             any kind with a threshold to obtain a
%                             stepfunction
%
%    Ver. 1.0.0
%
%    Created:       Alexander Kramlich        (17.07.2015)
%    Last modified: Alexander Kramlich        (17.07.2015)
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
%   COMPARE_WITH_THRESHOLD(input_signal, threshold_factor)
%   COMPARE_WITH_THRESHOLD(input_signal, threshold_factor,...
%                          'parameterName', 'parameterValue')
%
%                           --- Input Arguments ---
%   input_signal:
%       Single signal (e.g. NLEO of an ECG) or several signals in a matrix
%       of dimension MxN, where M is the total number of signals and N the
%       length of one single signal. The given function operates along the
%       second dimension of the matrix.
%
%   threshold_factor:
%       An more or less arbitrary factor, which is multiplied with the
%       standard deviation of the signal to obtain the threshold. 
%
%                           --- Output Arguments ---
%   stepfunction:
%       Stepfunction is a logic array with a length of the original signal.
%       1 (true) indicates the presence of activity, whereas 0 (false)
%       stands for inactivity.
%       The value of the stepfunction is set to 1 (true), if the input
%       signal (e.g. NLEO of an ECG) is above threshold. Otherwise the
%       value of the stepfunction is 0 (false).
%

%% Validation of input
p = inputParser;

addRequired(p, 'input_signal',...
            @(x) validateattributes(x, {'numeric'}, {'nonempty', '2d'}));
addRequired(p, 'threshold_factor',...
            @(x) validateattributes(x, {'numeric'}, {'nonempty', 'scalar', '>', 0}));

parse(p, input_signal, threshold_factor)

%% Function's local variables
% Number of samples in one single signal
N = size(p.Results.input_signal, 2);

%% Comparison with a threshold
% Threshold: T = k * sigma
%   k     - sort of arbitrary factor
%   sigma - standard deviation of the input signal

k = threshold_factor;
sigma = std(input_signal, 0, 2);    % 0 - no weighting,...
                                    % 2 - computation along second dimension 
T = k * sigma;

stepfunction = input_signal > T*ones(1, N);