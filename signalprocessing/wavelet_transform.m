function output_signal = wavelet_transform(input_signal, varargin)
% -------------------------------------------------------
%
%    WAVELET_TRANSFORM - Wavelet transform.
%
%    Ver. 2.1.0
%
%    Created:       Alexander Kramlich        (02.05.2015)
%    Last modified: Alexander Kramlich        (12.10.2015)
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
%   Description
%
%   WAVELET_TRANSFORM(signal)
%   WAVELET_TRANSFORM(signal, 'parameterName', 'parameterValue')
%
%                           --- Input Arguments ---
%   signal:
%       Single signal (e.g. EGM, ECG, etc.) (or several signals in a matrix
%       of dimension MxN, where M is the total number of signals and N the
%       length of one single signal).
%   
%   Name-Value Pair Arguments:
%       'WaveletShape'
%           ---
%           'db1' / 'db2' / 'db3' / 'db4' / 'db5' /
%           'coif1' / 'coif2' / 'coif3' / 'coif4' (default) / 'coif5' /
%           'bior1.1' / 'bior1.3' / 'bior1.5' /
%           'bior2.2' / 'bior2.4' / 'bior2.6' / 'bior2.8' /
%           'bior3.1' / 'bior3.3' / 'bior3.5' / 'bior3.7' / 'bior3.9' /
%           'bior4.4' / 'bior5.5' / 'bior6.8'
%           ---
%           Shape of the wavelet used for the wavelet transform.
%       'Levels'
%           Levels of the wavelet transform used for the computation of the
%           output.
%           Default value: [4,5,6]
%       'TransformMode'
%           --- 'phasefree' (default) / 'regular' ---
%           Parameter to switch between regular and phase free
%           wavelet transform technique.
%       'OutputMode'
%           --- 'sum_of_squares' (default) / 'sum' / 'product' ---
%           Parameter to choose the method of how the levels are put
%           together.
%

%% Validation of input
p = inputParser;
defaultWaveletShape = 'coif4';
defaultLevels = [4, 5, 6];
defaultTransformMode = 'phasefree';
defaultOutputMode = 'sum_of_squares';

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
expectedTransformModes = {'phasefree', 'regular'};
expectedOutputModes = {'sum_of_squares', 'sum', 'product'};
                      
addRequired(p, 'signal', @(x) validateattributes(x, {'numeric'}, {'nonempty', '2d'}));

addParameter(p, 'WaveletShape',...
                defaultWaveletShape,...
                @(x) any(validatestring(x, expectedWaveletShapes)));
addParameter(p, 'Levels',...
                defaultLevels,...
                @(x) validateattributes(x, {'numeric'}, {'nonempty', 'row', '>', 0, '<=', 10}));
addParameter(p, 'TransformMode',...
                defaultTransformMode,...
                @(x) any(validatestring(x, expectedTransformModes)));
addParameter(p, 'OutputMode',...
                defaultOutputMode,...
                @(x) any(validatestring(x, expectedOutputModes)));
            
parse(p, input_signal, varargin{:})

%% Function's local variables
[M, N] = size(input_signal);

% Wavelet shape
wavelet_shape = p.Results.WaveletShape;

% Levels
levels = p.Results.Levels;

% Mode
transform_mode = p.Results.TransformMode;
output_mode = p.Results.OutputMode;

%% Wavelet transform
% Extension of the EMG signal to the length equal to the next power of 2
N_ext = 2^nextpow2(N);
n_right = floor((N_ext-N)/2);
n_left = N_ext-N-n_right;

signal_ext = wextend('ac', 'sp0', input_signal, n_right, 'r');
signal_ext = wextend('ac', 'sp0', signal_ext, n_left, 'l');

t = N_ext:-1:1;

% Discrete Wavelet Transformation
%
% Sampling frequency: 2048 Hz
%
%   --- Detail coefficients ---
% -------------------------------
%  Level     Frequency range (Hz)
% -------------------------------
%      1            512-1024
%      2             256-512
%      3             128-256
%      4              64-128
%      5               32-64
%      6               16-32
%      7                8-16
%    ...                 ...
% ------------------------------
K = max(levels);

D_forward = zeros(K, M, N); 
D_backward = zeros(K, M, N);
D_output = zeros(K, M, N);

switch transform_mode
    case 'phasefree'
        for m=1:M
            % Forward transform
            [~,swd] = swt(signal_ext(m,:), K, wavelet_shape);
            D_forward(:,m,:) = swd(:, n_left+1:end-n_right); 

            % Backward transform
            [~,swd] = swt(fliplr(signal_ext(m,:)), K, wavelet_shape);
            swd = fliplr(swd);
            D_backward(:,m,:) = swd(:, n_left+1:end-n_right);

            % Addition of forward and backward transforms
            D_output = D_forward + D_backward;
        end
    case 'regular'
        for m=1:M
            % Forward transform
            [~,swd] = swt(signal_ext(m,:), K, wavelet_shape);
            D_output(:,m,:) = swd(:, n_left+1:end-n_right); 
        end
    otherwise
end

switch output_mode
    case 'sum_of_squares'
        % Sum of squares
        output_signal = zeros(M,N);
        for k=levels
            output_signal = output_signal + reshape(D_output(k,:,:), [M,N]).^2;
        end
    case 'sum'
        % Sum
        output_signal = zeros(M,N);
        for k=levels
            output_signal = output_signal + reshape(D_output(k,:,:), [M,N]);
        end
    case 'product'
        % Product
        % Note: Product of detail coefficients in time domain corresponds
        % to a convolution in frequency domain, which is not intuitive and
        % the outcome is difficult to predict. However the calculation of a
        % product may reveal signal information, which at a first
        % appearence are hidden.
        output_signal = ones(M,N);
        for k=levels
            output_signal = output_signal .* reshape(D_output(k,:,:), [M,N]);
        end
    otherwise
end