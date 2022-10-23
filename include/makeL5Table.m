function [L5IcodesTable, L5QcodesTable] = makeL5Table(settings)
%Function generates L5 codes for all 37 satellites based on the settings
%provided in the structure "settings". The codes are digitized at the
%sampling frequency specified in the settings structure.

%L5codesTable = makeL5Table(settings)
%
%   Inputs:
%       settings        - receiver settings
%   Outputs:
%       L5codesTable    - an array of arrays (matrix) containing L5 codes
%                       for all satellite PRN-s


codes_L5I = load ('C:\Users\Sahil\Desktop\Sterenn\Matlab\SDR\include\codes_L5I.mat');
codes_L5Q = load ('C:\Users\Sahil\Desktop\Sterenn\Matlab\SDR\include\codes_L5Q.mat');

%--- Find number of samples per spreading code ----------------------------
samplesPerCode = round(settings.samplingFreq / ...
                           (settings.codeFreqBasis / settings.codeLength));

%--- Prepare the output matrix to speed up function -----------------------
L5IcodesTable = zeros(32, samplesPerCode);
L5QcodesTable = zeros(32, samplesPerCode); 

%--- Find time constants --------------------------------------------------
ts = 1/settings.samplingFreq;   % Sampling period in sec
tc = 1/settings.codeFreqBasis;  % L5 chip period in sec
 
%=== For all satellite PRN-s ...
for PRN = 1:32
    %--- Generate L5 code for given PRN -----------------------------------
    L5Icode = codes_L5I.codes_L5I(:,PRN);
    L5Qcode = codes_L5Q.codes_L5Q(:,PRN);

    %=== Digitizing =======================================================
    
    %--- Make index array to read L5 code values -------------------------
    % The length of the index array depends on the sampling frequency -
    % number of samples per millisecond (because one L5 code period is one
    % millisecond).
    codeValueIndex = ceil((ts/tc) * (1:samplesPerCode));
    
    %--- Correct the last index (due to number rounding issues) -----------
    codeValueIndex(end) = 10230;
    
    %--- Make the digitized version of the L5 code -----------------------
    % The "upsampled" code is made by selecting values form the L5 code
    % chip array (L5code) for the time instances of each sample.
    L5IcodesTable(PRN, :) = L5Icode(codeValueIndex);
    L5QcodesTable(PRN, :) = L5Qcode(codeValueIndex);
    
%     for k=0:19
%         L5IcodesTable(PRN, k*samplesPerCode+1:(k+1)*samplesPerCode) = L5Icode(codeValueIndex);
%         L5QcodesTable(PRN, k*samplesPerCode+1:(k+1)*samplesPerCode) = L5Qcode(codeValueIndex);
%     end % for k=0:19
        
end % for PRN = 1:37
