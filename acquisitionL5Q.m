function acqResults = acquisitionL5Q(longSignal, settings)
%Function performs cold start acquisition on the collected "data". It
%searches for GPS signals of all satellites, which are listed in field
%"acqSatelliteList" in the settings structure. Function saves code phase
%and frequency of the detected signals in the "acqResults" structure.
%
%acqResults = acquisition(longSignal, settings)
%
%   Inputs:
%       longSignal    - 11 ms of raw signal from the front-end 
%       settings      - Receiver settings. Provides information about
%                       sampling and intermediate frequencies and other
%                       parameters including the list of the satellites to
%                       be acquired.
%   Outputs:
%       acqResults    - Function saves code phases and frequencies of the 
%                       detected signals in the "acqResults" structure. The
%                       field "carrFreq" is set to 0 if the signal is not
%                       detected for the given PRN number. 
 
%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
% 
% Copyright (C) Darius Plausinaitis and Dennis M. Akos
% Written by Darius Plausinaitis and Dennis M. Akos
% Based on Peter Rinder and Nicolaj Bertelsen
%--------------------------------------------------------------------------
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
%USA.
%--------------------------------------------------------------------------

%CVS record:
%$Id: acquisition.m,v 1.1.2.12 2006/08/14 12:08:03 dpl Exp $

%% Initialization =========================================================

% Find number of samples per spreading code
samplesPerCode = round(settings.samplingFreq / ...
                        (settings.codeFreqBasis / settings.codeLength));

secondaryCodePeriod = 20;

% Create two 1msec vectors of data to correlate with and one with zero DC
% signal1 = longSignal(1 : samplesPerCode);
% signal2 = longSignal(2*samplesPerCode+1 : 3*samplesPerCode);
% signal3 = longSignal(3*samplesPerCode+1 : 4*samplesPerCode);
signal1 = longSignal(1 : secondaryCodePeriod*samplesPerCode);

signal0DC = longSignal - mean(longSignal); %use to make the power ratio after acquisition

% Find sampling period
ts = 1 / settings.samplingFreq;

% Find phase points of the local carrier wave 
phasePoints_total = (0 : (secondaryCodePeriod*samplesPerCode-1)) * 2 * pi * ts;

% Number of the frequency bins for the given acquisition band (50Hz steps)
numberOfFrqBins = round(settings.acqSearchBand * 2) + 1;

% Generate all L5I codes and sample them according to the sampling freq.
[~,L5QCodesTable] = makeL5Table(settings);
% L5CodesTable = GNSScodegen(settings.acqSatelliteList, 'L5I', 0);
secondaryCodeQ = [-1 -1 -1 -1 -1 1 -1 -1 1 1 -1 1 -1 1 -1 -1 1 1 1 -1];

for i =1:secondaryCodePeriod
    L5QCodesTableTot(:,(i-1)*samplesPerCode+1:i*samplesPerCode)=secondaryCodeQ(i)*L5QCodesTable;
end

NHcode = makeL5QTableSecondary(settings);

%--- Initialize arrays to speed up the code -------------------------------
 % Search results of all frequency bins and code shifts (for one satellite)
results     = zeros(numberOfFrqBins, samplesPerCode);
resultsNH   = zeros(numberOfFrqBins, 50000);
% resultsAbs  = zeros(numberOfFrqBins, samplesPerCode);
resultsAbs = zeros(numberOfFrqBins, samplesPerCode*secondaryCodePeriod);


% Carrier frequencies of the frequency bins
frqBins     = zeros(1, numberOfFrqBins);


%--- Initialize acqResults ------------------------------------------------
% Carrier frequencies of detected signals
acqResults.carrFreq     = zeros(1, settings.acqSatelliteList(end));
% C/A code phases of detected signals
acqResults.codePhase    = zeros(1, settings.acqSatelliteList(end));
% Correlation peak ratios of the detected signals
acqResults.peakMetric   = zeros(1, settings.acqSatelliteList(end));

fprintf('(');

% Perform search for all listed PRN numbers ...
for PRN = settings.acqSatelliteList

    fd = [];
    cp = [];
    IQ = [];
    I = [];
%% Correlate signals ======================================================
% %____ Uncomment this part for the sum acq ________________________________
%     %--- Perform DFT of C/A code ------------------------------------------
%     
%     
%     for p = 0:1:49
%         signal1 = longSignal(p*samplesPerCode+1 : (p+1)*samplesPerCode);
%         phasePoints =  phasePoints_total(p*samplesPerCode+1 : (p+1)*samplesPerCode);
%         L5QCodeFreqDom = conj(fft(L5QCodesTable(PRN, :)));
% %         L5QCodeFreqDom = conj(fft(L5QCodesTable(PRN, :).*NHcode(PRN,p*samplesPerCode+1 : (p+1)*samplesPerCode)));
% 
%         %--- Make the correlation for whole frequency band (for all freq. bins)
%         for frqBinIndex = 1:numberOfFrqBins
% 
%             %--- Generate carrier wave frequency grid (0.5kHz step) -----------
%             frqBins(frqBinIndex) = settings.IF - ...
%                                    (settings.acqSearchBand/2) * 1000 + ...
%                                    0.5e3 * (frqBinIndex - 1);
%  
%             %--- Generate local sine and cosine -------------------------------
%             sinCarr = sin(frqBins(frqBinIndex) * phasePoints);
%             cosCarr = cos(frqBins(frqBinIndex) * phasePoints);
% %             sigCarr = exp(1j*frqBins(frqBinIndex) * phasePoints);
% 
%             %--- "Remove carrier" from the signal -----------------------------
%             I1      = sinCarr .* signal1;
%             Q1      = cosCarr .* signal1;
% %             I1 = real(sigCarr .* signal1);
% %             Q1 = imag(sigCarr .* signal1);
%             
% %              %--- First summation of the 20 reccord of 1ms -----------------
% %               Q1 = 0;
% %               I1 = 0;
% %               for p = 0:19
% %                   signalp = longSignal(p*samplesPerCode+1 : (p+1)*samplesPerCode);
% %                   phasePoints =  phasePoints_total(p*samplesPerCode+1 : (p+1)*samplesPerCode);
% %                   NHcodep = NHcode(PRN,p*samplesPerCode+1 : (p+1)*samplesPerCode);
% %                   sinCarr = sin(frqBins(frqBinIndex) * phasePoints) .* NHcodep;
% %                   cosCarr = cos(frqBins(frqBinIndex) * phasePoints) .* NHcodep;
% %                   
% %                   Ip = sinCarr .* signalp;
% %                   Qp = cosCarr .* signalp;
% %                   I1 = I1 + Ip;
% %                   Q1 = Q1 + Qp;
% %               end %for p = 0:19
% 
%              %--- Convert the baseband signal to frequency domain --------------
%             IQfreqDom1 = fft(I1 + 1j*Q1);
% 
%             %--- Multiplication in the frequency domain (correlation in time
%             %domain)
%             convCodeIQ = IQfreqDom1 .* L5QCodeFreqDom;
% 
%             %--- Perform inverse DFT and store correlation results ------------
%             acqRes11 = ifft(convCodeIQ);
% 
%             %--- Check which msec had the greater power and save that, will
%             %"blend" 1st and 2nd msec but will correct data bit issues
% 
%             results(frqBinIndex, :) = acqRes11;
% %             [peakSizeIQ, frequencyBinIndexIQ] = max(max(results, [], 2));
% % 
% %             %--- Find code phase of the same correlation peak ---------------------
% %             [peakSizeIQ, codePhaseIQ] = max(max(results));
% 
%         end % frqBinIndex = 1:numberOfFrqBins
%         eval(['CM' num2str(p,'%02d') '=results;']);
%         eval(['CAF' num2str(p,'%02d') '=abs(results);']);
% %         eval(['IQ' num2str(p,'%02d') '=results(frequencyBinIndex, codePhase);']);
% %         IQ = [IQ; results(frequencyBinIndexIQ, codePhaseIQ)];
% %         fd = [fd;frequencyBinIndexIQ];
% %         cp = [cp;codePhaseIQ];
%         
%      end %for p = 0:1:19
%      
% %      I = sign(imag(IQ));
% %      CorrNH = xcorr(I,secondaryCodeQ');
%     CCAF = CM00+CM01+CM02+CM03+CM04+CM05+CM06+CM07+CM08+CM09+CM10+CM11+CM12+CM13+CM14+CM15+CM16+CM17+CM18+CM19; 
% %     CCAF = CM00*secondaryCodeQ(1)+CM01*secondaryCodeQ(2)+CM02*secondaryCodeQ(3)+CM03*secondaryCodeQ(4)+CM04*secondaryCodeQ(5)...
% %         +CM05*secondaryCodeQ(6)+CM06*secondaryCodeQ(7)+CM07*secondaryCodeQ(8)+CM08*secondaryCodeQ(9)+CM09*secondaryCodeQ(10)...
% %         +CM10*secondaryCodeQ(11)+CM11*secondaryCodeQ(12)+CM12*secondaryCodeQ(13)+CM13*secondaryCodeQ(14)+CM14*secondaryCodeQ(15)...
% %         +CM15*secondaryCodeQ(16)+CM16*secondaryCodeQ(17)+CM17*secondaryCodeQ(18)+CM18*secondaryCodeQ(19)+CM19*secondaryCodeQ(20);
%     resultsAbs = abs(CCAF).^2;
% 
% %     resultsAbs = abs(results).^2;
    
%_____ UnComment this part to do the 20ms acq _____________________________

 %--- Perform DFT of CL5 code ------------------------------------------
    L5QCodeFreqDom = conj(fft(L5QCodesTableTot(PRN, :)));

        %--- Make the correlation for whole frequency band (for all freq. bins)
        for frqBinIndex = 1:numberOfFrqBins

            %--- Generate carrier wave frequency grid (0.5kHz step) -----------
            frqBins(frqBinIndex) = settings.IF - ...
                                   (settings.acqSearchBand/2) * 1000 + ...
                                   0.5e3 * (frqBinIndex - 1);

            %--- Generate local sine and cosine -------------------------------
            sinCarr = sin(frqBins(frqBinIndex) * phasePoints_total);
            cosCarr = cos(frqBins(frqBinIndex) * phasePoints_total);
%             sigCarr = exp(1j*frqBins(frqBinIndex) * phasePoints_total);

            %--- "Remove carrier" from the signal -----------------------------
            I1      = sinCarr .* signal1;
            Q1      = cosCarr .* signal1;
%             I1 = real(sigCarr .* signal1);
%             Q1 = imag(sigCarr .* signal1);

             %--- Convert the baseband signal to frequency domain --------------
            IQfreqDom1 = fft(I1 + 1j*Q1);

            %--- Multiplication in the frequency domain (correlation in time
            %domain)
            convCodeIQ = IQfreqDom1 .* L5QCodeFreqDom;

            %--- Perform inverse DFT and store correlation results ------------
            acqRes11 = abs(ifft(convCodeIQ)) .^2;

            %--- Check which msec had the greater power and save that, will
            %"blend" 1st and 2nd msec but will correct data bit issues

            resultsAbs(frqBinIndex, :) = acqRes11;
    
        end %for frqBinIndex = 1:


%% Look for correlation peaks in the results ==============================
    % Find the highest peak and compare it to the second highest peak
    % The second peak is chosen not closer than 1 chip to the highest peak

    %--- Find the correlation peak and the carrier frequency --------------
    [peakSize, frequencyBinIndex] = max(max(resultsAbs, [], 2));

    %--- Find code phase of the same correlation peak ---------------------
    [peakSize, codePhase] = max(max(resultsAbs));

    %--- Find 1 chip wide C/A code phase exclude range around the peak ----
    samplesPerCodeChip  = round(settings.samplingFreq / settings.codeFreqBasis);
    excludeRangeIndex11 = codePhase - samplesPerCodeChip;
    excludeRangeIndex12 = codePhase + samplesPerCodeChip;

    %--- Correct C/A code phase exclude range if the range includes array
    %boundaries
    if excludeRangeIndex11 < 2
        codePhaseRange = excludeRangeIndex12 : ...
                         (samplesPerCode + excludeRangeIndex11);
                         
    elseif excludeRangeIndex12 >= samplesPerCode
        codePhaseRange = (excludeRangeIndex12 - samplesPerCode) : ...
                         excludeRangeIndex11;
    else
        codePhaseRange = [1:excludeRangeIndex11, ...
                          excludeRangeIndex12 : samplesPerCode];
    end

    %--- Find the second highest correlation peak in the same freq. bin ---
    [secondPeakSize, codePhase2] = max(resultsAbs(frequencyBinIndex, codePhaseRange));
    if (peakSize/secondPeakSize) < settings.acqThreshold
        excludeRangeIndex21 = codePhase2 - samplesPerCodeChip;
        excludeRangeIndex22 = codePhase2 + samplesPerCodeChip;

        if codePhase2 >= codePhase
            codePhaseRange = [codePhaseRange(1:excludeRangeIndex11),codePhaseRange(excludeRangeIndex22:end)];          
        elseif codePhase2 <= codePhase
            codePhaseRange = [codePhaseRange(1:excludeRangeIndex21),codePhaseRange(excludeRangeIndex12:end)];                 
        end %if codePhase2 >= codePhase
    
    end %if (peakSize/secondPeakSize) < settings.acqThreshold
    
    thirdPeakSize  = max(resultsAbs(frequencyBinIndex, codePhaseRange));
    %--- Store result -----------------------------------------------------
    acqResults.peakMetric(PRN) = peakSize/thirdPeakSize;
    
    % If the result is above threshold, then there is a signal ...
    if (peakSize/thirdPeakSize) > settings.acqThreshold

%% Fine resolution frequency search =======================================
%%--- Plot FFTs of the signal acquistions if plotFFTs is high -------------  
        figure(PRN)
        subplot(1,2,1);
        plot(resultsAbs(frequencyBinIndex,:));
            title ('Acquisition Results, Code', 'FontSize',25);
            xlabel('Code delay (samples)', 'FontSize',20);
            ylabel('Correlation Absolute Value', 'FontSize',20);
        subplot(1,2,2);
        plot(resultsAbs(:,codePhase));
            title ('Acquisition results, Frequency', 'FontSize',25);
            xlabel('Frequency offset (bins)', 'FontSize',20);
            ylabel('Correlation Absolute Value', 'FontSize',20);          

        %--- Indicate PRN number of the detected signal -------------------
        fprintf('%02d ', PRN);
        
        %--- Generate 10msec long C/A codes sequence for given PRN --------
%         PRNL = PRN(ones(1,1));
        L5Code = GNSScodegen(PRN,'L5Q',0);
        
        codeValueIndex = floor((ts * (1:2*samplesPerCode)) / ...
                               (1/settings.codeFreqBasis));
                           
        longL5Code = L5Code((rem(codeValueIndex, 10230) + 1));
    
        %--- Remove C/A code modulation from the original signal ----------
        % (Using detected C/A code phase)
        xCarrier = ...
            signal0DC(codePhase:(codePhase + 2*samplesPerCode-1)) ...
            .* longL5Code;
        
%         %% Do the correalation with the NH code ===========================
%         NHCodeFreqDom = conj(fft(NHcode(PRN, :)));
%         for frqBinIndex = 1:numberOfFrqBins
%             %--- Generate carrier wave frequency grid (0.5kHz step) -----------
%             frqBins(frqBinIndex) = settings.IF - ...
%                                    (settings.acqSearchBand/2) * 1000 + ...
%                                    0.5e3 * (frqBinIndex - 1);
%             
%             NHcorr = xCarrier .* NHCodeFreqDom(1:50000);
%             Convprod = ifft(NHcorr);
%             resultsNH(frqBinIndex, :) = Convprod;
%         end %frqBinIndex = 1:numberOfFrqBins
%          Absres = abs(resultsNH).^2;
%         
%          %--- Find the correlation peak and the carrier frequency --------------
%         [NHpeakSize, NHfrequencyBinIndex] = max(max(Absres, [], 2));
% 
%         %--- Find code phase of the same correlation peak ---------------------
%         [NHpeakSize, NHcodePhase] = max(max(Absres));
%         
%         %% end NH correlation =============================================
         
        %--- Find the next highest power of two and increase by 8x --------
        fftNumPts = 8*(2^(nextpow2(length(xCarrier))));
        
        %--- Compute the magnitude of the FFT, find maximum and the
        %associated carrier frequency 
        fftxc = abs(fft(xCarrier, fftNumPts)); 
        
        
        uniqFftPts = ceil((fftNumPts + 1) / 2);
        [fftMax, fftMaxIndex] = max(fftxc);
        fftFreqBins = (0 : uniqFftPts-1) * settings.samplingFreq/fftNumPts;
        
        offsetposition=0;
         % find the average value among FFT peaks
        
        % define a searching area to find all FFT peaks
        searchingPoint =150;
        leftBoundary = fftMaxIndex - 150 ;
        rightBoundary = fftMaxIndex + 150;
        if leftBoundary <2
            searchingPoint =  150+ leftBoundary;
            leftBoundary =2;
        end
        if rightBoundary > (fftNumPts -1)
            rightBoundary = fftNumPts -1;
        end    
        searchArea=fftxc(leftBoundary:rightBoundary);% FFT result in searching area
        
        % find the position FFT peaks
        position = find (searchArea>fftMax*0.75);% 0.75 is a threshold which defines the minimum value of a "peak"
        midposition = ceil((max(position)+min(position))/2); %find the average frequnecy
        offsetposition = midposition - searchingPoint;
     
        if (fftMaxIndex > uniqFftPts) %and should validate using complex data
            if (rem(fftNumPts,2)==0)  %even number of points, so DC and Fs/2 computed
                fftFreqBinsRev=-fftFreqBins((uniqFftPts-1):-1:2);
                [fftMax, fftMaxIndex] = max(fftxc((uniqFftPts+1):length(fftxc)));
                acqResults.carrFreq(PRN)  = -fftFreqBinsRev(fftMaxIndex+offsetposition);
            else  %odd points so only DC is not included
                fftFreqBinsRev=-fftFreqBins((uniqFftPts):-1:2);
                [fftMax, fftMaxIndex] = max(fftxc((uniqFftPts+1):length(fftxc)));
                acqResults.carrFreq(PRN)  = fftFreqBinsRev(fftMaxIndex);
            end
        else
            acqResults.carrFreq(PRN)  = (-1)^(settings.fileType-1)*fftFreqBins(fftMaxIndex+offsetposition);
        end
        
%         % Normal search
%         if (fftMaxIndex > uniqFftPts) %and should validate using complex data
%             if (rem(fftNumPts,2)==0)  %even number of points, so DC and Fs/2 computed
%                 fftFreqBinsRev=-fftFreqBins((uniqFftPts-1):-1:2);
%                 [fftMax, fftMaxIndex] = max(fftxc((uniqFftPts+1):length(fftxc)));
%                 acqResults.carrFreq(PRN)  = -fftFreqBinsRev(fftMaxIndex);
%             else  %odd points so only DC is not included
%                 fftFreqBinsRev=-fftFreqBins((uniqFftPts):-1:2);
%                 [fftMax, fftMaxIndex] = max(fftxc((uniqFftPts+1):length(fftxc)));
%                 acqResults.carrFreq(PRN)  = fftFreqBinsRev(fftMaxIndex);
%             end
%         else
%             acqResults.carrFreq(PRN)  = (-1)^(settings.fileType-1)*fftFreqBins(fftMaxIndex);
%         end
        
        %--- Save properties of the detected satellite signal -------------
        acqResults.carrFreq(PRN)  = fftFreqBins(fftMaxIndex);
        acqResults.codePhase(PRN) = codePhase;
%         acqResults.codePhase2(PRN) = codePhase2;
    
    else
        %--- No signal with this PRN --------------------------------------
        fprintf('. ');
    end   % if (peakSize/secondPeakSize) > settings.acqThreshold
    
end    % for PRN = satelliteList

%=== Acquisition is over ==================================================
fprintf(')\n');
