function acqResults = postacquisition(longSignal, settings)

%% Initialization =========================================================
binsize = 10;
prn = settings.acqSatelliteList;

% carrFreq = 5001613.91198644;
% load('C:\Users\Sahil\Desktop\GNSS\Updated_Code\cleanStatic\carrfreq_cleanstatic.mat')
% load('C:\Users\Sahil\Desktop\GNSS\Updated_Code\Wengxiang data\carr_freq_wengxiang.mat')
% load('carrFreq.mat');
% Find number of samples per spreading code
samplesPerCode = round(settings.samplingFreq / ...
                        (settings.codeFreqBasis / settings.codeLength));
fd = [];
cp = [];
ps = [];
raw_data = [];
f = [];
remCarrPhase = 0;

% Find phase points of the local carrier wave 
ts = 1 / settings.samplingFreq;
phasePoints_total = (0 : (20*samplesPerCode-1)) * 2 * pi * ts;

signal0DC = longSignal - mean(longSignal);

%--- Initialize acqResults ------------------------------------------------
% Carrier frequencies of detected signals
acqResults.carrFreq     = zeros(1, settings.acqSatelliteList(end));
% C/A code phases of detected signals
acqResults.codePhase    = zeros(1, settings.acqSatelliteList(end));
% Correlation peak ratios of the detected signals
acqResults.peakMetric   = zeros(1, settings.acqSatelliteList(end));

fprintf('(');

for PRN = settings.acqSatelliteList

 for p = 0:1:19
%     for p = 0
signal1 = longSignal(p*samplesPerCode+1 : (p+1)*samplesPerCode);
raw_data = [raw_data;signal1];

% Find sampling period
phasePoints = phasePoints_total(p*samplesPerCode+1 : (p+1)*samplesPerCode);

% ts = 1 / settings.samplingFreq;
% % Find phase points of the local carrier wave 
% phasePoints = (0 : (samplesPerCode-1)) * 2 * pi * ts;


% Number of the frequency bins for the given acquisition band (500Hz steps)
numberOfFrqBins = round(settings.acqSearchBand * 2) + 1;

% Generate all C/A codes and sample them according to the sampling freq.
[~, L5QCodesTable] = makeL5Table(settings);


%--- Initialize arrays to speed up the code -------------------------------
% Search results of all frequency bins and code shifts (for one satellite)
results     = zeros(numberOfFrqBins, samplesPerCode);

CN1     = zeros(numberOfFrqBins, samplesPerCode);
% CN     = zeros(numberOfFrqBins, samplesPerCode);
% CN     = zeros(numberOfFrqBins, 2*samplesPerCode-1);


% Carrier frequencies of the frequency bins
frqBins     = zeros(1, numberOfFrqBins);
% prn = [3,4,6,7,10,13,16,19,23];

% prn = find(abs(acqResults.carrFreq) > 0);
% Perform search for all listed PRN numbers ...
% for s = 1:length(prn)
%     PRN = prn(s);

%% Correlate signals ======================================================   
    %--- Perform DFT of C/A code ------------------------------------------
    caCodeFreqDom = conj(fft(L5QCodesTable(PRN, :)));

    %--- Make the correlation for whole frequency band (for all freq. bins)
    for frqBinIndex = 1:numberOfFrqBins

        %--- Generate carrier wave frequency grid (0.5kHz step) -----------
%         frqBins(frqBinIndex) = acqResults.carrFreq(PRN);

        frqBins(frqBinIndex) = settings.IF - ...
                               (settings.acqSearchBand/2) + ...
                               15 * (frqBinIndex - 1);
%         frqBins(frqBinIndex) = settings.IF - 5*binsize* settings.acqSearchBand2 + binsize* (frqBinIndex +1);
%         frqBins  = -1860;                      
%           frqBins(frqBinIndex) = (carrFreq(PRN)-1750) + binsize* (frqBinIndex +1);
        %--- Generate local sine and cosine -------------------------------
        sigCarr = exp(1i*frqBins(frqBinIndex) * phasePoints);
%         sigCarr1 = exp(i*frqBins(frqBinIndex) * phasePoints1);
        
        %--- "Remove carrier" from the signal -----------------------------
        I1      = real(sigCarr .* signal1);
        Q1      = imag(sigCarr .* signal1);

        %--- Convert the baseband signal to frequency domain --------------
        IQfreqDom1 = fft(I1 + 1j*Q1);

        %--- Multiplication in the frequency domain (correlation in time
        %domain)
        convCodeIQ1 = IQfreqDom1 .* caCodeFreqDom;
%         convCodeIQ2 = IQfreqDom2 .* caCodeFreqDom;
%         correlation_I = conv(I1,caCodesTable(PRN, :),"same");
%         correlation_Q = conv(Q1,caCodesTable(PRN, :),"same");
%         correlation = correlation_I + j*correlation_Q;

%         correlation_I1 = xcorr(I1,caCodesTable(PRN, :));
%         correlation_Q1 = xcorr(Q1,caCodesTable(PRN, :));
%         correlation1 = correlation_I1 + j*correlation_Q1;
% 
%         conv_I1 = conv(I1,fliplr(caCodesTable(PRN, :)));
%         conv_Q1 = conv(Q1,fliplr(caCodesTable(PRN, :)));
%         convolution1 = conv_I1 + j*conv_Q1;
%         convCodeIQ11 = IQfreqDom11 .* caCodeFreqDom;
%         convCodeIQ22 = IQfreqDom22 .* caCodeFreqDom;

        %--- Perform inverse DFT and store correlation results ------------
        r1 = ifft(convCodeIQ1);
%         r2 = ifft(convCodeIQ2);
%         i1 = imag(ifft(convCodeIQ1));
%         i2 = imag(ifft(convCodeIQ2));
%         p1 = angle(ifft(convCodeIQ1));
%         p2 = angle(ifft(convCodeIQ2));
%         p1 = unwrap(p1,[],2);
%         p2 = unwrap(p2,[],2);
        acqRes1 = abs(ifft(convCodeIQ1)) .^ 2;
%         acqRes2 = abs(ifft(convCodeIQ2)) .^ 2;
        
%         acqRes11 = abs(ifft(convCodeIQ11)) .^ 2;
%         acqRes22 = abs(ifft(convCodeIQ22)) .^ 2;
%         acqRes1 = (r1+i1).^2;
%         acqRes2 = (r2+i2).^2;
%         acqRes1 = r1;
%         acqRes2 = r2;
        %--- Check which msec had the greater power and save that, will
        %"blend" 1st and 2nd msec but will correct data bit issues
%         if (max(acqRes1) > max(acqRes2))
            results(frqBinIndex, :) = acqRes1;
%             results1(frqBinIndex, :) = acqRes11;
            CN1(frqBinIndex, :) = r1;
%             CN(frqBinIndex, :) = correlation1;
%             CN(frqBinIndex, :) = correlation1;
%         else
%             results(frqBinIndex, :) = acqRes2;
% %             results1(frqBinIndex, :) = acqRes22;
%              CN1(frqBinIndex, :) = r2;
%         end
        
%         results = abs(results);
    end % frqBinIndex = 1:numberOfFrqBins

%% Look for correlation peaks in the results ==============================
    % Find the highest peak and compare it to the second highest peak
    % The second peak is chosen not closer than 1 chip to the highest peak
%     f = [f;frqBins];
    %--- Find the correlation peak and the carrier frequency --------------
    [peakSize, frequencyBinIndex] = max(max(CN1, [], 2));

    %--- Find code phase of the same correlation peak ---------------------
    [peakSize, codePhase] = max(max(CN1));
    P_max = peakSize;
    %--- Find 1 chip wide C/A code phase exclude range around the peak ----
%     samplesPerCodeChip   = round(settings.samplingFreq / settings.codeFreqBasis);
%     excludeRangeIndex1 = codePhase - samplesPerCodeChip;
%     excludeRangeIndex2 = codePhase + samplesPerCodeChip;
% %     ps = [ps;abs(max(max(CN1)))];
%     %--- Correct C/A code phase exclude range if the range includes array
%     %boundaries
%     if excludeRangeIndex1 < 2
%         codePhaseRange = excludeRangeIndex2 : ...
%                          (samplesPerCode + excludeRangeIndex1);
%                          
%     elseif excludeRangeIndex2 >= samplesPerCode
%         codePhaseRange = (excludeRangeIndex2 - samplesPerCode) : ...
%                          excludeRangeIndex1;
%     else
%         codePhaseRange = [1:excludeRangeIndex1, ...
%                           excludeRangeIndex2 : samplesPerCode];
%     end
    fd = [fd;frequencyBinIndex];
    cp = [cp;codePhase];
    eval(['CM' num2str(p,'%02d') '=CN1;']);
    eval(['CAF' num2str(p,'%02d') '=abs(CN1);']);
%     CP = abs(CN1(frequencyBinIndex,codePhase-150:codePhase+150));
    FD = abs(CN1(:,codePhase));
    CP = abs(CN1(frequencyBinIndex,:));
%     eval(['CP' num2str(p,'%02d') '=CP;']);
    eval(['FD' num2str(p,'%02d') '=FD;']);
    %--- Find the second highest correlation peak in the same freq. bin ---
%     secondPeakSize = max(results(frequencyBinIndex, codePhaseRange));
%     P_noise = mean(CN1(frequencyBinIndex, codePhaseRange));
%     P_var = var(CN1(frequencyBinIndex, codePhaseRange));
%     eval(['P' num2str(p,'%02d') '=P_noise;']);
%     eval(['V' num2str(p,'%02d') '=P_var;']);
    
%     figure(s)
%     CN = abs(CN1);
%     CP = CN(frequencyBinIndex,codePhase-75:codePhase+75);
%     DP = CN(:,codePhase);
%     CN2 = CN1(:,codePhase-75:codePhase+75); 
%     mesh(frqBins,linspace(1,1023,samplesPerCode),results')
% title(['PRN',num2str(PRN)]);
% xlabel('f_d (Hz)');
% ylabel('\tau (chips)');
% zlabel('abs(CCAF)');
 end % p= 0:1:19



% CCAF = CM00+CM01+CM02+CM03+CM04+CM05+CM06+CM07+CM08+CM09+CM10+CM11+CM12+CM13+CM14+CM15+CM16+CM17+CM18+CM19; 
CCAF_p = CM00;
CAF = abs(CCAF_p);
    [peakSize, frequencyBinIndex] = max(max(CCAF_p, [], 2));
    [peakSize, codePhase] = max(max(CCAF_p));
    f = frqBins(frequencyBinIndex);
    phase = angle(CCAF_p(frequencyBinIndex,codePhase));
%     CP = abs(CCAF(frequencyBinIndex,codePhase-150:codePhase+150));
%     FD = abs(CCAF(:,codePhase));
%     
%     frequency = frqBins(frequencyBinIndex);
% 
% F_CP = 0;
% for v = 1:1:20
%     CAF = eval(['CM' num2str(p,'%02d')]);
% C = CAF(fd(v),cp(v)-100:cp(v)+100);
% % C = CAF(fd(v),:);
% F_CP = F_CP + C;
% end
% A = abs(F_CP);
% 
% % CP = CAF(frequencyBinIndex,codePhase-300:codePhase+300);
% fd_bin = frqBins(fd); 
% 
% P = P00+P01+P02+P03+P04+P05+P06+P07+P08+P09+P10+P11+P12+P13+P14+P15+P16+P17+P18+P19;
% V = V00+V01+V02+V03+V04+V05+V06+V07+V08+V09+V10+V11+V12+V13+V14+V15+V16+V17+V18+V19;

%--- Find 1 chip wide C/A code phase exclude range around the peak ----
    samplesPerCodeChip   = round(settings.samplingFreq / settings.codeFreqBasis);
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
    [secondPeakSize, codePhase2] = max(CAF(frequencyBinIndex, codePhaseRange));
    if (peakSize/secondPeakSize) < settings.acqThreshold
        excludeRangeIndex21 = codePhase2 - samplesPerCodeChip;
        excludeRangeIndex22 = codePhase2 + samplesPerCodeChip;

        if codePhase2 >= codePhase
            codePhaseRange = [codePhaseRange(1:excludeRangeIndex11),codePhaseRange(excludeRangeIndex22:end)];          
        elseif codePhase2 <= codePhase
            codePhaseRange = [codePhaseRange(1:excludeRangeIndex21),codePhaseRange(excludeRangeIndex12:end)];                 
        end %if codePhase2 >= codePhase
    
    end %if (peakSize/secondPeakSize) < settings.acqThreshold
    
    thirdPeakSize  = max(CAF(frequencyBinIndex, codePhaseRange));
    %--- Store result -----------------------------------------------------
    acqResults.peakMetric(PRN) = peakSize/thirdPeakSize;
    
    % If the result is above threshold, then there is a signal ...
    if (peakSize/thirdPeakSize) > settings.acqThreshold

%% Fine resolution frequency search =======================================
%%--- Plot FFTs of the signal acquistions if plotFFTs is high ------------  
%         figure(PRN)
%         subplot(1,2,1);
%         plot(CAF(frequencyBinIndex,:));
%             title ('Acquisition Results, Code', 'FontSize',25);
%             xlabel('Code delay (samples)', 'FontSize',20);
%             ylabel('Correlation Absolute Value', 'FontSize',20);
%         subplot(1,2,2);
%         plot(CAF(:,codePhase));
%             title ('Acquisition results, Frequency', 'FontSize',25);
%             xlabel('Frequency offset (bins)', 'FontSize',20);
%             ylabel('Correlation Absolute Value', 'FontSize',20);          


        %--- Indicate PRN number of the detected signal -------------------
        fprintf('%02d ', PRN);
        
        %--- Generate 10msec long C/A codes sequence for given PRN --------
        PRNL = PRN(ones(1,1));
        L5Code = GNSScodegen(PRNL,'L5Q',0);
        
        codeValueIndex = floor((ts * (1:5*samplesPerCode)) / ...
                               (1/settings.codeFreqBasis));
                           
        longL5Code = L5Code((rem(codeValueIndex, 10230) + 1));
    
        %--- Remove C/A code modulation from the original signal ----------
        % (Using detected C/A code phase)
        xCarrier = ...
            signal0DC(codePhase:(codePhase + 5*samplesPerCode-1)) ...
            .* longL5Code;
        
        %--- Find the next highest power of two and increase by 8x --------
        fftNumPts = 8*(2^(nextpow2(length(xCarrier))));
        
        %--- Compute the magnitude of the FFT, find maximum and the
        %associated carrier frequency 
        fftxc = abs(fft(xCarrier, fftNumPts)); 
        
        uniqFftPts = ceil((fftNumPts + 1) / 2);
        [fftMax, fftMaxIndex] = max(fftxc(5 : uniqFftPts-5));
        
        fftFreqBins = (0 : uniqFftPts-1) * settings.samplingFreq/fftNumPts;
        
        %--- Save properties of the detected satellite signal -------------
        acqResults.carrFreq(PRN)  = fftFreqBins(fftMaxIndex);
        acqResults.codePhase1(PRN) = codePhase;
        acqResults.codePhase2(PRN) = codePhase2;
    
    else
        %--- No signal with this PRN --------------------------------------
        fprintf('. ');
    end   % if (peakSize/secondPeakSize) > settings.acqThreshold

 end  % PRN...
