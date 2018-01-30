function [ f_CENS ] = midi_chrama( midiInfo, FS, WS, HS, plotPattern )

%% 每0.1拍一個chrama
%   input  : midiInfo
%   output : chrama feature

 %%%%%%%%%%%%%%%%% 即使不同的拍號，也可以畫圖
 
 if nargin < 2, FS  = 0.01;       end
 if nargin < 3, WS  = 0.2;        end
 if nargin < 4, HS  = WS/2;       end
 if nargin < 5, plotPattern  = 1; end
 
%% testing code
%     clear all; close all; clc;
%     midiInfo    = midi_Preprocess('mz_545_1', 1);   % 需要 midiinfo 1(onset beat) 4(pitch)
%     plotPattern = 1;
%     WS = 0.2;    
%     HS = WS/2;
%     FS = 0.01;
    
%% pitch curve
    
    fs = 1/FS;
    newCurve = eps*ones(12, round(midiInfo(end,6)+midiInfo(end,7))*fs);
    
    for j = 1:size(midiInfo,1)
        Start = round(midiInfo(j,6)*fs); if Start==0; Start=1; end
        newCurve(mod(midiInfo(j,4)-60,12)+1, Start:round((midiInfo(j,6)+midiInfo(j,7))*fs) ) = 1;
    end

    ws        = WS*fs;
    hs        = HS*fs;
    hopNum    = ceil((size(newCurve,2)-ws)/hs);
    zeroPad_C = [ newCurve, eps*ones(12, hopNum*hs+ws-size(newCurve,2))];
    
    seg_num   = hopNum + 1;
    
    for j=1:seg_num
        overlapCurve(:, j) = sum(zeroPad_C(:,hs*(j-1)+1:hs*(j-1)+ws),2); 
        if sum( overlapCurve(:,j)>0.001 ) > 0
            overlapCurve(:, j) = overlapCurve(:, j)/max(overlapCurve(:, j));
        end
    end
    
    %% CENS 
    
    f_chroma_energy_distr  =  overlapCurve;
    parameter.quantSteps   = [40 20 10 5] / 100;
    parameter.quantWeights = [ 1 1 1 1]/4;
    % calculate a CENS feature

    % component-wise quantisation of the normalized chroma vectors
    f_stat_help = zeros(12,seg_num);

    for n=1:length(parameter.quantSteps)
        f_stat_help = f_stat_help + (f_chroma_energy_distr>parameter.quantSteps(n))*parameter.quantWeights(n);
    end

    paramCENS.winLenSmooth = 11;
    paramCENS.downsampSmooth = 5;

    % Temporal smoothing and downsampling
    [f_chroma_energy_stat,~] = smoothDownsampleFeature(f_stat_help,paramCENS);

    parameter.normThresh = 0.0001;
    % last step: normalize each vector with its l^2 norm
    f_CENS = normalizeFeature(f_chroma_energy_stat,2, parameter.normThresh);


    if plotPattern
        figure;
        x1 = ((1:length(overlapCurve)-1)*hs+ws)/fs;
        x2 = ((1:5:length(overlapCurve)-1)*hs+ws)/fs;
        y  = 1:12;
        subplot(2,1,1); imagesc(x1,y,overlapCurve);  xlabel('sec'); ylabel('chrama'); title(['f\_chrama : WS=' num2str(WS) '_s_e_c , HS=' num2str(HS) '_s_e_c']);
        subplot(2,1,2); imagesc(x2,y,f_CENS);  xlabel('sec'); ylabel('chrama');   title('f\_CENS');
    end
    
%     paramSM.smoothLenSM = 20;  % 改小一點 線條清楚
%     paramSM.tempoRelMin = 0.5;
%     paramSM.tempoRelMax = 2;
%     paramSM.tempoNum = 7;
%     paramSM.forwardBackward = 1;
%     paramSM.circShift = [0:11];
% 
%     [S,I] = features_to_SM(f_CENS,f_CENS,paramSM); 
%     % [S,I] = features_to_SM(overlapCurve,overlapCurve,paramSM); 
% 
%     paramVis.colormapPreset = 2;
%     visualizeSM(S,paramVis);
%     title('S');
% 
%     visualizeTransIndex(I,paramVis);
%     title('Transposition index');


end
