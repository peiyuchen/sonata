%% ANALYZING MUSIC RECORDINGS IN SONATA FORM

% case 1 : AUTOMATED METHODS FOR ANALYZING MUSIC RECORDINGS IN SONATA FORM (PAPER 2013) 
%          MIDI FILE -> CENS feature -> Self-similarity matrix -> Thumbnail
% case 2 : combine the "bar pattern"
%          MIDI FILE -> Bar pattern  -> Self-similarity matrix -> Thumbnail


clear all; close all; clc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  0. Loads a midi file & Choose method : case 1 or case 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
index = 17;%[1,2,3,4,5,8,16,17,20];
for song_no = 1:length(index)
    
close all; clc; clearvars -except index song_no
song_no
filename = ['b_' num2str(index(song_no)) '_1_CD'];
% filename = ['b_' num2str(index(song_no)) '_1'];

% midi_R = getNoteInWhichBar(filename, 1);
% midi_L = getNoteInWhichBar(filename, 2);

SAVE    =  0;

% [f_audio, sideinfo] = wav_to_audio('', '../midi2audio/', [filename '.wav']);
[f_audio, fs] = audioread(['../midi2audio/' filename '.m4a']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  1. add paths of toolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('toolbox/MATLAB_SM-Toolbox_1.0/SMtoolbox_functions');
addpath('toolbox/MATLAB_SM-Toolbox_1.0/MATLAB-Chroma-Toolbox_2.0');
addpath('toolbox/MIDI tool/miditoolbox');
addpath('toolbox/midi_lib/midi_lib');
addpath('toolbox/matlab-midi-master/src');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  2.  Computes chroma features (CENS variant with a feature resolution 
%      of 2 Hertz). The used functions are part of the Chroma Toolbox
%      http://www.mpi-inf.mpg.de/resources/MIR/chromatoolbox/ 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fs           = 0.01; % sec
% ws           =  0.2;
% hs           = ws/2;
% midi_Info    = readmidi_java( ['../midi/' filename '.mid'] );
% f_Chrama     = midi_chrama( midi_Info, fs, ws, hs );
% S_matrix     = similar_matrix( f_Chrama, f_Chrama, 'cosine', 1 );


paramPitch.winLenSTMSP      = 4410;
[f_pitch]   = audio_to_pitch_via_FB(f_audio,paramPitch);

paramCENS.winLenSmooth      = 11;
paramCENS.downsampSmooth    = 5;
[f_CENS]    = pitch_to_CENS(f_pitch,paramCENS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   3. Computes and visualizes an enhanced and thresholded similarity 
%      matrix. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
paramSM.smoothLenSM         =     20;
paramSM.tempoRelMin         =    0.5;
paramSM.tempoRelMax         =      2;
paramSM.tempoNum            =      7;
paramSM.forwardBackward     =      1;
paramSM.circShift           = [0:11];
[S,I] = features_to_SM(f_CENS,f_CENS,paramSM);
paramVis.colormapPreset     =      2;
visualizeSM(S,paramVis);
title('S');

visualizeTransIndex(I,paramVis);
title('Transposition index');

paramThres.threshTechnique  =      2;
paramThres.threshValue      =   0.15;
paramThres.applyBinarize    =      0;
paramThres.applyScale       =      1;
paramThres.penalty          =     -2;
[S_final] = threshSM(S,paramThres);  

paramVis.imagerange         = [-2,1];
paramVis.colormapPreset     =      3;
handleFigure = visualizeSM(S_final, paramVis);
title('Final S with thresholding for computing the scapeplot matrix');
 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   4. Computes and saves a fitness scape plot.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute fitness scape plot and save
parameter.dirFitness    = 'data_fitness/';
parameter.saveFitness   = 1;
parameter.title         = filename;
parameter.fitFileName   = [filename '_fit'];
parameter.fitnessSaveFileName = [parameter.dirFitness, filename,'_fit','.mat'];
%fitnessSaveFileName = [parameter.dirFitness,parameter.fitFileName,'.mat'];

%-----------!!IMPORTANT!!--------------------------------------------------%
% For fast computing of fitness scape plot, please enable parallel computing.
% To enable that, use command 'matlabpool open'.
% To disable that, use command 'matlabpool close'
%--------------------------------------------------------------------------%
if ~exist(parameter.fitnessSaveFileName)
    k=1
    [fitness_info,parameter] = SSM_to_scapePlotFitness(S_final, parameter);
    fitness_matrix           = fitness_info.fitness;
else % % instead of computing fitness, you can load a previously computed scape plot:
    k=2
    fitnessFile = load(parameter.fitnessSaveFileName);
    fitness_matrix = fitnessFile.fitness_info.fitness;
end

paramVisScp = [];
% paramVisScp.timeLineUnit = 'sample';
% paramVisScp.timeLineUnit = 'second'; paramVisScp.featureRate = ... 
[h_fig_scapeplot,x_axis,y_axis] = visualizeScapePlot(fitness_matrix,paramVisScp);
title('Fitness scape plot','Interpreter','none');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  5. Computes the thumbnail 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute thumbnail with length constraint
 parameter.len_min_seg_frame = ceil(size(S_final,1)/6);   % peiyu
% parameter.len_min_seg_frame= 20;
[thumb_frame] = scapePlotFitness_to_thumbnail(fitness_matrix,parameter);

% show corresponding thumbnail point in fitness scape plot
center_thumb_frame = floor((thumb_frame(1) + thumb_frame(2))/2);
length_thumb_frame = thumb_frame(2) - thumb_frame(1) + 1;

point_x_pos = x_axis(center_thumb_frame);
point_y_pos = y_axis(length_thumb_frame);

hold on;
plot(point_x_pos,point_y_pos,'o','LineWidth',2,'color',[1, 0.5, 0]);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   6. Computes optimal path family and induced segment family
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% find repetitions of thumbnail
[induced_frame,pathFamily] = thumbnailSSM_to_pathFamily(thumb_frame,S_final,parameter);

paramVisPathSSM                          =  [];
paramVisPathSSM.visualizeInducedSegments =  1;
paramVisPathSSM.visualizeWarpingpath     =  1;
paramVisPathSSM.dirFigure                = 'data_image/';
paramVisPathSSM.print                    = 1;
paramVisPathSSM.figureName               = [filename '_pathFamily'];
saveas(h_fig, ['//data_image/' filename '_result.jpg']);

visualizePathFamilySSM(S_final,pathFamily,paramVisPathSSM);
title('S, path family, and induced segment family');


% convert from frames to seconds
parameter.featureRate   = 10/paramCENS.downsampSmooth;
parameter.duration      = size(S_final,1)/parameter.featureRate;
induced_second          = convertSegment_frames_to_seconds(induced_frame,parameter.featureRate);
thumb_second            = convertSegment_frames_to_seconds(thumb_frame,parameter.featureRate);


% 這在幹嘛？
% attach audio file to SSMPathFamily
if isfield(parameter,'timeLineUnit') && (strcmp(parameter.timeLineUnit,'second'))
    parameterMPP.featureTimeResType = 'seconds';
else
    parameterMPP.featureTimeResType = 'features';
end
parameterMPP.featureRate = parameter.featureRate;
parameterMPP.fs = sideinfo.wav.fs;
h_fig = gcf;
makePlotPlayable(f_audio, h_fig, parameterMPP);



% assign label to each repetition and wrap up in segment struct
% computedSegments = wrapUpSegmentInStruct(induced_frame,thumb_frame);   % peiyu
computedSegments = wrapUpSegmentInStruct(induced_second,thumb_second);
paramVisSegFam = [];
% paramVisSegFam.duration = size(S_final,1);  % peiyu
paramVisSegFam.duration = parameter.duration;  % 畫圖x最大範圍

% 畫計算過後的結果
paramVisSegFam.showLabelText    = 1;
paramVisSegFam.segType          = 'computed';
visualizeSegFamily(computedSegments,paramVisSegFam);
title('Computed segmentation');


% attach audio file to segment family visualization
parameterMPP.featureRate = parameter.featureRate;
parameterMPP.fs = sideinfo.wav.fs;
parameterMPP.featureTimeResType = 'seconds';
h_fig = gcf;
makePlotPlayable(f_audio, h_fig, parameterMPP);
% by left clicking on the x-axis of the figure, the playback will
% jump to the clicked position.
% by right clicking on the x-axis of the figure, the player will stop.

%{
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   7. Loads ground truth segmentation and compares with computed 
%      segmentation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% reading ground truth from txt file:
dirAnnotation       = '../data_annotation/';
parameter.title     = filename;
groundTruth_struct  = paper_annotation([dirAnnotation 'pei_anno/' filename '_s.txt'],1);
% groundTruth_struct  = parseAnnotationFile([dirAnnotation parameter.title '_' unit{CASE} '.txt']);

% paramVisSegFam.segType = 'groundtruth';
% visualizeSegFamily(groundTruth_struct,paramVisSegFam);
% title('Ground truth segmentation');
% h_fig = gcf;
% makePlotPlayable(f_audio, h_fig, parameterMPP);

% show ground truth and computed result together
figure;
h_fig = subplot(2,1,1);
paramVisSegFam.segType = 'groundtruth';
visualizeSegFamily(groundTruth_struct,paramVisSegFam,h_fig);
title('Ground truth segmentation');
h_fig = subplot(2,1,2);
paramVisSegFam.segType = 'computed';
visualizeSegFamily(computedSegments,paramVisSegFam,h_fig);
title('Computed segmentation');
% makePlotPlayable(f_audio, h_fig, parameterMPP);

% saveas(h_fig, ['data_image/' filename '_result.jpg']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   8. Compute ans save Precision、Recall、F-score
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:size(groundTruth_struct, 2)
    GT_frame(i,1) = groundTruth_struct(i).start;
    GT_frame(i,2) = groundTruth_struct(i).end;
end
GT_frame = convertSegment_seconds_to_frames(GT_frame, parameter.featureRate);
% for i=1:size(computedSegments, 2)
%     COM(i).start = max(1,ceil(computedSegments(i).start/0.5));
%     COM(i).end   = max(1,ceil(computedSegments(i).end/0.5));
% end

Ground  = zeros(max([max(GT_frame), max(induced_frame),1]), 1);
Compute = zeros(max([max(GT_frame), max(induced_frame),1]), 1);

for i=1:size(groundTruth_struct, 2)
    Ground (GT_frame(i,1):GT_frame(i,2), 1) = 1;   
end
for i=1:size(computedSegments, 2)
    Compute(induced_frame(i,1) : induced_frame(i,2), 1) = 1;
end

Evaluation.true_positive  = sum( Ground .* Compute  );
Evaluation.false_negtive  = sum((Ground-Compute)== 1);
Evaluation.false_positive = sum((Ground-Compute)==-1);

Evaluation.Precision = round(Evaluation.true_positive/(Evaluation.true_positive+Evaluation.false_positive) , 4);
Evaluation.Recall    = round(Evaluation.true_positive/(Evaluation.true_positive+Evaluation.false_negtive)  , 4);
Evaluation.F_score   = round(2*Evaluation.Precision *Evaluation.Recall/(Evaluation.Precision+Evaluation.Recall), 4);

% eval(['save ../evaluation/' filename '_Coarse_result_PAPER computedSegments groundTruth_struct Evaluation']);

%{
    %% 精細分析
    % 紀錄pathFamily每個對應的bar的相似度從哪個similar matrix來的
    
    [GT_subject] = paper_annotation(['../data_annotation/pei_anno/' filename '_s.txt'], 2);
    for i=1:size(groundTruth_struct, 2)
        GT_fine_M1(i,:) = [GT_subject(i).M1_start, GT_subject(i).M1_end];
        GT_fine_M2(i,:) = [GT_subject(i).M2_start, GT_subject(i).M2_end];
    end
%     GT_fine_M1_frame = convertSegment_seconds_to_frames(GT_fine_M1, parameter.featureRate);
%     GT_fine_M2_frame = convertSegment_seconds_to_frames(GT_fine_M2, parameter.featureRate);
     
    f1 = figure;
    for p = 1:1%length(pathFamily)
%         subplot(2,2,p);
        pathFamily_sec{p,1} = (pathFamily{p,1}-1)/parameter.featureRate;
        
        for j=1:length(pathFamily{p,1})
            pathFamily{p,1}(3,j) = mod(12-I(pathFamily{p,1}(1,j),pathFamily{p,1}(2,j)),12); 
            
            plot(pathFamily_sec{p,1}(1,j), pathFamily{p,1}(3,j), '.', 'MarkerSize', 10, 'Color', 'b'); hold all;          
        end
        % 畫第一主題 開始、結束 ， 第二主題 開始、結束
        for gt = 1:length(groundTruth_struct)
            line([GT_fine_M1(gt,1) GT_fine_M1(gt,1)],[0 12],'Color','g','LineWidth',1); hold all;
            line([GT_fine_M1(gt,2) GT_fine_M1(gt,2)],[0 12],'Color','g','LineWidth',2); hold all;
            line([GT_fine_M2(gt,1) GT_fine_M2(gt,1)],[0 12],'Color','r','LineWidth',2); hold all;
            line([GT_fine_M2(gt,2) GT_fine_M2(gt,2)],[0 12],'Color','r','LineWidth',1); hold all;
        end
        xlabel('Time (sec)');   ylabel('Shift index');  
        axis([pathFamily_sec{p,1}(1,1), pathFamily_sec{p,1}(1,end), 0, 11]);  
        title(['no.' filename(3:4)]);
        
               
    end
    
%     saveas(f1, ['data_image/' filename '_tranIndex.jpg']); 
%}
%}
end