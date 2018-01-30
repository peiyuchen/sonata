%% ANALYZING MUSIC RECORDINGS IN SONATA FORM

% case 1 : AUTOMATED METHODS FOR ANALYZING MUSIC RECORDINGS IN SONATA FORM (PAPER 2013) 
%          MIDI FILE -> CENS feature -> Self-similarity matrix -> Thumbnail
% case 2 : combine the "bar pattern"
%          MIDI FILE -> Bar pattern  -> Self-similarity matrix -> Thumbnail


clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P = pattern
% S = similar matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  0. Loads a midi file & Choose method : case 1 or case 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
index = [1,2,3,4,8,16,17,20];
for song_no = 1%length(index)
clearvars -except index song_no
filename = ['b_' num2str(index(song_no)) '_1'];
% filename = 'b_1_1'; % b_17_1, mz_545_1

CASE    =  2;
SAVE    =  0;
NONO    =  1;

for NONO=3
S_no    =  NONO;
S_name  = {'max3', 'max2', 'right', 'left', 'com'};

% multiProcess   = 'lower pitch'; % old
multiProcess   = 'all path';    % new

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  1. add paths of toolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('toolbox/MATLAB_SM-Toolbox_1.0/SMtoolbox_functions');
addpath('toolbox/MATLAB_SM-Toolbox_1.0/MATLAB-Chroma-Toolbox_2.0');
addpath('toolbox/MIDI tool/miditoolbox');
addpath('toolbox/midi_lib/midi_lib');
addpath('toolbox/matlab-midi-master/src');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  2. CASE 1 : Computes chroma features (CENS variant with a feature resolution 
%              of 2 Hertz). The used functions are part of the Chroma Toolbox
%              http://www.mpi-inf.mpg.de/resources/MIR/chromatoolbox/ 
%
%     CASE 2 : Computes Bar pattern
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if CASE == 1   % case 1 : CENS feature


    fs           = 0.01; % sec
    ws           =  0.2;
    hs           = ws/2;
    midi_Info    = readmidi_java( ['../midi/' filename '.mid'] );
    f_Chrama     = midi_chrama( midi_Info, fs, ws, hs );
    % S_matrix     = similar_matrix( f_Chrama, f_Chrama, 'cosine', 1 );

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
    [S,I] = features_to_SM(f_Chrama,f_Chrama,paramSM);
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
    handleFigure = visualizeSM(S_final,paramVis);
    title('Final S with thresholding for computing the scapeplot matrix');
 
elseif CASE == 2  % case 2 : bar curve 

    [midi_right, time_sig]  = midi_Preprocess(filename, 1); % 讀右手note資訊
     midi_left              = midi_Preprocess(filename, 2); % 讀左手note資訊

    % 左右手的自相似矩陣
    
    if strcmp(multiProcess, 'lower pitch') % old : 取最低音
        P_right         = pitch_Pattern_v3( midi_right, time_sig);
        P_left          = pitch_Pattern_v3( midi_left , time_sig);
    
        S.right         = similar_matrix(P_left , P_right,'correlation',1);
        S.left          = similar_matrix(P_left , P_left ,'correlation',1);
        S.com           = similar_matrix(P_right, P_left ,'correlation',1);
        
    end
    
    if strcmp(multiProcess, 'all path')    % new : 所有路徑   
        P_right  = multiPitch_Pattern_v2( midi_right, time_sig );
        P_left   = multiPitch_Pattern_v2(  midi_left, time_sig );

        S.right = barPath_S_matrix(P_right, P_right, 'correlation');
        S.left  = barPath_S_matrix(P_left,  P_left , 'correlation');
        S.com   = barPath_S_matrix(P_right, P_left , 'correlation');
    end
    
    % 休止符
    rest_right      = rest_detection( midi_right, time_sig );
    rest_left       = rest_detection(  midi_left, time_sig );

    S.rest_right    = similar_matrix( rest_right, rest_right, 'hamming' );
    S.rest_left     = similar_matrix(  rest_left,  rest_left, 'hamming' );
    S.rest_com      = similar_matrix( rest_right,  rest_left, 'hamming' );
        
    matrix_len = max([length(S.right), length(S.left), length(S.rest_right), length(S.rest_left)]) + 1;
    
    S.right( matrix_len, matrix_len )      = 0;
    S.left ( matrix_len, matrix_len )      = 0;
    S.com  ( matrix_len, matrix_len )      = 0;      
    S.rest_right( matrix_len, matrix_len ) = 0;
    S.rest_left ( matrix_len, matrix_len ) = 0;
    S.rest_com  ( matrix_len, matrix_len ) = 0;
    
    % 左右手不變自相似矩陣    
    S.right = S.right .* S.rest_right;
    S.left  = S.left  .* S.rest_left ;
    S.com   = S.com   .* S.rest_com  ;

    S.max2 = max( S.right, S.left ); 
    S.max3 = max( S.max2 , S.com  );
    
    % 找 Smax index 
    if S_no == 1
        tmp_r = S.max3 - S.right; tmp_r(tmp_r==0)=-1; tmp_r(tmp_r~=-1)=0;
        tmp_l = S.max3 - S.left;  tmp_l(tmp_l==0)=2;  tmp_l(tmp_l~=2)=0;
        tmp_c = S.max3 - S.com;   tmp_c(tmp_c==0)=3;  tmp_c(tmp_c~=3)=0;
        S.index = zeros(size(S.max3)); S.index(find(tmp_c==3))=3; S.index(find(tmp_l==2))=2; S.index(find(tmp_r==-1))=1;
    elseif S_no == 2
        tmp_r = S.max3 - S.right; tmp_r(tmp_r==0)=-1; tmp_r(tmp_r~=-1)=0;
        tmp_l = S.max3 - S.left;  tmp_l(tmp_l==0)=2;  tmp_l(tmp_l~=2)=0;
        S.index = zeros(size(S.max2)); S.index(find(tmp_l==2))=2; S.index(find(tmp_r==-1))=1;
    elseif S_no == 3
        S.index = ones(size(S.right)); 
    elseif S_no == 4
        S.index = 2 * ones(size(S.left)); 
    elseif S_no == 5
        S.index = 3 * ones(size(S.com));
    end
    
    S_th = eval( [ 'S.' S_name{S_no}; ] );   
    figure, imagesc(S_th);    title( S_name{S_no} );   % colormap(flipud(gray)); 

  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   3. Computes and visualizes an thresholded similarity matrix. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    th                  = 0.85; % 取門檻值
    S_th( S_th < th )   = -2  ;
    S_th( S_th >= th )     = (S_th(S_th>=th)-min(min(S_th(S_th>=th)))) / ( max(max(S_th(S_th>=th))) - min(min(S_th(S_th>=th))) );
    figure, imagesc(S_th);  title(['S\_th = ' num2str(th)]);  colormap(flipud(gray)); 

    S_final = S_th;

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   4. Computes and saves a fitness scape plot.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% compute fitness scape plot and save
% parameter.dirFitness    = 'data_fitness/';
% parameter.saveFitness   = 1;
% parameter.title         = filename;

%-----------!!IMPORTANT!!--------------------------------------------------%
% For fast computing of fitness scape plot, please enable parallel computing.
% To enable that, use command 'matlabpool open'.
% To disable that, use command 'matlabpool close'
%--------------------------------------------------------------------------%
[fitness_info,parameter] = SSM_to_scapePlotFitness(S_final);
fitness_matrix           = fitness_info.fitness;

% % instead of computing fitness, you can load a previously computed scape plot:
% fitnessSaveFileName = ['data_fitness/',filename(1:end-4),'_fit','.mat'];
% fitnessFile = load(fitnessSaveFileName);
% fitness_matrix = fitnessFile.fitness_info.fitness;

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

paramVisPathSSM                          = [];
paramVisPathSSM.visualizeInducedSegments =  1;
paramVisPathSSM.visualizeWarpingpath     =  1;
visualizePathFamilySSM(S_final,pathFamily,paramVisPathSSM);
title('S, path family, and induced segment family');


% convert from frames to seconds
% parameter.featureRate = 10/paramCENS.downsampSmooth;
% parameter.duration = size(S_final,1)/parameter.featureRate;
% induced_second = convertSegment_frames_to_seconds(induced_frame,parameter.featureRate);
% thumb_second = convertSegment_frames_to_seconds(thumb_frame,parameter.featureRate);



% attach audio file to SSMPathFamily
% if isfield(parameter,'timeLineUnit') && (strcmp(parameter.timeLineUnit,'second'))
%     parameterMPP.featureTimeResType = 'seconds';
% else
%     parameterMPP.featureTimeResType = 'features';
% end
% parameterMPP.featureRate = parameter.featureRate;
% parameterMPP.fs = sideinfo.wav.fs;
% h_fig = gcf;
% makePlotPlayable(f_audio, h_fig, parameterMPP);



% assign label to each repetition and wrap up in segment struct
computedSegments = wrapUpSegmentInStruct(induced_frame,thumb_frame);   % peiyu
% computedSegments = wrapUpSegmentInStruct(induced_second,thumb_second);
paramVisSegFam = [];
paramVisSegFam.duration = size(S_final,1);
% paramVisSegFam.duration = parameter.duration;  % 畫圖x最大範圍

% 畫計算過後的結果
% paramVisSegFam.showLabelText = 1;
% paramVisSegFam.segType = 'computed';
% visualizeSegFamily(computedSegments,paramVisSegFam);
% title('Computed segmentation');


% attach audio file to segment family visualization
% parameterMPP.featureRate = parameter.featureRate;
% parameterMPP.fs = sideinfo.wav.fs;
% parameterMPP.featureTimeResType = 'seconds';
% h_fig = gcf;
% makePlotPlayable(f_audio, h_fig, parameterMPP);
% by left clicking on the x-axis of the figure, the playback will
% jump to the clicked position.
% by right clicking on the x-axis of the figure, the player will stop.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   7. Loads ground truth segmentation and compares with computed 
%      segmentation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% reading ground truth from txt file:
unit = {'sec','bar'};
dirAnnotation       = '../data_annotation/';
parameter.title     = filename;
groundTruth_struct  = parseAnnotationFile([dirAnnotation parameter.title '_' unit{CASE} '.txt']);

% paramVisSegFam.segType = 'groundtruth';
% visualizeSegFamily(groundTruth_struct,paramVisSegFam);
% title('Ground truth segmentation');
% h_fig = gcf;
% makePlotPlayable(f_audio, h_fig, parameterMPP);

for i=1:size(computedSegments, 2)
    if CASE == 1
        computedSegments(i).start = (computedSegments(i).start-1)*(hs*5);
        computedSegments(i).end   = computedSegments(i).end*(hs*5);
    end
    Compute(computedSegments(i).start:computedSegments(i).end,1) = 1;
end

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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   8. Compute ans save Precision、Recall、F-score
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% frame 轉 sec
if CASE == 1
    [(induced_frame(:,1)-1) (induced_frame(:,2))].*(hs*5);
end
Ground  = zeros(size(S_final,1), 1);
Compute = zeros(size(S_final,1), 1);

for i=1:size(groundTruth_struct, 2)
    Ground (groundTruth_struct(i).start:groundTruth_struct(i).end,1) = 1;   
end
for i=1:size(computedSegments, 2)
%     if CASE == 1
%         computedSegments(i).start = (computedSegments(i).start-1)*(hs*5);
%         computedSegments(i).end   = computedSegments(i).end*(hs*5);
%     end
    Compute(computedSegments(i).start:computedSegments(i).end,1) = 1;
end

Evaluation.true_positive  = sum( Ground .* Compute  );
Evaluation.false_negtive  = sum((Ground-Compute)== 1);
Evaluation.false_positive = sum((Ground-Compute)==-1);

Evaluation.Precision = round(Evaluation.true_positive/(Evaluation.true_positive+Evaluation.false_positive) , 4);
Evaluation.Recall    = round(Evaluation.true_positive/(Evaluation.true_positive+Evaluation.false_negtive)  , 4);
Evaluation.F_score   = round(2*Evaluation.Precision *Evaluation.Recall/(Evaluation.Precision+Evaluation.Recall), 4);

% save
if SAVE
    if CASE == 1
        fileID = fopen(['../evaluation/' filename '_CENS.txt'],'w');  
    elseif CASE == 2
        fileID = fopen(['../evaluation/' filename '_' S_name{S_no} '.txt'],'w');
    end

    word_method = 'Bar pattern';    if CASE == 1; word_method = 'CENS feature'; end
    fprintf(fileID,['Method        : ' word_method  '\n\n']);
    fprintf(fileID,['SSM           : ' S_name{S_no} '\n\n']);
    fprintf(fileID,['multi-process : ' multiProcess '\n\n']);

    fprintf(fileID,['\n']);
    fprintf(fileID,['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%']);
    fprintf(fileID,['\n\n']);

    for i=1:size(groundTruth_struct, 2)
        fprintf(fileID,['GT     :  %3d %3d \n'], groundTruth_struct(i).start, groundTruth_struct(i).end);
    end

    fprintf(fileID,['\n']);
    fprintf(fileID,['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%']);
    fprintf(fileID,['\n\n']);

    for i=1:size(computedSegments, 2)
        fprintf(fileID,['COMPUTE:  %3d %3d \n'], computedSegments(i).start, computedSegments(i).end);
    end

    fprintf(fileID,['\n']);
    fprintf(fileID,['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%']);
    fprintf(fileID,['\n\n']);

    fprintf(fileID,['\nTP : %d, FN : %d, FP : %d \n\n'],Evaluation.true_positive, Evaluation.false_negtive, Evaluation.false_positive);
    fprintf(fileID,['\nPrecision : %.4f, Recall : %.4f, F-score : %.4f \n\n'],Evaluation.Precision, Evaluation.Recall, Evaluation.F_score);

    fclose(fileID);
end

end

%% 精細分析
% 紀錄pathFamily每個對應的bar的相似度從哪個similar matrix來的
clear PATHFAMILY;
for p = 1:1%length(pathFamily)
    for j=1:length(pathFamily{p,1})
        pathFamily{p,1}(3,j) = S.index(pathFamily{p,1}(1,j),pathFamily{p,1}(2,j));
        
        % 比較路徑 p1 p2從哪來 right left or combine
        if pathFamily{p,1}(3,j)==1
            PATHFAMILY(p).p1(j,1) = P_right.down_PATH_Y( pathFamily{p,1}(1,j), 1 );
            PATHFAMILY(p).p2(j,1) = P_right.down_PATH_Y( pathFamily{p,1}(2,j), 1 );
        elseif pathFamily{p,1}(3,j)==2
            PATHFAMILY(p).p1(j,1) = P_left.down_PATH_Y( pathFamily{p,1}(1,j), 1 );
            PATHFAMILY(p).p2(j,1) = P_left.down_PATH_Y( pathFamily{p,1}(2,j), 1 );
        elseif pathFamily{p,1}(3,j)==3
            PATHFAMILY(p).p1(j,1) = P_right.down_PATH_Y( pathFamily{p,1}(1,j), 1 );
            PATHFAMILY(p).p2(j,1) = P_left.down_PATH_Y ( pathFamily{p,1}(2,j), 1 );
        end
    end
    if p==1
    [ PATHFAMILY_v2 ]      = ER_key_relation_v2( PATHFAMILY(p).p1, PATHFAMILY(p).p2, pathFamily{p,1}, filename );
    end
end


% [ PATHFAMILY1 ]      = ER_key_relation( P_right.down_PATH_Y, P_left.down_PATH_Y, pathFamily, filename );

% close all;
end
