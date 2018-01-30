%%%% motify : "Algorithms for Choral Analysis 2002"
%%%% date   : 18/01/22
%%%% content: �M�����R�t��k�A�H�p�`�����i����R�C
%%%% input  : midi data
%%%% output : chords in each bar

clear; close all; clc;

% include toolbox
addpath('toolbox/midi_lib/midi_lib');

% input midi data
filePath = '../midi';
fileName = 'mz_545_1_noRepeat';
[midiData, timeSig]  = midi_Preprocess(fileName);

midiData = normalize_midi_data(midiData);

% each bar's information
[barNote, barStart] = bar_note_data(midiData, timeSig);

% �j�p��
keyChord = make_key_chord;

% parameter
tempName    = {'maj','7','min','dim','xxx','X'};%,'Fully dim7','Half dim7','dim3','X'}; chordType
pitchName   = {'C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B'};
keyName     = {'C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B', ...
               'c', 'c#', 'd', 'd#', 'e', 'f', 'f#', 'g', 'g#', 'a', 'a#', 'b'};
%                 �j�T        �ݤC         �p�T       ��
template    = [ 0,4,7,-1; 0,4,7,10; 0,3,7,-1; 0,3,6,-1];  
evaI = 1;
evaChord = {'�p�`','���(onset)','�թ�','�M���W��','�M���s��','�Ƶ�'};% new




%%%%%%%%%%% �[�ק諸��k %%%%%%%%%%
isDownbeatAndDurWeight      = 0;
isLowWeight                 = 0;
isDownbeatWeight            = 0;
isSave                      = 1;
ratio                       = 1;


for j = 1:length(barNote)

    noteBar = barNote{j};
    
    % ���̤p���q�����I�����ǭ�
    [pitchClass, miniSeg] = min_seg_pitchclass(noteBar, barStart(j));   
    if isDownbeatAndDurWeight
        [pitchClass, lowestPClass, lowestP, miniSeg] = downbeat_pitchclass(noteBar, barStart(j), isDownbeatWeight, ratio);
    end
    %% �����Ҧ��i����q������ 
%     [maxScoreMatrix, scoreMatrixIdx, scoreMatrix] = record_min_seg_score(pitchClass);    

    templateNum  = size(template, 1);
    partitionNum = size(pitchClass, 1) + 1;
    rootS  = -inf * ones(partitionNum, partitionNum, templateNum * 12);
    P      = -inf * ones(partitionNum, partitionNum, templateNum * 12);
    N      = -inf * ones(partitionNum, partitionNum, templateNum * 12);
    M      = -inf * ones(partitionNum, partitionNum, templateNum * 12);
    S      = -inf * ones(partitionNum, partitionNum, templateNum * 12);
    
    lowestPitchP = zeros(partitionNum, partitionNum, templateNum * 12);
    lowestPitchN = -inf * ones(partitionNum, partitionNum, templateNum * 12);

    
     for ii = 1:partitionNum
         for jj = ii+1:partitionNum             
             if isLowWeight && min(lowestP(ii:jj-1, 1)) ~= Inf
                lowestPitchIdx = find(lowestP == min(lowestP(ii:jj-1, 1)));
             end
             for kk = 1:templateNum
                 for p = 1:12
                     templateNo = mod(template(kk,(template(kk,:)~=-1)) + p - 1 , 12);
                     Template = zeros(1,12); 
                     Template(1,templateNo + 1) = 1;
                     
                     rootS(ii, jj, (kk - 1) * 12 + p) = sum(sum(pitchClass(ii:jj-1, templateNo(1) + 1)));                 % �ڭ���score
                     P(ii, jj, (kk - 1) * 12 + p)     = sum(sum(pitchClass(ii:jj-1, templateNo + 1)));                  % ���q ��, �M�� ��
                     N(ii, jj, (kk - 1) * 12 + p)     = sum(sum(pitchClass(ii:jj-1, setdiff(1:12,templateNo+1))));    % ���q ��, �M���S��  
%                      N(ii,jj,(kk-1)*12+p)     = sum(sum(pitchClass(ii:jj-1,:)))-P(ii,jj,(kk-1)*12+p);    % ���q ��, �M���S��
                     M(ii, jj, (kk - 1) * 12 + p)     = sum((sum(pitchClass(ii:jj-1, :),1) == 0) .* Template);        % ���q�S��, �M�� ��
                     S(ii, jj, (kk - 1) * 12 + p)     = ...
                         P(ii, jj, (kk - 1) * 12 + p) - (N(ii, jj, (kk - 1) * 12 + p) + M(ii, jj, (kk-1) * 12 + p));
                     %%%%% �[�̧C��weight %%%%%%
                     if  isLowWeight && min(lowestP(ii:jj-1, 1)) ~= Inf
                         lowestPitchP(ii, jj, (kk - 1) * 12 + p)  = sum(sum(lowestPClass(lowestPitchIdx, templateNo + 1)));
                         lowestPitchN(ii, jj, (kk - 1) * 12 + p)  = sum(sum(lowestPClass(lowestPitchIdx, setdiff(1:12, templateNo + 1))));
                         S(ii, jj, (kk - 1) * 12 + p) = S(ii, jj, (kk - 1) * 12 + p) + ...
                             lowestPitchP(ii, jj, (kk - 1) * 12 + p) - lowestPitchN(ii, jj, (kk - 1) * 12 + p);
                     end
                     %%%%%%%%%%%%%%%%%%%%%%%%%%
                 end
             end
         end
     end

    [maxScoreMatrix, scoreMatrixIdx ] = max(S, [], 3);

    %% ���q�G�̨θ��|�t��k HarmAn
    % &&&&&&&&&&& .   downbeat  &&&&&&&&&&&&&&
    MARK = HarmAn(maxScoreMatrix, partitionNum);

%     partitionNum      = size(pitchClass,1)+1;
%     downbeat          = find(miniSeg(:,1)==fix(miniSeg(:,1)));
%     MARK = HarmAn(maxScoreMatrix, partitionNum, [downbeat; partitionNum]);
    % &&&&&&&&&&& .   downbeat  &&&&&&&&&&&&&&

    % ���� �̨θ��|��onset
    segOnset = [miniSeg(:,1)'-barStart(j), miniSeg(end,2)-barStart(j)];

    Ans.Mark    (j,1:length(MARK)) = MARK;
    Ans.segOnset(j,1:length(MARK)) = segOnset(MARK);

    %% �̨θ��|���M���аO
    Wtmp = 0;
    Stmp = 0;
    for k = 2:length(MARK)
        Ans.Chord{j,k-1}       = scoreMatrixIdx(MARK(k-1),MARK(k));
        Ans.sameScore(j,k-1)   = {''};
        Ans.Score(j,k-1)       = maxScoreMatrix(MARK(k-1),MARK(k));

        Wtmp                   = Wtmp + maxScoreMatrix(MARK(k-1),MARK(k));
        Stmp                   = Stmp + maxScoreMatrix(MARK(k-1),MARK(k)) * diff(miniSeg(k-1,:))/(miniSeg(end,end)-miniSeg(1,1));

        if k == length(MARK)
            Ans.ScoreSum(j,1) = Wtmp;                   % �p�`�M���`����
            Ans.ScorePar(j,1) = Wtmp/partitionNum;      % �p�`�M����������
            Ans.ScoreLen(j,1) = Stmp;                   % 
            Ans.ScoreLenPar(j,1) = Ans.ScoreLen(j,1)/partitionNum;
        end


        %% ���Ƥ@�˪���
        HighScoreChord = find(S(MARK(k-1),MARK(k),:)==maxScoreMatrix(MARK(k-1),MARK(k)));
        % step 1 : highest root scoreMatrix
        if numel(unique(HighScoreChord))~=1
            HighScoreRootScore  = rootS(MARK(k-1),MARK(k),HighScoreChord);
            MaxRootIdx          = find(HighScoreRootScore==max(HighScoreRootScore));
            HighScoreChord      = HighScoreChord(MaxRootIdx);
            step                = {'root'};
            % step 2 : highest probability of occurrence
            if numel(unique(HighScoreChord))~=1
                chordSpecies    = floor(HighScoreChord/12)+1;
                MaxProbIdx      = find(chordSpecies==min(chordSpecies));
                HighScoreChord  = HighScoreChord(MaxProbIdx); 
                step            = {'probability'};
                if numel(unique(HighScoreChord))~=1        
                    HighScoreChord = 72
                    step     = {'no solve'};
                    bar      = j
                end
            end

            Ans.Chord{j,k-1}     = HighScoreChord;
            Ans.sameScore(j,k-1) = step;
        end
        
        Ans.templ_no(j,k-1)  = ceil(Ans.Chord{j,k-1}(1)/12);
        Ans.pitch_no{j,k-1}  = ~(ceil(mod(Ans.Chord{j,k-1},12)/12))*12 + mod(Ans.Chord{j,k-1},12);
        Ans.ChordName{j,k-1} = strcat(pitchName(Ans.pitch_no{j,k-1}),':',tempName(Ans.templ_no(j,k-1)));
        
        evaI = evaI + 1;
        evaChord{evaI,1} = j;
        evaChord{evaI,2} = segOnset(MARK(k-1));
        evaChord(evaI,4) = Ans.ChordName{j,k-1};
        evaChord{evaI,5} = HighScoreChord;
        

    end

            %%
            %{
            
    % �⦹�p�`���M���b24�դ� �ݩ�ӽժ��M���ƶq
    for key_i=1:24
        ansBarChord = Ans.Chord(j,:);
        DIFF = numel(setdiff(ansBarChord(ansBarChord~=0),keyChord(key_i,:)));
        num = numel(ansBarChord(ansBarChord~=0));
        SCORE(j,key_i) = (num-DIFF)/num;
%         if key_i==2 && j ==8
%             setdiff(ansBarChord(ansBarChord~=0),keyChord(key_i,:))
%             ansBarChord
%             keyChord(key_i,:)
%             DIFF
%             num
% 
%         end
    end
    Cost          = [cellstr(keyName)' num2cell(SCORE(j,:))'];        
    [~,s]         = sort(SCORE(j,:)*-1);                          % which cost to using to rank
    Rank(:,:,j)   = [cellstr(keyName(s))' num2cell(SCORE(j,s))' num2cell(s)']; % rank, key name & cost in each bar
    costRank(j,:) = cell2mat(reshape(Rank(:,2,j),[],1));          % ranked cost in each bar
    nameRank(j,:) = cellstr (reshape(Rank(:,1,j),[],1))';         % ranked name in each bar
    nameIdxRank(j,:) = cell2mat(reshape(Rank(:,3,j),[],1));  
    costKey(j,:)  = SCORE(j,:);                                   % key cost    in each bar

    TMPnameRank = nameRank(j,:);
    TMPnameRank(costRank(j,:)==0) = [];
    nameRank_v2(j,1:length(TMPnameRank)) = TMPnameRank;
                
    score_rank = sort(unique(costRank(j,:)),'descend');
    if score_rank(1)~=0
        TMPnameRank = nameRank(j,:);
        tmp = TMPnameRank(costRank(j,:)==score_rank(1));
        nameRank1(j,1:length(tmp)) = tmp;
        keyRank(j,1:24) = 0;
        keyRank(j,SCORE(j,:)==score_rank(1)) = 1;

        keyRatio(j,:) = KeyConsistRatio(noteBar);               % new
        keyRatioRank1(j,:) = keyRatio(j,:).*keyRank(j,:);       % new
        keySortRatio(j,:) = keyRatioRank1(j,nameIdxRank(j,:));  % new
        % �B�z ��max�d�U�� ��L���n
        MAX = max(keySortRatio(j,:));
        No1Idx = find(keySortRatio(j,:)==max(keySortRatio(j,:)));
        FinalNo1(j,1:length(No1Idx)) = nameRank1(j,No1Idx);
    end
                
    if length(score_rank)>=2
        if score_rank(2)~=0
%                         tmp = TMPnameRank(costRank(j,:)==score_rank(2));
%                         nameRank2(j,1:length(tmp)) = tmp;

%                         keyRank(j,SCORE(j,:)==score_rank(2)) = 2;
        end
    end
    %}
end

if isSave
    cell2csv(['chordEva/eva_' fileName '.csv'], evaChord)
end




%% function 


function [barNote, barSTART] = bar_note_data(midiData, timeSig)
% %               1    2    4+8    4+16     4    8+16      8      16    
% noteValue   = [ 1; 0.5; 0.375; 0.3125; 0.25; 0.1875; 0.125; 0.0625 ];
% nearDur         = dsearchn(noteValue, Note{bar_i}(:,2)/2^timeSig(1,2));
% W_dur           = noteValue(nearDur);  
    % init
    addBeat     = 0; 
    addBar      = 0; 
    midiSize    = size(midiData,2);
    
    for i = 1:size(timeSig,1)

        NoOfBeats       =   timeSig(i,1)*(4/(2^timeSig(i,2)));
        if i ~= size(timeSig,1)
            NoOfBar     =   max(ceil( ((timeSig(i+1,5)-1) - timeSig(i,5) ) / NoOfBeats), 1);
        else
            NoOfBar     =   max(ceil( ( sum(midiData(end,1:2)) - timeSig(i,5) ) / NoOfBeats), 1);
        end

        for j = 1:NoOfBar
            barStart   = (addBeat+(j-1)*NoOfBeats);
            barEnd     = addBeat+j*NoOfBeats;
            notesOfBar = intersect(find(midiData(:,1)>=barStart), find(midiData(:,1)<barEnd));

            currentBar = addBar + j;
            midiData(notesOfBar, midiSize+1) = currentBar;
            barSTART(currentBar) = barStart;


            % �e�Ӥp�`��ƽu�ܦ��p�`������
            innerOffset     = intersect(find(sum(midiData(:,1:2),2)<barEnd) , find(sum(midiData(:,1:2),2)>barStart));
            overOnset       = intersect(find(midiData(:,1)<barStart), innerOffset);
            overNote        = midiData(overOnset, :);
            overNote(:,2)   = overNote(:,2) + barStart - overNote(:,1);
            overNote(:,1)   = barStart;

            noteBar         = midiData(notesOfBar, :);      % ���p�`��������
            noteBar         = [noteBar; overNote];          % ��onset�令���p�`�@�}�l

            if ~isempty(noteBar)

                noteBar       = sortrows(noteBar,1);        % �̾ڲĤ@��Ƨ� �p��j
%                 tmp           = length(noteBar);
                noteBar       = triDetection(noteBar);      % tri �B�z
    %                 noteBar(noteBar(:,2)<0.2,:) = [];         % ��p��duration 0.25�����ų��R���A�{���Ӥp���Ť��|��M�����Ӥj���v�T�C
                barNote{currentBar,1}      = noteBar;
            end
        end

        addBar  = addBar  + NoOfBar;
        addBeat = addBeat + NoOfBeats*NoOfBar;

    end

end
%%
% description : �إߤj�p�ժ��M��
% output : �j�p�ժ������M��
function [keyChord] = make_key_chord
    load Key_DiatonicTriad_V7.mat keyFusion

    % MAJOR    I IV V V7 ii vi viio iii
    tempMaj = [1 1 1 2 3 3 4 3];         % ��3
    % MINOR    i iv V V7 iio VI viio iii+
    tempMin = [3 3 1 2 4 1 4 1];

    keyScale = [1 4 5 5 2 6 7 3];
    keyChord = zeros(24, length(tempMaj));
    for i=1:12
        keyChord(i, :) = (tempMaj - 1) * 12 + keyFusion(i).diatonicTriad(1, keyScale) + 1;
        keyChord(i + 12, :) = (tempMin - 1) * 12 + keyFusion(i + 12).diatonicTriad(1, keyScale) + 1;
    end
end

%%
% description: �H���I��partition point�A�����C��minSeg���������v���A�H��minSeg���̧C��
% input:  noteBar      -> �p�`����
%         barStart     -> �p�`�}�l��beat
% output: pitchClass,  -> �C��minSeg���Ҧ����A�Bweight�����Ū��סA�j�p��N*12 (N��minSeg�ƶq,12���@�ӤK�׭����ƶq)
%         lowestPitchC -> �C��minSeg���̧C���A�Bweight�����Ū��סA�j�p��N*12 (N��minSeg�ƶq,12���@�ӤK�׭����ƶq)
%         lowestP      -> �C��minSeg���̧C���A�B���H�@�ӤK�רӬ�(��mod12)
%         miniSeg      -> �C��minSeg��onset offset��beat
% �䥪��̧C���v�����j�A���H���륪�k��A�ҥH����q�̧C���B�P�ɦ���ӭ�

function [pitchClass, lowestPitchC, lowestP, miniSeg] = downbeat_pitchclass(noteBar, barStart, isDownbeatWeight,ratio)


% &&&&&&&&&&& .   downbeat  &&&&&&&&&&&&&&
    beatInBar   = noteBar(:,9);
    miniSeg     = [[0:beatInBar-1]' [1:beatInBar]']; % �H���I���̤p���q (�� 0 1 2 3 ��)

% &&&&&&&&&&& .   �C�Ӥ��q�����ǭ��A���Ū��׬�weight  &&&&&&&&&&&&&&
    pitchClass    = zeros(size(miniSeg,1),12);
    lowestPitchC  = zeros(size(miniSeg,1),12);
    lowestP       = inf*ones(size(miniSeg,1),1);
    for P = 1:size(miniSeg,1)
        unit        = (1/2)^4;
        On          = zeros(size(noteBar,1),1);
        Off         = zeros(size(noteBar,1),1);
        onCutFlag   = zeros(size(noteBar,1),1);
        lowestFlag  = zeros(ceil((miniSeg(P,2)-miniSeg(P,1))/unit), 1);     % ����minSeg���Cunit���סA���X�ӭ�
        lowestNO    = zeros(ceil((miniSeg(P,2)-miniSeg(P,1))/unit) ,1);     % �C��unit�̧C���O���ӭ�
        lowestPitch = inf*ones(ceil((miniSeg(P,2)-miniSeg(P,1))/unit),1);   % �C��unit�̧C��������
        for note = 1:size(noteBar,1)
           
            On(note)  = round(noteBar(note,1)-barStart,4);
            Off(note) = round(sum(noteBar(note,1:2))-barStart,4);
            onsetCut  = 0; % isDownbeatWeight
            if On(note) < miniSeg(P,1)
                onsetCut = miniSeg(P,1)-On(note);  % isDownbeatWeight
                On(note) = miniSeg(P,1); 
            end
            if Off(note)> miniSeg(P,2); Off(note)= miniSeg(P,2); end
                        
            onCutFlag(note) = (onsetCut == fix(onsetCut)); % isDownbeatWeight

            if  Off(note)-On(note)>0   % �C�ӭ����S�����b��minSeg��
                pitchNo = mod(noteBar(note,4),12)+1;
                pitchClass(P,pitchNo) = pitchClass(P,pitchNo) + Off(note)-On(note);
                
                if isDownbeatWeight && On(note) == fix(On(note)) && onCutFlag(note)
                    pitchClass(P,pitchNo) = pitchClass(P,pitchNo) + Off(note)-On(note);
                end
                % ��on�Moff�O��minSeg��������unit
                On_num  = floor((On(note) - miniSeg(P,1))/unit) + 1;
                Off_num = ceil((Off(note) - miniSeg(P,1))/unit);

                [minPitch,idx] = min([lowestPitch(On_num:Off_num) ones(Off_num-On_num+1,1)*noteBar(note,4)],[],2);
                lowestFlag(On_num:Off_num) = lowestFlag(On_num:Off_num)+1;
                if ~isempty(find(idx==2, 1))
                    lowestNO(find(idx==2)+On_num-1) = note;
                    lowestPitch(On_num:Off_num)     = minPitch;
                end
            end
            
            if note == size(noteBar,1) && ~isempty(find(lowestFlag>1, 1))
                lowestP(P)  = min(lowestPitch(lowestFlag>1));
                lowIDX      = find(lowestPitch==lowestP(P));
                lowNUM      = numel(lowIDX);
                if isDownbeatWeight
                    weight      = 0;
                    intSect     = intersect(lowestNO, find(onCutFlag==1));
                    for inS = 1:length(intSect)
                        weight = weight + numel(find(lowestNO(lowIDX)==intSect(inS)));
                    end
                    lowNUM = lowNUM + weight*ratio;
                end
                lowestPitchC(P,mod(lowestP(P),12)+1) = lowNUM * unit;
                lowNUM * unit
            end
             
        end       
    end
    miniSeg = miniSeg + barStart;
end

function [pitchClass, miniSeg] = min_seg_pitchclass(noteBar, barStart)
% &&&&&&&&&&& .   downbeat  &&&&&&&&&&&&&&
% %     beatInBar = noteBar(:,9);
% % %     miniSeg = [[0:beatInBar-1]' [1:beatInBar]'];
% %     miniSeg       = unique(round([noteBar(:,1); sum(noteBar(:,1:2),2); [0:beatInBar-1]'+barStart] ,4));
% %     miniSeg(:,2)  = [ miniSeg(2:end,1);   beatInBar(1)+barStart ];
% &&&&&&&&&&& .   downbeat  &&&&&&&&&&&&&&

    
    % �̤p���q��onset offset
    miniSeg       = unique(round([noteBar(:,1); sum(noteBar(:,1:2),2)],4));
    miniSeg(:,2)  = [ miniSeg(2:end,1);   round(sum(noteBar(end,1:2)),4) ];
    
    % ��miniSeg�Ӥp���h��
    miniTime      = 0.025;
    if miniSeg(end,1)==miniSeg(end,2); miniSeg(end,:)=[]; end      
    miniSeg(diff(miniSeg,1,2)<=miniTime,:) = [];

    % �����C�ӳ̤p���q�����ǭ�
    pitchClass    = zeros(length(miniSeg),12);
    parpointNum   = zeros(length(miniSeg), 1);

    for k=1:size(noteBar, 1)
        s               = find(abs(round(noteBar(k,1),4) - miniSeg(:,1)) <= miniTime);   % �o�ӭ��q����seg�}�l
        e               = find(round(sum(noteBar(k,1:2)),4)>=miniSeg(:,2));              % �o�ӭ��q����seg����
        parpointNum(s)  = parpointNum(s) + 1;
        pitchNum        = mod(noteBar(k,4),12);
        pitchClass(s:e(end),pitchNum+1) = pitchClass(s:e(end),pitchNum+1)+1; % miniSeg * 12
    end
end

function [maxScoreMatrix, I, S] = record_min_seg_score(pitchClass)
    %                            �j�T        �ݤC         �p�T       ��
    template            = [ 0,4,7,-1;  0,4,7,10;   0,3,7,-1; 0,3,6,-1];  
    templateNum             = size(template,1);
    partitionNum      = size(pitchClass,1)+1;
     rootS  = -inf*ones(partitionNum, partitionNum, templateNum*12);
     P      = -inf*ones(partitionNum, partitionNum, templateNum*12);
     N      = -inf*ones(partitionNum, partitionNum, templateNum*12);
     M      = -inf*ones(partitionNum, partitionNum, templateNum*12);
     S      = -inf*ones(partitionNum, partitionNum, templateNum*12);

     for ii=1:partitionNum
         for jj=ii+1:partitionNum
            for kk=1:templateNum
                for p = 1:12
                    templateNo = mod(template(kk,(template(kk,:)~=-1)) + p - 1 , 12);
                    Template = zeros(1,12); Template(1,templateNo + 1) = 1;
                    rootS(ii,jj,(kk-1)*12+p) = sum(sum(pitchClass(ii:jj-1,templateNo(1)+1))); % �ڭ���score
                    P(ii,jj,(kk-1)*12+p)     = sum(sum(pitchClass(ii:jj-1,templateNo+1  )));                  % ���q ��, �M�� ��
                    N(ii,jj,(kk-1)*12+p)     = sum(sum(pitchClass(ii:jj-1,:)))-P(ii,jj,(kk-1)*12+p);    % ���q ��, �M���S��            
                    M(ii,jj,(kk-1)*12+p)     = sum((sum(pitchClass(ii:jj-1,:),1)==0).*Template);        % ���q�S��, �M�� ��
                    S(ii,jj,(kk-1)*12+p)     = P(ii,jj,(kk-1)*12+p) - ( N(ii,jj,(kk-1)*12+p) + M(ii,jj,(kk-1)*12+p) );
                end
            end
         end
     end

    [maxScoreMatrix,I] = max(S,[],3);
end

function [keyRatio] = KeyConsistRatio(Note)
% �p�`�� �C�ӽժ��զ������
    load Key_DiatonicTriad_V7.mat

    consist     = zeros( 1,12);
    keyConsist  = zeros(24,12);
    barPitch    = mod(Note(:,4),12)+1;
    tab = tabulate(barPitch);
    consist(1:size(tab,1)) = tab(:,2);
    
    for i=1:24
        keyConsist(i, keyH(i).consistNumber+1) = 1;
        keyRatio(i) = sum(keyConsist(i,:).* consist)/size(Note,1);
    end

end

%%
% desciption: �Nmidi��duration(unit:beat)���W�Ʀ����Ф@��
% input     : midiData           -> ��l   ��midi data
% output    : normalizedMidiData -> ���W�ƫ᪺midi data
function [midiData] = normalize_midi_data (midiData)
    %                  1    2    4+8    4+16     4    8+16      8      16    
    noteValue       = [1  0.5  0.375  0.3125  0.25  0.1875  0.125  0.0625];
    repNoteValue    = repmat(noteValue, length(midiData), 1);
    tmp             = midiData(:, 2) - floor(midiData(:, 2)) - repNoteValue; 
    tmp(tmp > 0)    = -10;
    [~,idx]         = max(tmp, [], 2);
    midiData(:,2) = noteValue(idx)' + floor(midiData(:,2));
    % nearDur         = dsearchn(newNoteValue, normalizedMidiData(:,2)-floor(normalizedMidiData(:,2)));
    % W_dur           = noteValue(nearDur) + floor(normalizedMidiData(:,2)); 
end