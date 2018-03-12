%%%% motify : "Algorithms for Choral Analysis 2002"
%%%% date   : 18/01/22
%%%% content: 和弦分析演算法，以小節為單位進行分析。
%%%% input  : midi data
%%%% output : chords in each bar

clear; close all; clc;

% include toolbox
addpath('toolbox/midi_lib/midi_lib');

% input midi data
filePath = '../midi';
fileName = 'm_7_cut';%m_7_cut, mz_545_1_noRepeat
[midiData, timeSig]  = midi_Preprocess(fileName);

% each bar's information
[barNote, barOnset] = bar_note_data(midiData, timeSig);

% 大小調
keyChord = make_key_chord;

% parameter
tempName    = {'maj','7','min','dim','xxx','X'};%,'Fully dim7','Half dim7','dim3','X'}; chordType
pitchName   = {'C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B'};
keyName     = {'C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B', ...
               'c', 'c#', 'd', 'd#', 'e', 'f', 'f#', 'g', 'g#', 'a', 'a#', 'b'};
%                 大三        屬七         小三       減
template    = [ 0,4,7,-1; 0,4,7,10; 0,3,7,-1; 0,3,6,-1];  
evaluationIdx = 1;
evaluationChord = {'小節','拍數(onset)','調性','和弦名稱','和弦編號','備註'};% new




%%%%%%%%%%% 加修改的方法 %%%%%%%%%%
isDownbeatAndDurWeight      = 1;
isLowWeight                 = 0;
isDownbeatWeight            = 0;
isSave                      = 0;


for j = 1:length(barNote)
    noteBar = barNote{j};

    if ~isempty(noteBar)
        % 找到最小片段分割點有哪些音
        [pitchClass, minSeg] = min_seg_pitchclass(noteBar, barOnset(j));   
        if isDownbeatAndDurWeight
            [pitchClass, lowestPClass, lowestP, minSeg] = downbeat_pitchclass(noteBar, barOnset(j), isDownbeatWeight);
        end
        %% 紀錄所有可能分段的分數 
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
                 for kk = 1:templateNum
                     for p = 1:12
                         templateNo = mod(template(kk,(template(kk,:)~=-1)) + p - 1 , 12);
                         Template = zeros(1,12); 
                         Template(1,templateNo + 1) = 1;

                         rootS(ii, jj, (kk - 1) * 12 + p) = sum(sum(pitchClass(ii:jj-1, templateNo(1) + 1)));                 % 根音的score
                         P(ii, jj, (kk - 1) * 12 + p)     = sum(sum(pitchClass(ii:jj-1, templateNo + 1)));                  % 片段 有, 和弦 有
                         N(ii, jj, (kk - 1) * 12 + p)     = sum(sum(pitchClass(ii:jj-1, setdiff(1:12,templateNo+1))));    % 片段 有, 和弦沒有  
    %                      N(ii,jj,(kk-1)*12+p)     = sum(sum(pitchClass(ii:jj-1,:)))-s(ii,jj,(kk-1)*12+p);    % 片段 有, 和弦沒有
                         M(ii, jj, (kk - 1) * 12 + p)     = sum((sum(pitchClass(ii:jj-1, :),1) == 0) .* Template);        % 片段沒有, 和弦 有
                         S(ii, jj, (kk - 1) * 12 + p)     = ...
                             P(ii, jj, (kk - 1) * 12 + p) - (N(ii, jj, (kk - 1) * 12 + p) + M(ii, jj, (kk-1) * 12 + p));
                         %%%%% 加最低音weight %%%%%%
                         if  isLowWeight && min(lowestP(ii:jj-1, 1)) ~= Inf
                             lowestPitchIdx = ii + find(lowestP(ii:jj-1, 1) == min(lowestP(ii:jj-1, 1))) - 1;
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

        %% 分段：最佳路徑演算法 HarmAn
        % &&&&&&&&&&& .   downbeat  &&&&&&&&&&&&&&
        optimalPath = HarmAn(maxScoreMatrix, partitionNum);

    %     partitionNum      = size(pitchClass,1)+1;
    %     downbeat          = find(minSeg(:,1)==fix(minSeg(:,1)));
    %     optimalPath = HarmAn(maxScoreMatrix, partitionNum, [downbeat; partitionNum]);
        % &&&&&&&&&&& .   downbeat  &&&&&&&&&&&&&&

        % 紀錄 最佳路徑的onset
        segOnset = [minSeg(:,1)' - barOnset(j), minSeg(end,2) - barOnset(j)];

        Ans.optimalPath(j, 1:length(optimalPath)) = optimalPath;
        Ans.segOnset(j, 1:length(optimalPath)) = segOnset(optimalPath);

        %% 最佳路徑的和弦標記
        for k = 2:length(optimalPath)
            Ans.chord{j,k-1}       = scoreMatrixIdx(optimalPath(k-1),optimalPath(k));
            Ans.tieBreaking(j,k-1) = {''};
            Ans.score(j,k-1)       = maxScoreMatrix(optimalPath(k-1),optimalPath(k));

            %% 分數一樣的話
            highScoreChord = find(S(optimalPath(k - 1), optimalPath(k), :) == maxScoreMatrix(optimalPath(k - 1), optimalPath(k)));
            % step 1 : highest root scoreMatrix
            if numel(unique(highScoreChord)) ~= 1
                highScoreRootScore  = rootS(optimalPath(k-1), optimalPath(k), highScoreChord);
                maxRootIdx          = find(highScoreRootScore == max(highScoreRootScore));
                highScoreChord      = highScoreChord(maxRootIdx);
                step                = {'root'};
                % step 2 : highest probability of occurrence
                if numel(unique(highScoreChord))~=1
                    chordSpecies    = floor(highScoreChord / 12) + 1;
                    maxProbIdx      = find(chordSpecies == min(chordSpecies));
                    highScoreChord  = highScoreChord(maxProbIdx); 
                    step            = {'probability'};
                    if numel(unique(highScoreChord)) ~= 1        
                        highScoreChord = 72
                        step     = {'no solve'};
                        bar      = j
                    end
                end

                Ans.chord{j, k - 1} = highScoreChord;
                Ans.tieBreaking(j, k - 1) = step;
            end

            Ans.templateNo(j, k - 1) = ceil(Ans.chord{j, k-1}(1) / 12);
            Ans.pitchNo{j, k - 1} = ~(ceil(mod(Ans.chord{j, k - 1}, 12) / 12)) * 12 + mod(Ans.chord{j, k - 1}, 12);
            Ans.chordName{j, k - 1} = strcat(pitchName(Ans.pitchNo{j, k - 1}), ':', tempName(Ans.templateNo(j, k - 1)));

            evaluationIdx = evaluationIdx + 1;
            evaluationChord{evaluationIdx, 1} = j;
            evaluationChord{evaluationIdx, 2} = segOnset(optimalPath(k - 1));
            evaluationChord(evaluationIdx, 4) = Ans.chordName{j, k - 1};
            evaluationChord{evaluationIdx, 5} = highScoreChord;


        end

                %%
                %{

        % 算此小節的和弦在24調中 屬於該調的和弦數量
        for key_i=1:24
            ansBarChord = Ans.chord{j,:};
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

            keyRatio(j,:) = keyConsistRatio(noteBar);               % new
            keyRatioRank1(j,:) = keyRatio(j,:).*keyRank(j,:);       % new
            keySortRatio(j,:) = keyRatioRank1(j,nameIdxRank(j,:));  % new
            % 處理 把max留下來 其他不要
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
end

if isSave
    cell2csv(['chordEva/eva_' fileName '.csv'], evaluationChord);
end
%%
clc
GTFile  = 'trans_chord_k309_m1';%trans_chord_k309_m1, trans_chord_k545
[number, text,  GTdata] = xlsread(['chordGT/' GTFile '.xlsx']);
[score, keyRatio] = key_detection(evaluationChord, barNote);%evaluationChord
finalScore = score.chordScoreNormalize .* keyRatio;

function [score, keyRatio] = key_detection(evaluationChord, barNote)
    keyName     = {'C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B', ...
                   'c', 'c#', 'd', 'd#', 'e', 'f', 'f#', 'g', 'g#', 'a', 'a#', 'b'};
    evaluationChord(1, :) = [];
    keyChord = make_key_chord;
    for barI = 1:cell2mat(evaluationChord(end, 1))
        barIdx = cell2mat(evaluationChord(:, 1)) == barI;
        barChord = cell2mat(evaluationChord(barIdx, 5));
        for keyI = 1:24
            chordIntersect = intersect(barChord, keyChord(keyI, :));
            chordScore = 0;
            for c = 1:length(chordIntersect)
                chordScore = chordScore + numel(find(barChord == chordIntersect(c)));
            end
            score.chordScoreSum(barI, keyI) = chordScore;
            score.chordScoreNormalize(barI, keyI) = chordScore/numel(barChord);
        end
        keyRatio(barI, :) = keyConsistRatio(barNote{barI}); % 組成音
        
        finalScore = score.chordScoreNormalize(barI, :).* keyRatio(barI, :);
        scoreSort = sort(unique(finalScore), 'descend');
        if length(scoreSort) > 1 && scoreSort(1) ~= 0
        chordNameNo1 = keyName(finalScore == scoreSort(1))';
        score.chordNo1(1:length(chordNameNo1), barI) = chordNameNo1;
        end
        if length(scoreSort) > 2 && scoreSort(2) ~= 0
            chordNameNo2 = keyName(finalScore == scoreSort(2))';
            score.chordNo2(1:length(chordNameNo2), barI) = chordNameNo2;
        end
    end
end


%%
% description : 整理每個小節有哪些音，及音符資訊
% input : mididata -> midi資訊
%         timeSig  -> 拍號
% output : barNote -> 每個小節的midi資訊
%          barOnset-> 每個小節開始的beat
function [barNote, onsetBar] = bar_note_data(midiData, timeSig)
    addBeat     = 0; 
    addBar      = 0; 
    midiSize    = size(midiData, 2);
    
    for i = 1:size(timeSig, 1)
        time = timeSig(i, 1) * (4 / (2^timeSig(i, 2))); % 拍數
        if i ~= size(timeSig, 1)
            barNo = max(ceil( ((timeSig(i + 1, 5) - 1) - timeSig(i, 5) ) / time), 1);
        else
            barNo = max(ceil( ( sum(midiData(end, 1:2)) - timeSig(i, 5) ) / time), 1);
        end

        for j = 1:barNo
            barOnset = addBeat + (j - 1) * time;
            barOffset = addBeat + j * time;
            barNoteIdx = intersect(find(midiData(:, 1) >= barOnset), find(midiData(:, 1) < barOffset));

            currentBar = addBar + j;
            midiData(barNoteIdx, midiSize + 1) = currentBar;
            onsetBar(currentBar) = barOnset;

            % 前個小節圓滑線至此小節的音符slur
            offsetInBarIdx = intersect(find(sum(midiData(:, 1:2), 2) < barOffset), find(sum(midiData(:, 1:2), 2) > barOnset));
            onsetbeforeBarIdx = intersect(find(midiData(:, 1) < barOnset), offsetInBarIdx);
            slurNote = midiData(onsetbeforeBarIdx, :);
            slurNote(:, 2) = slurNote(:, 2) + barOnset - slurNote(:, 1);
            slurNote(:, 1) = barOnset; % 把onset改成此小節一開始

            noteBar = [midiData(barNoteIdx, :); slurNote]; % 此小節內的音符

            if ~isempty(noteBar)
                noteBar = sortrows(noteBar, 1);
                noteBar = trill_detection(noteBar); % tri 處理
                noteBar = normalize_midi_data(noteBar);

                barNote{currentBar, 1} = noteBar;
            end
        end
        
        addBar  = addBar  + barNo;
        addBeat = addBeat + time * barNo;
    end
end
%%
% description : 建立大小調的和弦
% output : 大小調的順階和弦
function [keyChord] = make_key_chord
    load Key_DiatonicTriad_V7.mat keyFusion

    % MAJOR    I IV V V7 ii vi viio iii
    tempMaj = [1 1 1 2 3 3 4 3];         % 有3
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
% description: 以拍點為partition point，紀錄每個minSeg內的音符權重，以及minSeg內最低音
% input:  noteBar      -> 小節音符
%         barOnset     -> 小節開始的beat
% output: pitchClass,  -> 每個minSeg有所有音，且weight為音符長度，大小為N*12 (N為minSeg數量,12為一個八度音的數量)
%         lowestPitchClass -> 每個minSeg的最低音，且weight為音符長度，大小為N*12 (N為minSeg數量,12為一個八度音的數量)
%         lowestP      -> 每個minSeg的最低音，且不以一個八度來看(不mod12)
%         minSeg       -> 每個minSeg之onset offset的beat
% 找左手最低音權重較大，難以分辨左右手，所以找片段最低音且同時有兩個音
function [pitchClass, lowestPitchClass, lowestP, minSeg] = downbeat_pitchclass(noteBar, barOnset, isDownbeatWeight)
    beatInBar = noteBar(:,9);
    minSeg = [0:beatInBar-1; 1:beatInBar]'; % 以拍點為最小片段 (第 0 1 2 3 拍)
%     minSeg = [0:0.5:beatInBar; 0:0.5:beatInBar]; % 以拍點為最小片段 (第 0 0.5 1 1.5 2 2.5 3 3.5 拍)
%     minSeg([2,17])=[]; minSeg = reshape(minSeg,2,length(minSeg)/2)';
    
    pitchClass = zeros(size(minSeg, 1), 12);
    lowestPitchClass = zeros(size(minSeg, 1), 12);

    lowestP = inf * ones(size(minSeg, 1), 1);
    
    for s = 1:size(minSeg, 1)
        unit = (1/2)^4;
        noteOn = zeros(size(noteBar,1),1);
        noteOff = zeros(size(noteBar,1),1);
        isOnsetInt = zeros(size(noteBar, 1), 1);
        noteNum = zeros(ceil((minSeg(s, 2) - minSeg(s, 1)) / unit), 1); % 紀錄minSeg內每unit長度，有幾個音
        noteLowestIdx = zeros(ceil((minSeg(s, 2) - minSeg(s, 1)) / unit), 1); % 每個unit最低音是哪個音
        noteLowestPitch = inf * ones(ceil((minSeg(s, 2) - minSeg(s, 1)) / unit), 1); % 每個unit最低音的音高

        for n = 1:size(noteBar, 1)
            noteOn(n) = noteBar(n, 1) - barOnset;
            noteOff(n) = sum(noteBar(n,1:2)) - barOnset;
            isOnsetInt(n) = (noteOn(n) == fix(noteOn(n))); % isDownbeatWeight
            
            if noteOn(n) < minSeg(s, 1); noteOn(n) = minSeg(s, 1); end
            if noteOff(n) > minSeg(s, 2); noteOff(n) = minSeg(s, 2); end
            
            if  noteOff(n) - noteOn(n) > 0   % 判斷落在minSeg內的音
                pitchNo = mod(noteBar(n, 4), 12) + 1;
                pitchClass(s, pitchNo) = pitchClass(s, pitchNo) + noteOff(n) - noteOn(n);
                
                if isDownbeatWeight && isOnsetInt(n)
                    pitchClass(s, pitchNo) = pitchClass(s, pitchNo) + noteOff(n) - noteOn(n);
                end

                onUnitNum  = floor((noteOn(n) - minSeg(s,1)) / unit) + 1; % 音on和off是佔minSeg內的哪些unit
                offUnitNum = ceil((noteOff(n) - minSeg(s,1)) / unit);
                repPitch = ones(offUnitNum - onUnitNum + 1, 1) * noteBar(n, 4);

                [minPitch, idx] = min([noteLowestPitch(onUnitNum:offUnitNum) repPitch], [], 2);
                noteNum(onUnitNum:offUnitNum) = noteNum(onUnitNum:offUnitNum) + 1;

                if ~isempty(find(idx == 2, 1))
                    noteLowestIdx(find(idx == 2) + onUnitNum - 1) = n;
                    noteLowestPitch(onUnitNum:offUnitNum) = minPitch;
                end
            end
            
            %% 兩個音同時，取C4以下的音
            if n == size(noteBar, 1)
                if ~isempty(find(noteNum == 2, 2)) %%  兩個音同時取最低音,c4以下
                    lowestP(s) = min(noteLowestPitch(noteNum == 2));
                    if lowestP(s) < 60
                        idx = find(noteLowestPitch == lowestP(s)); %要改
                        lowestDur = numel(idx);
                        if isDownbeatWeight
                            weight = 0;
                            intSect = intersect(noteLowestIdx, find(isOnsetInt == 1));
                            for inS = 1:length(intSect)
                                weight = weight + numel(find(noteLowestIdx(idx) == intSect(inS)));
                            end
                            lowestDur = lowestDur + 0.5*weight;
                        end
                        lowestPitchClass(s, mod(lowestP(s), 12) + 1) = lowestDur * unit;
                    end
                end
                if ~isempty(find(noteNum > 2, 2)) %%  三個音同時取最低音，優先
                    lowestP(s) = min(noteLowestPitch(noteNum > 2));
                    idx = find(noteLowestPitch == lowestP(s));
                    lowestDur = numel(idx);
                    if isDownbeatWeight
                        weight = 0;
                        intSect = intersect(noteLowestIdx, find(isOnsetInt == 1));
                        for inS = 1:length(intSect)
                            weight = weight + numel(find(noteLowestIdx(idx) == intSect(inS)));
                        end
                        lowestDur = lowestDur + 0.5*weight;
                    end
                    lowestPitchClass(s, mod(lowestP(s), 12) + 1) = lowestDur * unit;
                end 
            end
            %%  兩個音同時取最低音
%             if n == size(noteBar, 1) && ~isempty(find(noteNum > 1, 1))
%                 lowestP(s) = min(noteLowestPitch(noteNum > 1));
%                 idx = find(noteLowestPitch == lowestP(s));
%                 lowestDur = numel(idx);
%                 if isDownbeatWeight
%                     weight = 0;
%                     intSect = intersect(noteLowestIdx, find(isOnsetInt == 1));
%                     for inS = 1:length(intSect)
%                         weight = weight + numel(find(noteLowestIdx(idx) == intSect(inS)));
%                     end
%                     lowestDur = lowestDur + weight;
%                 end
%                 lowestPitchClass(s, mod(lowestP(s), 12) + 1) = lowestDur * unit;
%             end
        end     
    end

    minSeg = minSeg + barOnset;
end

%%
% description: 所有片段的pitchClass
% input     : noteBar -> 小節內的 midi data
%             barOnset-> 小節開始的beat
% output    : pitchClass -> 紀錄小節所有片段的pitch class
%             minSeg -> 小節所有最小片段onset offset
function [pitchClass, minSeg] = min_seg_pitchclass(noteBar, barOnset)

    % downbeat 為最小片段
%     beatInBar = noteBar(:, 9);
% %     minSeg = [0:beatInBar-1; 1:beatInBar]';
%     minSeg = unique(round([noteBar(:, 1); sum(noteBar(:, 1:2), 2); (0:beatInBar-1)' + barOnset], 4));
%     minSeg(:, 2) = [minSeg(2:end,1); beatInBar(1) + barOnset];

    % 最小片段的onset offset
    minSeg = unique(round([noteBar(:, 1); sum(noteBar(:, 1:2), 2)], 4));
    minSeg(:, 2) = [minSeg(2:end, 1); round(sum(noteBar(end, 1:2)), 4)];
    
    % 把minSeg太小的去掉
    miniTime = 0.025;
    if minSeg(end, 1) == minSeg(end, 2); minSeg(end, :)=[]; end
    minSeg(diff(minSeg, 1, 2) <= miniTime, :) = [];

    % 紀錄每個最小片段有哪些音
    pitchClass = zeros(length(minSeg), 12);
    for k=1:size(noteBar, 1)
        s = find(abs(round(noteBar(k, 1), 4) - minSeg(:, 1)) <= miniTime); % 這個音從哪個seg開始
        e = find(round(sum(noteBar(k, 1:2)), 4) >= minSeg(:, 2)); % 這個音從哪個seg結束
        pitchNum = mod(noteBar(k, 4), 12);
        pitchClass(s:e(end), pitchNum + 1) = pitchClass(s:e(end), pitchNum + 1) + 1;
    end
end

%%
function [maxScoreMatrix, I, S] = record_min_seg_score(pitchClass)
    %                            大三        屬七         小三       減
    template            = [ 0,4,7,-1;  0,4,7,10;   0,3,7,-1; 0,3,6,-1];  
    templateNum             = size(template,1);
    partitionNum      = size(pitchClass,1)+1;
     rootS  = -inf*ones(partitionNum, partitionNum, templateNum*12);
     s      = -inf*ones(partitionNum, partitionNum, templateNum*12);
     N      = -inf*ones(partitionNum, partitionNum, templateNum*12);
     M      = -inf*ones(partitionNum, partitionNum, templateNum*12);
     S      = -inf*ones(partitionNum, partitionNum, templateNum*12);

     for ii=1:partitionNum
         for jj=ii+1:partitionNum
            for kk=1:templateNum
                for p = 1:12
                    templateNo = mod(template(kk,(template(kk,:)~=-1)) + p - 1 , 12);
                    Template = zeros(1,12); Template(1,templateNo + 1) = 1;
                    rootS(ii,jj,(kk-1)*12+p) = sum(sum(pitchClass(ii:jj-1,templateNo(1)+1))); % 根音的score
                    s(ii,jj,(kk-1)*12+p)     = sum(sum(pitchClass(ii:jj-1,templateNo+1  )));                  % 片段 有, 和弦 有
                    N(ii,jj,(kk-1)*12+p)     = sum(sum(pitchClass(ii:jj-1,:)))-s(ii,jj,(kk-1)*12+p);    % 片段 有, 和弦沒有            
                    M(ii,jj,(kk-1)*12+p)     = sum((sum(pitchClass(ii:jj-1,:),1)==0).*Template);        % 片段沒有, 和弦 有
                    S(ii,jj,(kk-1)*12+p)     = s(ii,jj,(kk-1)*12+p) - ( N(ii,jj,(kk-1)*12+p) + M(ii,jj,(kk-1)*12+p) );
                end
            end
         end
     end

    [maxScoreMatrix,I] = max(S,[],3);
end

%%
% description: 計算每個小節所有的大小調組成音比例
% input     : note -> midi data
% output    : keyRatio -> 所有調組成音比例
function [keyRatio] = keyConsistRatio(note)
    load Key_DiatonicTriad_V7.mat keyFusion

    pitchConsist = zeros(1, 12);
    keyConsist = zeros(24, 12);
    keyRatio = zeros(1, 12);
    
    barPitch = mod(note(:,4), 12) + 1;
    pitchTabulate = tabulate(barPitch);
    pitchConsist(1:size(pitchTabulate, 1)) = pitchTabulate(:, 2);
    
%     for i=1:24
%         keyConsist(i, keyFusion(i).consistNumber(:) + 1) = 1;
%         keyRatio(i) = sum(keyConsist(i,:) .* pitchConsist) / size(note, 1);
%     end
    %% new
    MAX = 0;
     for k=1:24
         note1 = note; 
         for j=1:7
             note1(mod(note1(:,4),12)==keyFusion(k).consistNumber(j),:)=[];
         end
         keyRatio(k) = (size(note,1) - length(unique(note1(:,4))))/size(note,1);
%          keyRatio(k) = size(note1,1);
%          keyRatio(k) = length(unique(note1(:,4)));
%          MAX = max(MAX, keyRatio(k));
%                  barBARbar{i,k} = note1{i,1};
     end
%     keyRatio = (MAX-keyRatio)/MAX;
end
