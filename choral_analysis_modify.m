%%%% motify : "Algorithms for Choral Analysis 2002"
%%%% date   : 18/01/22
%%%% content: 和弦分析演算法，以小節為單位進行分析。
%%%% input  : filename
%%%% output : chords in each bar

% description: 修改 "Algorithms for Choral Analysis 2002":和弦分析演算法，以小節為單位進行分析
% date   : 18/01/22
% input     : filename -> midi檔名(必須把音檔放在midi資料夾中)
%             parameter-> 一些條件的參數
%                 isDownbeatAndDurWeight : 以拍點為切割點，音符長度為權重
%                 isLowWeight : 左手最低音權重調大(兩倍)
%                 isDownbeatWeight : 拍點上的音權重調大(兩倍)
%             isSave-> 和弦結果是否存檔
% output    : evaluationChord -> 和弦結果
function evaluationChord = choral_analysis_modify(filename, parameter, isSave)

    if nargin < 2, parameter = []; end
    if nargin < 3, isSave = 0; end
    if isfield(parameter,'isDownbeatAndDurWeight')==0
        parameter.isDownbeatAndDurWeight = 1;
    end
    if isfield(parameter,'isLowWeight')==0
        parameter.isLowWeight = 0;
    end
    if isfield(parameter,'isDownbeatWeight')==0
        parameter.isDownbeatWeight = 0;
    end
   
    % include toolbox
    addpath('toolbox/midi_lib/midi_lib');

% % input midi data
% filePath = '../midi';
% filename = 'h_37_cut';%m_7_cut, mz_545_1_noRepeat, 
    [midiData, timeSig]  = midi_Preprocess(filename);

    % each bar's information
    [barNote, barOnset] = bar_note_data(midiData, timeSig);

    % parameter
    tempName    = {'maj','7','min','dim','xxx','X'};%,'Fully dim7','Half dim7','dim3','X'}; chordType
    pitchName   = {'C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B'};
    keyName     = {'C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B', ...
                   'c', 'c#', 'd', 'd#', 'e', 'f', 'f#', 'g', 'g#', 'a', 'a#', 'b'};
    %                 大三        屬七         小三       減
    template    = [ 0,4,7,-1; 0,4,7,10; 0,3,7,-1; 0,3,6,-1];  
    evaluationIdx = 1;
    evaluationChord = {'小節','拍數(onset)','調性','和弦名稱','和弦編號','備註'};% new




    for j = 1:length(barNote)
        noteBar = barNote{j};
        if ~isempty(noteBar)
            % 找到最小片段分割點有哪些音
            [pitchClass, minSeg] = min_seg_pitchclass(noteBar, barOnset(j));   
            if parameter.isDownbeatAndDurWeight
                [pitchClass, lowestPClass, lowestP, minSeg] = downbeat_pitchclass(noteBar, barOnset(j), parameter.isDownbeatWeight);
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
                             M(ii, jj, (kk - 1) * 12 + p)     = sum((sum(pitchClass(ii:jj-1, :),1) == 0) .* Template);        % 片段沒有, 和弦 有
                             S(ii, jj, (kk - 1) * 12 + p)     = ...
                                 P(ii, jj, (kk - 1) * 12 + p) - (N(ii, jj, (kk - 1) * 12 + p) + M(ii, jj, (kk-1) * 12 + p));
                             %%%%% 加最低音weight %%%%%%
                             if  parameter.isLowWeight && min(lowestP(ii:jj-1, 1)) ~= Inf
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
            optimalPath = HarmAn(maxScoreMatrix, partitionNum);

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
        end 
    end

    if isSave
        cell2csv(['chordEva/eva_' filename '.csv'], evaluationChord);
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


