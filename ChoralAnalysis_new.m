 % implement : Algorithms for Choral Analysis 2002
    addpath('toolbox/midi_lib/midi_lib');
    clear; close all; clc;

    tempName = {'maj','7','min','dim','xxx','X'};%,'Fully dim7','Half dim7','dim3','X'};
    pitchName = {'C','C#','D','D#','E','F','F#','G','G#','A','A#','B'};
    keyName = {'C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B', ...
           'c', 'c#', 'd', 'd#', 'e', 'f', 'f#', 'g', 'g#', 'a', 'a#', 'b'};
 
%../midi/
    filename = 'mz_545_1_noRepeat';%beethoven_string_quartet_135, mz_16_1 O48N2AOM quartet_16_3 dim7_1 
    [midiINFO, timeSig]  = midi_Preprocess(filename);
    
    [Note, barStart] = BarInfo(midiINFO, timeSig); % each bar's information
    unit = 0.25; 
    twoNeipoint = [floor(midiINFO(:,2)/unit)*unit, ceil(midiINFO(:,2)/unit)*unit];
    [~, I] = min(abs(twoNeipoint-midiINFO(:,2)),[],2);
    k = twoNeipoint(:,I);
    
    
    [keyCHORD] = keyChord;

    for j = 1:length(Note)

    noteBar = Note{j};
    
    %% 找到最小片段分割點有哪些音

            % 最小片段的onset offset
            miniSeg       = unique(round([noteBar(:,1); sum(noteBar(:,1:2),2)],4));
            miniSeg(:,2)  = [ miniSeg(2:end,1); round(sum(noteBar(end,1:2)),4)];

            if miniSeg(end,1)==miniSeg(end,2); miniSeg(end,:)=[]; end                      % 把miniSeg太小的去掉
            miniTime      = 0.025;
            miniSeg(diff(miniSeg,1,2)<=miniTime,:) = [];

            miniSeg2Onset = tabulate(round(noteBar(:,1),4));
            Onset2        = find(miniSeg2Onset(:,2)>=2);

            % 紀錄每個最小片段有哪些音
            pitchClass    = zeros(length(miniSeg),12);
            parpointNum   = zeros(length(miniSeg), 1);

            for k=1:size(noteBar, 1)
                s               = find(abs(round(noteBar(k,1),4) - miniSeg(:,1)) <= miniTime);   % 這個音從哪個seg開始
                e               = find(round(sum(noteBar(k,1:2)),4)>=miniSeg(:,2));              % 這個音從哪個seg結束
                parpointNum(s)  = parpointNum(s) + 1;
                pitchNum        = mod(noteBar(k,4),12);
                pitchClass(s:e(end),pitchNum+1) = pitchClass(s:e(end),pitchNum+1)+1;
            end

    %% 紀錄所有可能分段的分數 
                                    % 大三        屬七         小三      減
            template            = [ 0,4,7,-1;  0,4,7,10;   0,3,7,-1; 0,3,6,-1];  

            tempNUM             = size(template,1);
            partitionPoint      = length(miniSeg)+1;

             P      = -inf*ones(partitionPoint, partitionPoint, tempNUM*12);
             N      = -inf*ones(partitionPoint, partitionPoint, tempNUM*12);
             M      = -inf*ones(partitionPoint, partitionPoint, tempNUM*12);
             S      = -inf*ones(partitionPoint, partitionPoint, tempNUM*12);

             for ii=1:partitionPoint
                 for jj=ii+1:partitionPoint
                    for kk=1:tempNUM
                        for p = 1:12
                            temp = mod(template(kk,(template(kk,:)~=-1)) + p - 1 , 12);
                            Template = zeros(1,12); Template(1,temp + 1) = 1;

                            rootS(ii,jj,(kk-1)*12+p) = sum(sum(pitchClass(ii:jj-1,temp(1)+1))); % 根音的score
                            P(ii,jj,(kk-1)*12+p)     = sum(sum(pitchClass(ii:jj-1,temp+1  )));                  % 片段 有, 和弦 有
                            N(ii,jj,(kk-1)*12+p)     = sum(sum(pitchClass(ii:jj-1,:)))-P(ii,jj,(kk-1)*12+p);    % 片段 有, 和弦沒有            
                            M(ii,jj,(kk-1)*12+p)     = sum((sum(pitchClass(ii:jj-1,:),1)==0).*Template);        % 片段沒有, 和弦 有
                            S(ii,jj,(kk-1)*12+p)     = P(ii,jj,(kk-1)*12+p) - ( N(ii,jj,(kk-1)*12+p) + M(ii,jj,(kk-1)*12+p) );
                        end
                    end
                 end
             end

            [weight,I] = max(S,[],3);


    %% 分段：最佳路徑演算法 HarmAn
   %             MARK = HarmAn( [1;Onset2;partitionPoint], partitionPoint, weight );
            MARK                        = HarmAn( 1:partitionPoint,partitionPoint, weight );
            Ans.Mark(j,1:length(MARK))  = MARK;

            % 紀錄 最佳路徑的onset
            segOnset = miniSeg(:,1)' - barStart(j);
            segOnset = [segOnset miniSeg(end,2)-barStart(j)];
            Ans.segOnset(j,1:length(MARK)) = segOnset(MARK);

            Dim7Time                    = 1;

    %% 最佳路徑的和弦標記
            Wtmp = 0;
            Stmp = 0;
            for k = 2:length(MARK)
                Ans.Chord(j,k-1)       = I(MARK(k-1),MARK(k));
                Ans.sameScore(j,k-1)   = -1;
                Ans.Score(j,k-1)       = weight(MARK(k-1),MARK(k));

                Wtmp                   = Wtmp + weight(MARK(k-1),MARK(k));
                
                Stmp                   = Stmp + weight(MARK(k-1),MARK(k))*diff(miniSeg(k-1,:))/(miniSeg(end,end)-miniSeg(1,1));
                if k == length(MARK)
                    Ans.ScoreSum(j,1) = Wtmp;
                    Ans.ScorePar(j,1) = Wtmp/partitionPoint;
                    Ans.ScoreLen(j,1) = Stmp;
                    Ans.ScoreLenPar(j,1) = Ans.ScoreLen(j,1)/partitionPoint;
                end


        %% 分數一樣的話
                HighScoreChord = find(S(MARK(k-1),MARK(k),:)==weight(MARK(k-1),MARK(k)));

                if numel(unique(HighScoreChord))~=1
        % step 1 : highest root weight
                    HighScoreRootScore = rootS(MARK(k-1),MARK(k),HighScoreChord);
                    MaxRootIdx = find(HighScoreRootScore==max(HighScoreRootScore));

                    if numel(unique(MaxRootIdx))==1
                        AnsChord = HighScoreChord(MaxRootIdx);
                        step     = 1;
                    else
        % step 2 : highest probability of occurrence
                        chordSpecies = floor(HighScoreChord/12)+1;
                        MaxProbIdx   = find(chordSpecies==min(chordSpecies));

                        if numel(unique(MaxProbIdx))==1
                            AnsChord = HighScoreChord(MaxProbIdx);
                            step     = 2;
                        else          
                            HighScoreChord(MaxProbIdx)
                            AnsChord = size(tempName,2)*12;%97; %Ans.Chord(j,k-1);
                            step     = 4
                            bar      = j
                        end
                    end

                    Ans.Chord(j,k-1) = AnsChord;
                    Ans.sameScore(j,k-1) = step;
                end

                Ans.templ_no(j,k-1)  = ceil(Ans.Chord(j,k-1)/12);
                Ans.pitch_no(j,k-1)  = ~(ceil(mod(Ans.Chord(j,k-1),12)/12))*12 + mod(Ans.Chord(j,k-1),12);
                Ans.ChordName{j,k-1} = strcat(pitchName(Ans.pitch_no(j,k-1)),':',tempName(Ans.templ_no(j,k-1)));

            end
            
            % 算此小節的和弦在24調中 屬於該調的和弦數量
            for key_i=1:24
                ansBarChord = Ans.Chord(j,:);
                DIFF = numel(setdiff(ansBarChord(ansBarChord~=0),keyCHORD(key_i,:)));
                num = numel(ansBarChord(ansBarChord~=0));
                SCORE(j,key_i) = (num-DIFF)/num;
                if key_i==2 && j ==8
                    setdiff(ansBarChord(ansBarChord~=0),keyCHORD(key_i,:))
                    ansBarChord
                    keyCHORD(key_i,:)
                    DIFF
                    num
                   
                end
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
    end
    %%
%     [aaaaaaa] = postProcessing(keyRank);
    
%     function [a] = postProcessing(keyRank)
%     a=[];
%     rankmat = keyRank;
%     
%         for i=1:size(keyRank,2)
%             zero = find(keyRank(:,i)==0);
%             for j=1:length(zero)
%                 if zero(j)>2 && zero(j)<size(keyRank,1)
%                     rankmat(zero(j));
%                 end
%             end
%         end
%     end

%     imagePre = SCORE';
%     imagePre(imagePre~=1) = 0;
%     imshow(imagePre, 'InitialMagnification', 800)
%     pause;
%     imagePre = SCORE';
%     imagePre(imagePre~=0) = 1;
%     imshow(imagePre, 'InitialMagnification', 800)
    
%%
% for i=1:length(Note)
%     KeyConsistRatio(Note{58})
% end
    
function [keyRatio] = KeyConsistRatio(Note)
% 小節內 每個調的組成音比例
    load    Key_DiatonicTriad_V7.mat

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

function [Note, barSTART] = BarInfo(midiINFO, timeSig)
    
    addBeat     = 0; 
    addBar      = 0; 
    midiSize    = size(midiINFO,2);
    

    for i = 1:size(timeSig,1)

        NoOfBeats       =   timeSig(i,1)*(4/(2^timeSig(i,2)));
        if i ~= size(timeSig,1)
            NoOfBar     =   max(ceil( ((timeSig(i+1,5)-1) - timeSig(i,5) ) / NoOfBeats), 1);
        else
            NoOfBar     =   max(ceil( ( sum(midiINFO(end,1:2)) - timeSig(i,5) ) / NoOfBeats), 1);
        end

        for j = 1:NoOfBar
            barStart   = (addBeat+(j-1)*NoOfBeats);
            barEnd     = addBeat+j*NoOfBeats;
            notesOfBar = intersect(find(midiINFO(:,1)>=barStart), find(midiINFO(:,1)<barEnd));

            currentBar = addBar + j;
            barSTART(currentBar) = barStart;
            midiINFO(notesOfBar, midiSize+1) = currentBar;

            % 前個小節圓滑線至此小節的音符
            innerOffset     = intersect(find(sum(midiINFO(:,1:2),2)<barEnd) , find(sum(midiINFO(:,1:2),2)>barStart));
            overOnset       = intersect(find(midiINFO(:,1)<barStart), innerOffset);
            overNote        = midiINFO(overOnset, :);
            overNote(:,2)   = overNote(:,2) + barStart - overNote(:,1);
            overNote(:,1)   = barStart;

            noteBar         = midiINFO(notesOfBar, :);      % 此小節內的音符
            noteBar         = [noteBar; overNote];          % 把onset改成此小節一開始

            if ~isempty(noteBar)

                noteBar       = sortrows(noteBar,1);        % 依據第一行排序 小到大
                tmp = length(noteBar);
                noteBar       = triDetection(noteBar);      % tri 處理
    %                 noteBar(noteBar(:,2)<0.2,:) = [];         % 把小於duration 0.25的音符都刪掉，認為太小音符不會對和弦有太大的影響。
                Note{currentBar,1}      = noteBar;
            end
        end

        addBar  = addBar  + NoOfBar;
        addBeat = addBeat + NoOfBeats*NoOfBar;

    end

end


function [keyCHORD] = keyChord
    load    Key_DiatonicTriad_V7.mat

    tempName = {'maj','7','min','dim','xxx','X'};%,'Fully dim7','Half dim7','dim3','X'};
    pitchName = {'C','C#','D','D#','E','F','F#','G','G#','A','A#','B'};
    
    % major   I IV V V7 ii vi viio iii
    tempNo = [1 1 1 2 3 3 4 3];  % 有3
    for i=1:12
        keyCHORD(i,:) = (tempNo-1)*12 + keyFusion(i).diatonicTriad(1,[1 4 5 5 2 6 7 3])+1;
    end
    % minor   i iv V V7 iio VI viio iii+
    tempNo = [3 3 1 2 4 1 4 1];
    for i=13:24
        keyCHORD(i,:) = (tempNo-1)*12 + keyFusion(i).diatonicTriad(1,[1 4 5 5 2 6 7 3])+1;
    end
end