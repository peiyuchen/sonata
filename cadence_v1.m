% implement : Algorithms for Choral Analysis 2002
 
    clear all; close all; clc;

    template  = [ 0,4,7,-1;  0,4,7,10;   0,3,7,-1;  0,3,6, 9;  0,3,6,10;   0,3,6,-1 ];
    tempName  = {'maj','7','min','減七°7','半減七?','°','ERROR'};
    pitchName = {'C :','C#:','D :','D#:','E :','F :','F#:','G :','G#:','A :','A#:','B :'};
    %bee_14_1
    filename             = 'bee_14_1';%Hob._XVI15_2nd_movement';% mz_545_1_exposition, beethoven_string_quartet_135, mz_16_1, quartet_16_3, dim7_1, bee_14_1
    [midiINFO, timeSig]  = midi_Preprocess(['../midi/' filename]);
    roundUnit = 4;
    midiINFO(:,1:2)      = round(midiINFO(:,1:2),roundUnit);
    midiSize             = size(midiINFO,2);
    addBeat = 0; addBar = 0; refNextBarRoot = 0; Untreated = []; markLen = 0; markSec = []; Tri_bar = [];
    

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

                noteBar       = sortrows(noteBar,1);      % 依據第一行排序 小到大
                triBefore     = size(noteBar,1);
                noteBar       = triDetection(noteBar);    % tri 處理
                triAfter      = size(noteBar,1);
                if triBefore~=triAfter; Tri_bar = [Tri_bar currentBar]; end
                noteBar(noteBar(:,2)<0.2,:) = [];             % 把小於duration 0.25的音符都刪掉，認為太小音符不會對和弦有太大的影響。

    %% 找到最小片段分割點有哪些音

                % 最小片段的onset offset
%                 miniSeg       = 
                [miniSeg,ia1,~]     = unique(round([noteBar(:,1); sum(noteBar(:,1:2),2)],roundUnit));
                miniSeg(:,2)        = [ miniSeg(2:end,1); round(sum(noteBar(end,1:2)),roundUnit)];
                    % % 單位秒的也要做
                    allSeg_Sec       = round([noteBar(:,6); sum(noteBar(:,6:7),2)],roundUnit);
                    miniSeg_Sec      = allSeg_Sec(ia1);
%                     miniSeg_Sec      = unique(round([noteBar(:,6); sum(noteBar(:,6:7),2)],roundUnit));
                    miniSeg_Sec(:,2) = [ miniSeg_Sec(2:end,1); round(sum(noteBar(end,6:7)),roundUnit)];
                   
                % 把miniSeg太小的去掉
                if miniSeg(end,1)==miniSeg(end,2); miniSeg(end,:)=[]; miniSeg_Sec(end,:)=[];end
                
                miniTime      = 0.035;
                    miniSeg_Sec(diff(miniSeg,1,2)<=miniTime,:) = []; % % 單位秒的也要做
                miniSeg(diff(miniSeg,1,2)<=miniTime,:) = [];
                    
                    
                
                % 一次有兩個音才考慮它可能是路徑
                miniSeg2Onset = tabulate(round(noteBar(:,1),roundUnit));
                Onset2        = find(miniSeg2Onset(:,2)>=2);

                % 紀錄每個最小片段有哪些音
                pitchClass    = zeros(length(miniSeg),12);
                pitchCell{currentBar,size(miniSeg, 1)}     = [];

                for k=1:size(noteBar, 1)
                    s               = find(abs(round(noteBar(k,1),roundUnit) - miniSeg(:,1)) <= miniTime);   % 這個音從哪個seg開始
                    e               = find(round(sum(noteBar(k,1:2)),roundUnit)>=miniSeg(:,2));              % 這個音從哪個seg結束
                    pitchNum        = mod(noteBar(k,4),12);
%                     pitchCell{currentBar,s:e(end)} = [pitchCell{currentBar,s:e(end)} noteBar(k,4)];
                    for t=s:e(end)
                        pitchCell{currentBar,t} = [pitchCell{currentBar,t} noteBar(k,4)];
                    end
%                     pitchCell{currentBar,s:e(end)} = cat(2, pitchCell{currentBar,s:e(end)}, noteBar(k,4));
                    pitchClass(s:e(end),pitchNum+1) = pitchClass(s:e(end),pitchNum+1)+1;
                end

    %% 紀錄所有可能分段的分數 

                tempNUM             = size(template,1);
                partitionPoint      = size(miniSeg,1)+1;

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
                MARK     = HarmAn( 1:partitionPoint,partitionPoint, weight );
                Ans.Mark(currentBar,1:length(MARK))  = MARK;
                markLen(currentBar) = length(MARK)-1;
                markSec  = [markSec; miniSeg_Sec(MARK(1:end-1),1)]; % mark的sec
                % 紀錄 最佳路徑的onset
                segOnset = miniSeg(:,1)' - barStart;
                Ans.segOnset(currentBar,1:length(MARK)-1)  = segOnset(MARK(1:end-1));
                Ans.segOnset1(currentBar,1:length(MARK)-1) = miniSeg(MARK(1:end-1))';
                

    %% 最佳路徑的和弦標記
                % initial
                Dim7Time = 1;

                for k = 2:length(MARK)
                    Ans.Chord(currentBar,k-1)       = I(MARK(k-1),MARK(k));
                    Ans.sameScore(currentBar,k-1)   = -1;
                    Ans.Score(currentBar,1)         = weight(MARK(k-1),MARK(k));
            %% 最佳路徑pitch，PAC用
                    Ans.segPitch{currentBar,k-1}      = unique(cat(2,pitchCell{currentBar,MARK(k-1):MARK(k)-1}));
                    Ans.segPitchFirst{currentBar,k-1} = unique(pitchCell{currentBar,MARK(k-1)});
                    
            %% 分數一樣的話
                    HighScoreChord = find(S(MARK(k-1),MARK(k),:)==weight(MARK(k-1),MARK(k)));
currentBar
                    if numel(unique(HighScoreChord))~=1
            % step 1 : highest root weight
            HighScoreChord
                        HighScoreRootScore = rootS(MARK(k-1),MARK(k),HighScoreChord)
                        MaxRootIdx = find(HighScoreRootScore==max(HighScoreRootScore))

                        if numel(unique(MaxRootIdx))==1
                            AnsChord = HighScoreChord(MaxRootIdx);
                            step     = 1
                        else
            % step 2 : highest probability of occurrence
                            chordSpecies = floor(HighScoreChord(MaxRootIdx)/12)+1;
                            MaxProbIdx   = find(chordSpecies==min(chordSpecies));

                            if numel(unique(MaxProbIdx))==1
                                AnsChord = HighScoreChord(MaxProbIdx);
                                step     = 2
                            else                        
            % step 3 : Diminished 7th resolution
                                if ceil(HighScoreChord(MaxProbIdx)/12) == 4
                                    Resolution_Dim7.idx(currentBar,Dim7Time)  = k-1;
                                    Resolution_Dim7.root{currentBar,k-1}      = HighScoreChord;
                                    
                                    step     = 3
                                    Dim7Time = Dim7Time + 1;
                                    AnsChord = Ans.Chord(currentBar,k-1);
                                else
                                    AnsChord = HighScoreChord(MaxProbIdx);
%                                     AnsChord = 73; %Ans.Chord(currentBar,k-1);

                                    step     = 4
                                    Untreated = [Untreated currentBar];
                                    warning(['第 ' num2str(currentBar) '小節 - 即使使用規則還是有相同分數的template']);
                                end
                            end
                        end
                        %^^
                            if ceil(HighScoreChord/12) == 4
                                Resolution_Dim7.idx(currentBar,Dim7Time)  = k-1;
                                Resolution_Dim7.root{currentBar,k-1}      = HighScoreChord;

                                step     = 3
                                Dim7Time = Dim7Time + 1;
                                AnsChord = Ans.Chord(currentBar,k-1);
                            else
                                AnsChord = HighScoreChord(MaxProbIdx);
    %                                     AnsChord = 73; %Ans.Chord(currentBar,k-1);

                                step     = 4
                                Untreated = [Untreated currentBar];
                                warning(['第 ' num2str(currentBar) '小節 - 即使使用規則還是有相同分數的template']);
                            end
                            %^^

                        Ans.Chord(currentBar,k-1) = AnsChord;
                        Ans.sameScore(currentBar,k-1) = step;
                    end

                    Ans.templNo(currentBar,k-1)   = ceil(Ans.Chord(currentBar,k-1)/12);
                    Ans.pitchNo(currentBar,k-1)   = ~(ceil(mod(Ans.Chord(currentBar,k-1),12)/12))*12 + mod(Ans.Chord(currentBar,k-1),12);
                    Ans.ChordName{currentBar,k-1} = strcat(pitchName(Ans.pitchNo(currentBar,k-1)), tempName(Ans.templNo(currentBar,k-1)));

                    if refNextBarRoot
                        tmp = (Resolution_Dim7.idx(currentBar-1)~=0);
                        idx = Resolution_Dim7.idx(currentBar-1,end);
%                         idx               = find(Resolution_Dim7.idx(currentBar-1)~=0);
                        nowRoot           = ~(ceil(mod(Resolution_Dim7.root{currentBar-1,idx(end)},12)/12))*12 + mod(Resolution_Dim7.root{currentBar-1,idx(end)},12);
                        nextRoot          = Ans.pitchNo(currentBar,k-1);
                        reRoot            = ~(ceil(mod((nextRoot-1),12)/12))*12 + mod((nextRoot-1),12);
                        refNextBarRoot    = ~refNextBarRoot;
                        %^^
                        rootShift    = [nowRoot, circshift(nowRoot,3), circshift(nowRoot,2), circshift(nowRoot,1)]  %^^
                        nextRepresentation = (mod(template(Ans.templNo(currentBar,k-1),(template(Ans.templNo(currentBar,k-1),:)~=-1)) + Ans.pitchNo(currentBar,k-1)-1,12)+1)';
                        len = length(nextRepresentation);
                        if length(nextRepresentation)==3; nextRepresentation(4) = NaN; len = 3; end
                        nextRepresentation
                        semiDis    = min(mod((rootShift+12)-nextRepresentation,12),mod((nextRepresentation+12)-rootShift,12));
                        sumSemiDis = sum(semiDis(1:len,:))
                        [~,reRoot_idx] = min(sumSemiDis);
                        reRoot_idx
                        %end
%                         if  any(nowRoot == reRoot)
                            Ans.Chord(currentBar-1,idx(end))     = Resolution_Dim7.root{currentBar-1,idx(end)}(reRoot_idx);
%                             Ans.Chord(currentBar-1,idx(end))     = Resolution_Dim7.root{currentBar-1,idx(end)}(nowRoot == reRoot);
                            Ans.templNo(currentBar-1,idx(end))   = ceil(Ans.Chord(currentBar-1,idx(end))/12);
                            Ans.pitchNo(currentBar-1,idx(end))   = ~(ceil(mod(Ans.Chord(currentBar-1,idx(end)),12)/12))*12 + mod(Ans.Chord(currentBar-1,idx(end)),12);
                            Ans.ChordName{currentBar-1,idx(end)} = strcat(pitchName(Ans.pitchNo(currentBar-1,idx(end))), tempName(Ans.templNo(currentBar-1,idx(end))));    

%                         else
                            Untreated = [Untreated currentBar];
                            warning(['第 ' num2str(currentBar) '小節 - step3 : Fully Diminished 7th 沒有解決']);
%                         end 
                    end
                end     

            % step 3 : Diminished 7th resolution
                resolution_idx = find( Ans.sameScore(currentBar,:) == 3 );
                if resolution_idx
                    for t=1:numel(resolution_idx)        
                        nowRoot      = ~(ceil(mod(Resolution_Dim7.root{currentBar,resolution_idx(t)},12)/12))*12 + mod(Resolution_Dim7.root{currentBar,resolution_idx(t)},12);
                        rootShift    = [nowRoot, circshift(nowRoot,3), circshift(nowRoot,2), circshift(nowRoot,1)]  %^^
                        if Ans.Chord(currentBar,resolution_idx(t)+1) ~= 0 % 不是小節最後一個chord
                            %^^
                            nextRepresentation = (mod(template(Ans.templNo(currentBar,resolution_idx(t)+1),(template(Ans.templNo(currentBar,resolution_idx(t)+1),:)~=-1)) + Ans.pitchNo(currentBar,resolution_idx(t)+1)-1,12)+1)';
                            len = length(nextRepresentation);
                            if length(nextRepresentation)==3; nextRepresentation(4) = NaN; len = 3; end
                            nextRepresentation
                            semiDis    = min(mod((rootShift+12)-nextRepresentation,12),mod((nextRepresentation+12)-rootShift,12));
                            sumSemiDis = sum(semiDis(1:len,:))
                            [~,reRoot_idx] = min(sumSemiDis);
                            reRoot_idx
                            %end
                            nextRoot = ~(ceil(mod(Ans.Chord(currentBar,resolution_idx(t)+1),12)/12))*12 + mod(Ans.Chord(currentBar,resolution_idx(t)+1),12);
                            reRoot   = ~(ceil(mod((nextRoot-1),12)/12))*12 + mod((nextRoot-1),12);

%                             if  any(nowRoot == reRoot)
                                AnsChord = Resolution_Dim7.root{currentBar,resolution_idx(t)}(reRoot_idx)  %^^
%                                 AnsChord = Resolution_Dim7.root{currentBar,resolution_idx(t)}(nowRoot == reRoot);

                                Ans.Chord(currentBar,resolution_idx(t))     = AnsChord
                                Ans.templNo(currentBar,resolution_idx(t))   = ceil(AnsChord/12); 
                                Ans.pitchNo(currentBar,resolution_idx(t))   = ~(ceil(mod(AnsChord,12)/12))*12 + mod(AnsChord,12);
                                Ans.ChordName{currentBar,resolution_idx(t)} = strcat(pitchName(Ans.pitchNo(currentBar,resolution_idx(t))), tempName(Ans.templNo(currentBar,resolution_idx(t))));    

%                             else
%                                 Untreated = [Untreated currentBar];
%                                 warning(['第 ' num2str(currentBar) '小節 - step3 : Fully Diminished 7th 沒有解決']);
%                             end 

                        else
                            refNextBarRoot = 1;
                        end   

                    end

                end 


            end

        end
        Untreated

        addBar  = addBar  + NoOfBar;
        addBeat = addBeat + NoOfBeats*NoOfBar;

    end
    
    %% 存標記檔

%     fid = fopen(['chord_label/' filename '.txt'],'w');
%     accu=0; markSec = repmat(markSec,1,2);
%     for b = 1:addBar
%         for m = 1:markLen(b)
%             fprintf(fid, ['%d\t %d\t %d %s\n'], markSec(accu+m,:), b, char(Ans.ChordName{b,m}));
%         end
%         accu = accu + markLen(b);
%     end
%     fclose(fid);

%     a = authenticCadence(Ans.segPitchFirst, Ans.segOnset, Ans.Chord);
%     b = harmonicAnalysis( midiINFO, Ans.segOnset1, Ans.Chord );

