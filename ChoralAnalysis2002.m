 % implement : Algorithms for Choral Analysis 2002
 
    clear; close all; clc;

    tempName = {'maj','7','min','Fully dim7','Half dim7','dim3','X'};
    pitchName = {'C :','C#:','D :','D#:','E :','F :','F#:','G :','G#:','A :','A#:','B :'};

    filepath = '../midi/';
    filename = 'mz_545_1_noRepeat';%beethoven_string_quartet_135, mz_16_1 O48N2AOM quartet_16_3 dim7_1
    [midiINFO, timeSig]  = midi_Preprocess([filepath filename]);

    midiINFO = normalize_midi_data(midiINFO);

    addBeat = 0; addBar = 0; refNextBarRoot = 0;

    eva_i = 1;
    evaChord = {'小節','拍數(onset)','調性','和弦名稱','和弦編號','備註'};% new

    isSave = 1;
    
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
            midiINFO(notesOfBar, 12) = currentBar;

            % 前個小節圓滑線至此小節的音符
            innerOffset     = intersect(find(sum(midiINFO(:,1:2),2)<barEnd) , find(sum(midiINFO(:,1:2),2)>barStart));
            overOnset       = intersect(find(midiINFO(:,1)<barStart), innerOffset);
            overNote        = midiINFO(overOnset, :);
            overNote(:,2)   = overNote(:,2) + barStart - overNote(:,1);
            overNote(:,1)   = barStart;

            noteBar         = midiINFO(notesOfBar, :);   % 此小節內的音符

            noteBar         = [noteBar; overNote];          % 把onset改成此小節一開始


            if ~isempty(noteBar)

                noteBar       = sortrows(noteBar,1);      % 依據第一行排序 小到大
                noteBar       = triDetection(noteBar);    % tri 處理

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

                template            = [ 0,4,7,-1;  0,4,7,10;   0,3,7,-1;  
                                        0,3,6, 9;  0,3,6,10;   0,3,6,-1 ];

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
                MARK                        = HarmAn( weight, partitionPoint, 1:partitionPoint );
                Ans.Mark(j,1:length(MARK))  = MARK;

                % 紀錄 最佳路徑的onset
                segOnset = miniSeg(:,1)' - barStart;
                Ans.segOnset(j,1:length(MARK)-1) = segOnset(MARK(1:end-1));

                Dim7Time                    = 1;

    %% 最佳路徑的和弦標記

                for k = 2:length(MARK)
                    Ans.Chord(j,k-1)       = I(MARK(k-1),MARK(k));
                    Ans.sameScore(j,k-1)   = -1;
                    Ans.Score(j,1)          = weight(MARK(k-1),MARK(k));


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
            % step 3 : Diminished 7th resolution
                                if ceil(HighScoreChord/12) == 4
                                    Resolution_Dim7.idx(j,Dim7Time)  = k-1;
                                    Resolution_Dim7.root{j,k-1}     = HighScoreChord;
                                    step     = 3
                                    Dim7Time = Dim7Time + 1;
                                    AnsChord = Ans.Chord(j,k-1);
                                else
                                    AnsChord = 84%97; %Ans.Chord(j,k-1);
                                    step     = 4
                                    bar      = j
%                                     error('???');
                                end
                            end
                        end

                        Ans.Chord(j,k-1) = AnsChord;
                        Ans.sameScore(j,k-1) = step;
                    end

                    Ans.templ_no(j,k-1)  = ceil(Ans.Chord(j,k-1)/12);
                    Ans.pitch_no(j,k-1)  = ~(ceil(mod(Ans.Chord(j,k-1),12)/12))*12 + mod(Ans.Chord(j,k-1),12);
                    Ans.ChordName{j,k-1} = strcat(pitchName(Ans.pitch_no(j,k-1)), tempName(Ans.templ_no(j,k-1)));

                    if refNextBarRoot
                        idx               = find(Resolution_Dim7.idx(j-1)~=0);
                        nowRoot           = ~(ceil(mod(Resolution_Dim7.root{j-1,idx(end)},12)/12))*12 + mod(Resolution_Dim7.root{j-1,idx(end)},12);
                        nextRoot          = Ans.pitch_no(j,k-1);
                        reRoot            = ~(ceil(mod((nextRoot-1),12)/12))*12 + mod((nextRoot-1),12);
                        refNextBarRoot    = ~refNextBarRoot;

                        if  any(nowRoot == reRoot)
                            Ans.Chord(j-1,idx(end))     = Resolution_Dim7.root{j-1,idx(end)}(nowRoot == reRoot);
                            Ans.templ_no(j-1,idx(end))  = ceil(Ans.Chord(j-1,idx(end))/12);
                            Ans.pitch_no(j-1,idx(end))  = ~(ceil(mod(Ans.Chord(j-1,idx(end)),12)/12))*12 + mod(Ans.Chord(j-1,idx(end)),12);
                            Ans.ChordName{j-1,idx(end)} = strcat(pitchName(Ans.pitch_no(j-1,idx(end))), tempName(Ans.templ_no(j-1,idx(end))));    
                            
                            evaChord(eva_i,4) = Ans.ChordName{j-1,idx(end)};
                            evaChord{eva_i,5} = Ans.Chord(j-1,idx(end));
                        else
                            error('step3 : Fully Diminished 7th 沒有解決');
                        end 
                    end
                    
                    eva_i = eva_i + 1;
                    evaChord{eva_i,1} = j;
                    evaChord{eva_i,2} = segOnset(MARK(k-1));
                    evaChord(eva_i,4) = Ans.ChordName{j,k-1};
                    evaChord{eva_i,5} = Ans.Chord(j,k-1);
                end

            % step 3 : Diminished 7th resolution
                resolution_idx = find( Ans.sameScore(j,:) == 3 );
                if resolution_idx
                    for t=1:numel(resolution_idx)        
                        nowRoot      = ~(ceil(mod(Resolution_Dim7.root{j,resolution_idx(t)},12)/12))*12 + mod(Resolution_Dim7.root{j,resolution_idx(t)},12);

                        if Ans.Chord(j,resolution_idx(t)+1) ~= 0
                            nextRoot = ~(ceil(mod(Ans.Chord(j,resolution_idx(t)+1),12)/12))*12 + mod(Ans.Chord(j,resolution_idx(t)+1),12);
                            reRoot   = ~(ceil(mod((nextRoot-1),12)/12))*12 + mod((nextRoot-1),12);

                            if  any(nowRoot == reRoot)
                                AnsChord = Resolution_Dim7.root{j,resolution_idx(t)}(nowRoot == reRoot);

                                Ans.Chord(j,resolution_idx(t))     = AnsChord;
                                Ans.templ_no(j,resolution_idx(t))  = ceil(AnsChord/12); 
                                Ans.pitch_no(j,resolution_idx(t))  = ~(ceil(mod(AnsChord,12)/12))*12 + mod(AnsChord,12);
                                Ans.ChordName{j,resolution_idx(t)} = strcat(pitchName(Ans.pitch_no(j,resolution_idx(t))), tempName(Ans.templ_no(j,resolution_idx(t))));    
                                
                                nowBarIdx = find(cell2mat(evaChord(2:end,1))==j);
                                evaChord(nowBarIdx(resolution_idx(t))+1,4) = Ans.ChordName{j,resolution_idx(t)};
                                evaChord{nowBarIdx(resolution_idx(t))+1,5} = Ans.Chord(j,resolution_idx(t));
                            else
                                error('step3 : Fully Diminished 7th 沒有解決');
                            end 

                        else
                            refNextBarRoot = 1;
                        end   
                    end
                end 
            end
        end

        addBar  = addBar  + NoOfBar;
        addBeat = addBeat + NoOfBeats*NoOfBar;

    end
    
    if isSave
        cell2csv(['chordEva/eva_' filename '.csv'], evaChord)
    end


