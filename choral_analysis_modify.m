%%%% motify : "Algorithms for Choral Analysis 2002"
%%%% date   : 18/01/22
%%%% content: �M�����R�t��k�A�H�p�`�����i����R�C
%%%% input  : filename
%%%% output : chords in each bar

% description: �ק� "Algorithms for Choral Analysis 2002":�M�����R�t��k�A�H�p�`�����i����R
% date   : 18/01/22
% input     : filename -> midi�ɦW(�����⭵�ɩ�bmidi��Ƨ���)
%             parameter-> �@�Ǳ��󪺰Ѽ�
%                 isDownbeatAndDurWeight : �H���I�������I�A���Ū��׬��v��
%                 isLowWeight : ����̧C���v���դj(�⭿)
%                 isDownbeatWeight : ���I�W�����v���դj(�⭿)
%             isSave-> �M�����G�O�_�s��
% output    : evaluationChord -> �M�����G
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
    %                 �j�T        �ݤC         �p�T       ��
    template    = [ 0,4,7,-1; 0,4,7,10; 0,3,7,-1; 0,3,6,-1];  
    evaluationIdx = 1;
    evaluationChord = {'�p�`','���(onset)','�թ�','�M���W��','�M���s��','�Ƶ�'};% new




    for j = 1:length(barNote)
        noteBar = barNote{j};
        if ~isempty(noteBar)
            % ���̤p���q�����I�����ǭ�
            [pitchClass, minSeg] = min_seg_pitchclass(noteBar, barOnset(j));   
            if parameter.isDownbeatAndDurWeight
                [pitchClass, lowestPClass, lowestP, minSeg] = downbeat_pitchclass(noteBar, barOnset(j), parameter.isDownbeatWeight);
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
                     for kk = 1:templateNum
                         for p = 1:12
                             templateNo = mod(template(kk,(template(kk,:)~=-1)) + p - 1 , 12);
                             Template = zeros(1,12); 
                             Template(1,templateNo + 1) = 1;

                             rootS(ii, jj, (kk - 1) * 12 + p) = sum(sum(pitchClass(ii:jj-1, templateNo(1) + 1)));                 % �ڭ���score
                             P(ii, jj, (kk - 1) * 12 + p)     = sum(sum(pitchClass(ii:jj-1, templateNo + 1)));                  % ���q ��, �M�� ��
                             N(ii, jj, (kk - 1) * 12 + p)     = sum(sum(pitchClass(ii:jj-1, setdiff(1:12,templateNo+1))));    % ���q ��, �M���S��  
                             M(ii, jj, (kk - 1) * 12 + p)     = sum((sum(pitchClass(ii:jj-1, :),1) == 0) .* Template);        % ���q�S��, �M�� ��
                             S(ii, jj, (kk - 1) * 12 + p)     = ...
                                 P(ii, jj, (kk - 1) * 12 + p) - (N(ii, jj, (kk - 1) * 12 + p) + M(ii, jj, (kk-1) * 12 + p));
                             %%%%% �[�̧C��weight %%%%%%
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

            %% ���q�G�̨θ��|�t��k HarmAn
            optimalPath = HarmAn(maxScoreMatrix, partitionNum);

            % ���� �̨θ��|��onset
            segOnset = [minSeg(:,1)' - barOnset(j), minSeg(end,2) - barOnset(j)];

            Ans.optimalPath(j, 1:length(optimalPath)) = optimalPath;
            Ans.segOnset(j, 1:length(optimalPath)) = segOnset(optimalPath);

            %% �̨θ��|���M���аO
            for k = 2:length(optimalPath)
                Ans.chord{j,k-1}       = scoreMatrixIdx(optimalPath(k-1),optimalPath(k));
                Ans.tieBreaking(j,k-1) = {''};
                Ans.score(j,k-1)       = maxScoreMatrix(optimalPath(k-1),optimalPath(k));

                %% ���Ƥ@�˪���
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
% description: �H���I��partition point�A�����C��minSeg���������v���A�H��minSeg���̧C��
% input:  noteBar      -> �p�`����
%         barOnset     -> �p�`�}�l��beat
% output: pitchClass,  -> �C��minSeg���Ҧ����A�Bweight�����Ū��סA�j�p��N*12 (N��minSeg�ƶq,12���@�ӤK�׭����ƶq)
%         lowestPitchClass -> �C��minSeg���̧C���A�Bweight�����Ū��סA�j�p��N*12 (N��minSeg�ƶq,12���@�ӤK�׭����ƶq)
%         lowestP      -> �C��minSeg���̧C���A�B���H�@�ӤK�רӬ�(��mod12)
%         minSeg       -> �C��minSeg��onset offset��beat
% �䥪��̧C���v�����j�A���H���륪�k��A�ҥH����q�̧C���B�P�ɦ���ӭ�
function [pitchClass, lowestPitchClass, lowestP, minSeg] = downbeat_pitchclass(noteBar, barOnset, isDownbeatWeight)
    beatInBar = noteBar(:,9);
    minSeg = [0:beatInBar-1; 1:beatInBar]'; % �H���I���̤p���q (�� 0 1 2 3 ��)
    pitchClass = zeros(size(minSeg, 1), 12);
    lowestPitchClass = zeros(size(minSeg, 1), 12);
    lowestP = inf * ones(size(minSeg, 1), 1);
    
    for s = 1:size(minSeg, 1)
        unit = (1/2)^4;
        noteOn = zeros(size(noteBar,1),1);
        noteOff = zeros(size(noteBar,1),1);
        isOnsetInt = zeros(size(noteBar, 1), 1);
        noteNum = zeros(ceil((minSeg(s, 2) - minSeg(s, 1)) / unit), 1); % ����minSeg���Cunit���סA���X�ӭ�
        noteLowestIdx = zeros(ceil((minSeg(s, 2) - minSeg(s, 1)) / unit), 1); % �C��unit�̧C���O���ӭ�
        noteLowestPitch = inf * ones(ceil((minSeg(s, 2) - minSeg(s, 1)) / unit), 1); % �C��unit�̧C��������

        for n = 1:size(noteBar, 1)
            noteOn(n) = noteBar(n, 1) - barOnset;
            noteOff(n) = sum(noteBar(n,1:2)) - barOnset;
            isOnsetInt(n) = (noteOn(n) == fix(noteOn(n))); % isDownbeatWeight
            
            if noteOn(n) < minSeg(s, 1); noteOn(n) = minSeg(s, 1); end
            if noteOff(n) > minSeg(s, 2); noteOff(n) = minSeg(s, 2); end
            
            if  noteOff(n) - noteOn(n) > 0   % �P�_���bminSeg������
                pitchNo = mod(noteBar(n, 4), 12) + 1;
                pitchClass(s, pitchNo) = pitchClass(s, pitchNo) + noteOff(n) - noteOn(n);
                
                if isDownbeatWeight && isOnsetInt(n)
                    pitchClass(s, pitchNo) = pitchClass(s, pitchNo) + noteOff(n) - noteOn(n);
                end

                onUnitNum  = floor((noteOn(n) - minSeg(s,1)) / unit) + 1; % ��on�Moff�O��minSeg��������unit
                offUnitNum = ceil((noteOff(n) - minSeg(s,1)) / unit);
                repPitch = ones(offUnitNum - onUnitNum + 1, 1) * noteBar(n, 4);

                [minPitch, idx] = min([noteLowestPitch(onUnitNum:offUnitNum) repPitch], [], 2);
                noteNum(onUnitNum:offUnitNum) = noteNum(onUnitNum:offUnitNum) + 1;

                if ~isempty(find(idx == 2, 1))
                    noteLowestIdx(find(idx == 2) + onUnitNum - 1) = n;
                    noteLowestPitch(onUnitNum:offUnitNum) = minPitch;
                end
            end
            
            %% ��ӭ��P�ɡA��C4�H�U����
            if n == size(noteBar, 1)
                if ~isempty(find(noteNum == 2, 2)) %%  ��ӭ��P�ɨ��̧C��,c4�H�U
                    lowestP(s) = min(noteLowestPitch(noteNum == 2));
                    if lowestP(s) < 60
                        idx = find(noteLowestPitch == lowestP(s)); %�n��
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
                if ~isempty(find(noteNum > 2, 2)) %%  �T�ӭ��P�ɨ��̧C���A�u��
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
            %%  ��ӭ��P�ɨ��̧C��
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


% description: �Ҧ����q��pitchClass
% input     : noteBar -> �p�`���� midi data
%             barOnset-> �p�`�}�l��beat
% output    : pitchClass -> �����p�`�Ҧ����q��pitch class
%             minSeg -> �p�`�Ҧ��̤p���qonset offset
function [pitchClass, minSeg] = min_seg_pitchclass(noteBar, barOnset)

    % downbeat ���̤p���q
%     beatInBar = noteBar(:, 9);
% %     minSeg = [0:beatInBar-1; 1:beatInBar]';
%     minSeg = unique(round([noteBar(:, 1); sum(noteBar(:, 1:2), 2); (0:beatInBar-1)' + barOnset], 4));
%     minSeg(:, 2) = [minSeg(2:end,1); beatInBar(1) + barOnset];

    % �̤p���q��onset offset
    minSeg = unique(round([noteBar(:, 1); sum(noteBar(:, 1:2), 2)], 4));
    minSeg(:, 2) = [minSeg(2:end, 1); round(sum(noteBar(end, 1:2)), 4)];
    
    % ��minSeg�Ӥp���h��
    miniTime = 0.025;
    if minSeg(end, 1) == minSeg(end, 2); minSeg(end, :)=[]; end
    minSeg(diff(minSeg, 1, 2) <= miniTime, :) = [];

    % �����C�ӳ̤p���q�����ǭ�
    pitchClass = zeros(length(minSeg), 12);
    for k=1:size(noteBar, 1)
        s = find(abs(round(noteBar(k, 1), 4) - minSeg(:, 1)) <= miniTime); % �o�ӭ��q����seg�}�l
        e = find(round(sum(noteBar(k, 1:2)), 4) >= minSeg(:, 2)); % �o�ӭ��q����seg����
        pitchNum = mod(noteBar(k, 4), 12);
        pitchClass(s:e(end), pitchNum + 1) = pitchClass(s:e(end), pitchNum + 1) + 1;
    end
end


function [maxScoreMatrix, I, S] = record_min_seg_score(pitchClass)
    %                            �j�T        �ݤC         �p�T       ��
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
                    rootS(ii,jj,(kk-1)*12+p) = sum(sum(pitchClass(ii:jj-1,templateNo(1)+1))); % �ڭ���score
                    s(ii,jj,(kk-1)*12+p)     = sum(sum(pitchClass(ii:jj-1,templateNo+1  )));                  % ���q ��, �M�� ��
                    N(ii,jj,(kk-1)*12+p)     = sum(sum(pitchClass(ii:jj-1,:)))-s(ii,jj,(kk-1)*12+p);    % ���q ��, �M���S��            
                    M(ii,jj,(kk-1)*12+p)     = sum((sum(pitchClass(ii:jj-1,:),1)==0).*Template);        % ���q�S��, �M�� ��
                    S(ii,jj,(kk-1)*12+p)     = s(ii,jj,(kk-1)*12+p) - ( N(ii,jj,(kk-1)*12+p) + M(ii,jj,(kk-1)*12+p) );
                end
            end
         end
     end

    [maxScoreMatrix,I] = max(S,[],3);
end


