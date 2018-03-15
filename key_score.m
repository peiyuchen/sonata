% description : �p��C�Ӥp�`���C�@�ӽժ�����
% input : chord -> �M����T
%         barNote -> �p�`���Ÿ�T
% output : score -> �C�Ӥp�`�C�ӽժ�����(�M��*����)
function [score] = key_score(evaluationChord, barNote)
    keyName     = {'C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B', ...
                   'c', 'c#', 'd', 'd#', 'e', 'f', 'f#', 'g', 'g#', 'a', 'a#', 'b'};
    evaluationChord(1, :) = [];
    keyChord = make_key_chord;
    keyRatio = zeros(cell2mat(evaluationChord(end, 1)), 24);
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
        
        score.keyScaleRatio(barI, :) = key_scale_ratio(barNote{barI}); % �զ���
        score.finalScore = score.chordScoreNormalize(barI, :).* score.keyScaleRatio(barI, :);
        
        % ��p�`��key�Ĥ@�G�W
        scoreSort = sort(unique(score.finalScore), 'descend');
        if length(scoreSort) > 1 && scoreSort(1) ~= 0
            chordNameNo1 = keyName(score.finalScore == scoreSort(1))';
            score.chordNo1(1:length(chordNameNo1), barI) = chordNameNo1;
        end
        if length(scoreSort) > 2 && scoreSort(2) ~= 0
            chordNameNo2 = keyName(score.finalScore == scoreSort(2))';
            score.chordNo2(1:length(chordNameNo2), barI) = chordNameNo2;
        end  
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
% description: �p��C�Ӥp�`�Ҧ����j�p�ղզ������
% input     : note -> midi data
% output    : keyRatio -> �Ҧ��ղզ������
function [keyRatio] = key_scale_ratio(note)
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
    % new
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