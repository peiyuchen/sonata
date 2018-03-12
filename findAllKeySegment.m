% clear student				% 清除 student 變數
% 	student.name = '洪鵬翔';			% 加入 name 欄位
% 	student.id = 'mr871912';			% 加入 id 欄位
% 	student.scores = [58, 75, 62];			% 加入 scores 欄位
% 	% 以下是新加入的第二筆資料
% 	student(2).name = '邱中人';
% 	student(2).id = 'mr872510';
% 	student(2).scores = [25, 36, 92];
clc
keyName     = {'C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B', ...
               'c', 'c#', 'd', 'd#', 'e', 'f', 'f#', 'g', 'g#', 'a', 'a#', 'b'};

data = finalScore';
evaChord = evaluationChord;

keyNum = size(data, 1);
barNum = size(data, 2);

% case 1 : find 0 to 1 -> start
zero = find(data == 0);
one = find(data ~= 0);
zero(zero > keyNum * (barNum - 1)) = []; % 最後一行的0刪掉
segOnset = intersect(zero + keyNum, one); % 把0後移一行後跟1比較，位置一樣的地方為0->1
% case 2 : bar = 1, have 1
segOnset = sort([segOnset; one(one <= keyNum)]);

% 
bar = cell2mat(evaChord(2:end,1));
keyChord = make_key_chord;
clear segment

% 紀錄每個片段
for i = 1:length(segOnset)
    start = ceil(segOnset(i) / keyNum);
    
    segment(i).keyIdx = segOnset(i) - (start - 1) * keyNum;
    segment(i).key = keyName{segment(i).keyIdx};
    segment(i).start = ceil(segOnset(i) / keyNum);
    segment(i).end = segment(i).start;
    segment(i).lens = 1;
    nowI = segment(i).keyIdx;
    nowJ = segment(i).start;
    segment(i).score = data(nowI, nowJ);
    segment(i).chordI_IV_V = 0;
    segment(i).chordI_V = 0;
    segment(i).chordI = 0;
    segment(i).chordV = 0;
    chordNum = 0;
    
    evaChordBarCol = cell2mat(evaChord(2:end, 1));
    barIdx = find(evaChordBarCol == nowJ) + 1;
    barChord = cell2mat(evaChord(barIdx, 5));
    barChordLens = cell2mat(evaChord(barIdx, 2));
    chordNum = chordNum + 4;%+ length(barChord);
    for j = 1:length(barChord)
        if j == length(barChord)
            lens = 4 - barChordLens(j);
        else
            lens = barChordLens(j+1) - barChordLens(j);
        end
        
        segment(i).chordI_IV_V = segment(i).chordI_IV_V + numel(intersect(keyChord(segment(i).keyIdx, 1:4), barChord(j)))*lens;
        segment(i).chordI_V = segment(i).chordI_V + numel(intersect(keyChord(segment(i).keyIdx, [1 3 4]), barChord(j)))*lens;
        segment(i).chordI = segment(i).chordI + numel(intersect(keyChord(segment(i).keyIdx, 1), barChord(j)))*lens;
        segment(i).chordV = segment(i).chordV + numel(intersect(keyChord(segment(i).keyIdx, [3 4]), barChord(j)))*lens;
    end
    
    while data(nowI, nowJ) && nowJ ~= barNum
        nowJ = nowJ + 1;
        if data(nowI, nowJ)
            segment(i).lens = segment(i).lens + 1;
            segment(i).score = segment(i).score + data(nowI, nowJ);
            segment(i).end = nowJ;
            
            evaChordBarCol = cell2mat(evaChord(2:end, 1));
            barIdx = find(evaChordBarCol == nowJ) + 1;
            barChord = cell2mat(evaChord(barIdx, 5));
            barChordLens = cell2mat(evaChord(barIdx, 2));
%             chordNum = chordNum + length(barChord);
            for j = 1:length(barChord)
                if j == length(barChord)
                    lens = 4 - barChordLens(j);
                else
                    lens = barChordLens(j+1) - barChordLens(j);
                end
                segment(i).chordI_IV_V = segment(i).chordI_IV_V + numel(intersect(keyChord(segment(i).keyIdx, 1:4), barChord(j)))*lens;
                segment(i).chordI_V = segment(i).chordI_V + numel(intersect(keyChord(segment(i).keyIdx, [1 3 4]), barChord(j)))*lens;
                segment(i).chordI = segment(i).chordI + numel(intersect(keyChord(segment(i).keyIdx, 1), barChord(j)))*lens;
                segment(i).chordV = segment(i).chordV + numel(intersect(keyChord(segment(i).keyIdx, [3 4]), barChord(j)))*lens;
            end
        end
    end
    segment(i).score = segment(i).score / segment(i).lens;
    chordNum = (segment(i).end - segment(i).start + 1)*4;
    segment(i).chordI_IV_V = segment(i).chordI_IV_V / chordNum;
    segment(i).chordI_V = segment(i).chordI_V / chordNum;
    segment(i).chordI = segment(i).chordI / chordNum;
    segment(i).chordV = segment(i).chordV / chordNum;
end

% 依照長度排序
[~, idx] = sort(cat(1, segment.lens), 'descend');
segment = segment(idx);

% 找該小節有的segment取最長
lensMore2 = find(cat(1, segment.lens) > 1);
scoreMore8 = find(cat(1, segment.score) > 0.8);
segIdx = intersect(lensMore2, scoreMore8);

candiKey = cell(length(segIdx), barNum);
candiKeyIdx = zeros(length(segIdx), barNum);
candiKeyScore = zeros(length(segIdx), barNum);

for i = 1:length(segIdx)
    for j = segment(segIdx(i)).start:segment(segIdx(i)).end
        candiKey{i, j} = segment(segIdx(i)).key;
    end
    candiKeyIdx(i, segment(segIdx(i)).start:segment(segIdx(i)).end) = 1;
    candiKeyScore(i, segment(segIdx(i)).start:segment(segIdx(i)).end) = segment(segIdx(i)).score;

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