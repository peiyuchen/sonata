% % Major scale consists of notes
% M_C         = {'C' , 'D' , 'E' , 'F' , 'G' , 'A' , 'B' };
% M_Csharp    = {'C#', 'D#', 'F' , 'F#', 'G#', 'A#', 'C' };
% M_D         = {'D' , 'E' , 'F#', 'G' , 'A' , 'B' , 'C#'};
% M_Dsharp    = {'D#', 'F' , 'G' , 'G#', 'A#', 'C' , 'D' };
% M_E         = {'E' , 'F#', 'G#', 'A' , 'B' , 'C#', 'D#'};
% M_F         = {'F' , 'G' , 'A' , 'A#', 'C' , 'D' , 'E' };
% M_Fsharp    = {'F#', 'G#', 'A#', 'B' , 'C#', 'D#', 'F'};
% M_G         = {'G' , 'A' , 'B' , 'C' , 'D' , 'E' , 'F#'};
% M_Gsharp    = {'G#', 'A#', 'C' , 'C#', 'D#', 'F' , 'G' };
% M_A         = {'A' , 'B' , 'C#', 'D' , 'E' , 'F#', 'G#'};
% M_Asharp    = {'A#', 'C' , 'D' , 'D#', 'F' , 'G' , 'A' };
% M_B         = {'B' , 'C#', 'D#', 'E' , 'F#', 'G#', 'A#'};
% 
% 
% % NATURE Minor scale consists of notes
% Nm_C        = {'C' , 'D' , 'D#', 'F' , 'G' , 'G#', 'A#'};
% Nm_Csharp   = {'C#', 'D#', 'E' , 'F#', 'G#', 'A' , 'B' };
% Nm_D        = {'D' , 'E' , 'F' , 'G' , 'A' , 'A#', 'C' };
% Nm_Dsharp   = {'D#', 'F', 'F#', 'G#', 'A#', 'B' , 'C#'};
% Nm_E        = {'E' , 'F#', 'G' , 'A' , 'B' , 'C' , 'D' };
% Nm_F        = {'F' , 'G' , 'G#', 'A#', 'C' , 'C#', 'D#'};
% Nm_Fsharp   = {'F#', 'G#', 'A' , 'B' , 'C#', 'D' , 'E' };
% Nm_G        = {'G' , 'A' , 'A#', 'C' , 'D' , 'D#', 'F' };
% Nm_Gsharp   = {'G#', 'A#', 'B' , 'C#', 'D#', 'E' , 'F#'};
% Nm_A        = {'A' , 'B' , 'C' , 'D' , 'E' , 'F' , 'G' };
% Nm_Asharp   = {'A#', 'C' , 'C#', 'D#', 'F' , 'F#', 'G#'};
% Nm_B        = {'B' , 'C#', 'D' , 'E' , 'F#', 'G' , 'A' };
% 
% % HARMONIC Minor scale consists of notes
% Hm_C        = {'C' , 'D' , 'D#', 'F' , 'G' , 'G#', 'B' };
% Hm_Csharp   = {'C#', 'D#', 'E' , 'F#', 'G#', 'A' , 'C' };
% Hm_D        = {'D' , 'E' , 'F' , 'G' , 'A' , 'A#', 'C#'};
% Hm_Dsharp   = {'D#', 'F', 'F#', 'G#', 'A#', 'B' , 'D'  };
% Hm_E        = {'E' , 'F#', 'G' , 'A' , 'B' , 'C' , 'D#'};
% Hm_F        = {'F' , 'G' , 'G#', 'A#', 'C' , 'C#', 'E' };
% Hm_Fsharp   = {'F#', 'G#', 'A' , 'B' , 'C#', 'D' , 'F' };
% Hm_G        = {'G' , 'A' , 'A#', 'C' , 'D' , 'D#', 'F#'};
% Hm_Gsharp   = {'G#', 'A#', 'B' , 'C#', 'D#', 'E' , 'G' };
% Hm_A        = {'A' , 'B' , 'C' , 'D' , 'E' , 'F' , 'G#'};
% Hm_Asharp   = {'A#', 'C' , 'C#', 'D#', 'F' , 'F#', 'A' };
% Hm_B        = {'B' , 'C#', 'D' , 'E' , 'F#', 'G' , 'A#'};

pitchName   = {'C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B'};
pitchNumber = [  0,    1,   2,    3,   4,   5,    6,   7,    8,   9,   10,  11];

% KeyN = { M_C{:};  M_Csharp{:};  M_D{:};  M_Dsharp{:};  M_E{:};  M_F{:};  M_Fsharp{:};  M_G{:};  M_Gsharp{:};  M_A{:};  M_Asharp{:};  M_B{:}; ...
%         Nm_C{:}; Nm_Csharp{:}; Nm_D{:}; Nm_Dsharp{:}; Nm_E{:}; Nm_F{:}; Nm_Fsharp{:}; Nm_G{:}; Nm_Gsharp{:}; Nm_A{:}; Nm_Asharp{:}; Nm_B{:}};
% 
% KeyH = { M_C{:};  M_Csharp{:};  M_D{:};  M_Dsharp{:};  M_E{:};  M_F{:};  M_Fsharp{:};  M_G{:};  M_Gsharp{:};  M_A{:};  M_Asharp{:};  M_B{:}; ...
%         Hm_C{:}; Hm_Csharp{:}; Hm_D{:}; Hm_Dsharp{:}; Hm_E{:}; Hm_F{:}; Hm_Fsharp{:}; Hm_G{:}; Hm_Gsharp{:}; Hm_A{:}; Hm_Asharp{:}; Hm_B{:}};
% cat(1,cat(1,keyN(:).consistNote))
% clearvars -except pitchName pitchNumber KeyN KeyH

% init
% KEY             = struct();
% keyNumber       = zeros(1,7);
% diatonicTriad   = -1*ones(4,7);

% keyTmp          = {KeyN, KeyH, KeyH};
load('Key_DiatonicTriad_V7.mat')
KEYN        = cat(1,keyN(:).consistNote);
KEYH        = cat(1,keyH(:).consistNote);
KEYF        = cat(1,keyFusion(:).consistNote);

%% 
xlsFile = 'chord_k309_m1';
[number, text, rawData] = xlsread(['chordGT/' xlsFile '.xlsx']);
scale   = {'I', 'II', 'III','IV', 'V', 'VI', 'VII','V7','viio','iio'};
degree  = [  1,    2,     3,   4,   5,    6,     7,   5,    7 ,    2, ]; 
tempName = {'maj','7','min','dim','xxx','X'};
newChord = rawData; newChord(1,5) = {'和弦編號'};

for i=2:length(rawData)
    
    % key
    pNo = pitchNumber(strcmpi(pitchName, rawData{i,3}))+1;
    
    if rawData{i,3}>='a'&& rawData{i,3}<='z'
        pNo = pNo + 12;
    end

    % degree
    scaleDegree = regexp(rawData{i,4}, '/', 'split');
    
    seven       = cell2mat(strfind(scaleDegree(1), '7'));
    dim         = cell2mat(strfind(scaleDegree(1), 'o'));
    scaleDegree = strrep(scaleDegree, {'o'}, {''});   % 把減和弦刪掉
    scaleDegree = strrep(scaleDegree, {'7'}, {''});   % 把7和弦刪掉
    
    %% 和弦根音
    
    Pde = degree(strcmpi(scale, scaleDegree(end)));   % 級數
    special = '';
    
    % 拿坡里和弦＆德國六和弦
    if strcmpi('N', scaleDegree(end))
        Pde = 2;
        special = 'N';
    elseif strcmpi('Ger', scaleDegree(end))
        Pde = 6;
        special = 'Ger';
    end
    
    
%     newChord{i,4} = KeyH(pNo, Pde);
    newChord{i,4} = KEYF(pNo, Pde);
    
    de = find(strcmpi(scale, scaleDegree(end)));
    if ~isempty(de) && de == 7 && pNo>12         % 如果是VII代表用自然小調音階
        newChord{i,4} = KEYN(pNo, Pde);
    end
    newChord{i,5} = find(strcmpi(newChord{i,4}, pitchName))-1;
    
    % 拿坡里和弦＆德國六和弦
    if ~isempty(special)
        newChord{i,5} = mod(newChord{i,5} - 1, 12);
        newChord{i,4} = pitchName(pitchNumber==newChord{i,5});
    end
    
    % 副 屬/導 和弦
    if length(scaleDegree) == 2
        tmp = degree(strcmpi(scale, scaleDegree(1)));
        if tmp == 5; Pde2 = 7;  end
        if tmp == 7; Pde2 = 11; end
        newChord{i,5} = mod(newChord{i,5} + Pde2, 12);
        newChord{i,4} = pitchName(pitchNumber==newChord{i,5});
    end
    tempName = {'maj','7','min','dim','xxx','X'};
    
    %% 和弦樣版

    add = [];
    if scaleDegree{1}(1)>='a' && scaleDegree{1}(1)<='z'
        temp = 3;
    else
        temp = 1;
    end
    
    if dim; temp = 4; end
    if seven
        if strcmpi(scaleDegree{1},'V')
            temp = 2;
        else
            add = '7';
        end
    end
    newChord{i,4} = [newChord{i,4}{:} ':' tempName{temp} add];
    newChord{i,5} = newChord{i,5} + 1 + 12 * (temp-1);
    
    if ~isempty(special)
        newChord{i,6} = special;
    end
    if length(scaleDegree) == 2
        newChord{i,6} = 'second';
    end
 
end
cell2csv(['chordGT/trans_' xlsFile '.csv'],newChord)



timeSig = 4;
barNum = newChord{end,1};
unit = 0.5;

GT_tmpNum   = zeros(barNum, timeSig/unit);
GT_tmpName  = repmat({'-'}, barNum, timeSig/unit);%cell (barNum, timeSig/unit);
for i=2:length(newChord)
    on = newChord{i,2}/unit+1; 
    if i~=length(newChord)
        off = newChord{i+1,2}/unit;
        if newChord{i+1,2}==0; off = 4/unit; end
    else
        off = 4/unit;
    end
    GT_tmpNum (newChord{i,1}, on:off) = newChord{i,5};
    GT_tmpName{newChord{i,1}, on} = newChord{i,4};
end
