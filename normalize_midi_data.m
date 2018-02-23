% description: 將midi的duration(unit:beat)正規化成跟譜一樣
% input     : midiData           -> 原始   的midi data
% output    : normalizedMidiData -> 正規化後的midi data
function [midiData] = normalize_midi_data (midiData)
    %                  1   2+4   2    4+8    4+16     4    8+16      8      16    
    noteValue       = [1  0.75 0.5  0.375  0.3125  0.25  0.1875  0.125  0.0625];
    repNoteValue    = repmat(noteValue, size(midiData, 1), 1);

    tmp             = midiData(:, 2) - floor(midiData(:, 2)) - repNoteValue;
    tmp(tmp > 0)    = -10;
    [~,idx]         = max(tmp, [], 2);
    
    midiData(:,2) = noteValue(idx)' + floor(midiData(:,2));
    
    % onset 也要改變
    noteValue       = [noteValue 0];
    repNoteValue    = repmat(noteValue, size(midiData,1), 1);

    tmp             = midiData(:, 1) - floor(midiData(:, 1)) - repNoteValue;
    tmp(tmp > 0)    = -10;
    [~,idx]         = max(tmp, [], 2);
    
    midiData(:,1) = noteValue(idx)' + floor(midiData(:,1));
end