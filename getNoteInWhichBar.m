function [ midiInfo ] = getNoteInWhichBar( filename, track )

% testing code
%     clear all; close all; clc;
%     filename = 'b_16_1.mid';
%     track = 1;
    [ midiInfo, timeSig ] = midi_Preprocess(filename, track);    % »Ý­n midiinfo 1(onset beat) 4(pitch)

    
% initial
    addBar          =   0;
    addBeat         =   0;
    

    for i = 1:size(timeSig,1)
        % set
        NoOfBeats       =   timeSig(i,1)*(4/(2^timeSig(i,2)));
        
        if i ~= size(timeSig,1)
            NoOfBar     =   ceil( (timeSig(i+1,5) - timeSig(i,5) ) /NoOfBeats);
        else
            NoOfBar     =   ceil( ( midiInfo(end,1) + midiInfo(end,2) - timeSig(i,5) ) / NoOfBeats);
        end
        
        for j = 1:NoOfBar
            
            note_X = []; note_Y = [];
            bar_start   = (addBeat+(j-1)*NoOfBeats);
            bar_end     = addBeat+j*NoOfBeats;
            notesOfBar  = intersect(find(midiInfo(:,1)>=bar_start), find(midiInfo(:,1)<bar_end));
            
            midiInfo(notesOfBar, 11) = addBar + j;
            
        end

        addBar  = addBar + NoOfBar;
        addBeat = addBeat + NoOfBeats*NoOfBar;

    end
     
    
end
    