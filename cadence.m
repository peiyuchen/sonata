 clear all; close all; clc;
 
 addBar          = 0;
 addBeat         = 0;
 
 filename = ['../midi/b_20_1'];
 [midiALL, timeSig]  = midi_Preprocess(filename);
 
 for i = 1:size(timeSig,1)
     
        NoOfBeats       =   timeSig(i,1)*(4/(2^timeSig(i,2)));
        
        if i ~= size(timeSig,1)
            NoOfBar     =   max(ceil( ((timeSig(i+1,5)-1) - timeSig(i,5) ) / NoOfBeats), 1);
        else
            NoOfBar     =   max(ceil( ( (midiALL(end,1)+midiALL(end,2)) - timeSig(i,5) ) / NoOfBeats), 1);
        end
        
        currentBeat = addBeat;
        
        for j = 1:NoOfBar
            
            barStart   = addBeat + (j-1) * NoOfBeats;
            barEnd     = addBeat +  j    * NoOfBeats;
            notesOfBar = intersect(find(midiALL(:,1)>=barStart), find(midiALL(:,1)<barEnd));
            currentBar = addBar  +  j;
            currentBeat= currentBeat +  ~(j==1)*NoOfBeats;
            
            if ~isempty(notesOfBar) % 此小節內有音符                
                Note(currentBar).X          = midiALL(notesOfBar, 1) - currentBeat;
                Note(currentBar).pitch      = midiALL(notesOfBar, 4);
                Note(currentBar).pitchHist  = hist(mod(midiALL(notesOfBar, 4), 12),0:11);

            end
            
        end
        
        addBar  = addBar  + NoOfBar;
        addBeat = addBeat + NoOfBeats*NoOfBar;
        
    end
 
 
 
 
