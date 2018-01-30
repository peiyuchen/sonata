addpath('toolbox/matlab-midi-master/src');

index = [1,2,3,4,5,8,16,17,20];
for song_no = 1:length(index)
    clearvars -except song_no index
    
    filename    = ['b_' num2str(index(song_no)) '_1'];
    
    midi_info   = readmidi(['../midi/' filename '.mid']);
    note_info   = midiInfo(midi_info, 0); % 0 不要印outputFormat資訊
    [audio, Fs] = midi2audio(note_info);
    
    audiowrite(['../midi2audio/' filename '.wav'], audio, Fs);
    
end