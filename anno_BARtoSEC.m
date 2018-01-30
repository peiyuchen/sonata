clc;
clear;
fclose all;

midi_fp = '../midi/';

dir = '../data_annotation/pei_anno/';
index = [1,2,3,4,5,8,16,17,20];
for i=1:length(index)
    
    name = ['b_' num2str(index(i)) '_1'];
    midiInfo = readmidi_java([midi_fp '/' name '.mid'], 1);
    filename = [name '_muse_s.txt'];


    fidin=fopen([dir filename]); % ���}test2.txt���
    fidout=fopen([dir name '_s.txt'],'w'); % �Ы�MKMATLAB.txt���
    
    j = 1; k = 1;
    
    while ~feof(fidin) % �P�_�O�_����󥽧�

        tline(j,1) = {fgetl(fidin)};% �q���Ū��
        if length(tline{j,1})>=9  && strcmp(tline{j,1}(1:9), 'structure')
            data_ori(k,1) = tline(j,1);
            data_split = []; data{k,5} = [];
            data_split = textscan(tline{j,1},'%s ');
            data(k,1:length(data_split{:})) = data_split{:}';
            time(k,:) = data_split{:}(2:3)';
            for t = 1:2
                tmp = []; 
                tmp = strsplit(time{k,t}, ':');
                tmp_sec = strsplit(tmp{3}, '.');
                TIME(k,t) = str2num(tmp{1})*60 + str2num(tmp{2}) + ( str2num(tmp_sec{1}) + str2num(tmp_sec{2})/80 )/24;            
            end
             if TIME(k,2) > max(midiInfo(:,6)+midiInfo(:,7))
                TIME(k,t) = max(midiInfo(:,6)+midiInfo(:,7));
             end
            
            data_merge(k,:) = [ data(k,1) TIME(k,1) TIME(k,2) data(k,4)];
            fprintf(fidout,'%s\t%s\t%s\t%s\n',data{k,1}, num2str(TIME(k,1)), num2str(TIME(k,2)), data{k,4});
            
            k = k + 1;
        end
        j = j + 1;
        
    % if double(tline(1))>=48&&double(tline(1))<=57 % �P�_���r�ŬO�_�O�ƭ�
     % �p�G�O�Ʀr��A�⦹��ƾڼg�J���MKMATLAB.txt
    % continue % �p�G�O�D�Ʀr�~��U�@���`��
    end
    
    % fclose(fidout);
    % MK=importdata('MKMATLAB.txt'); % �N�ͦ���MKMATLAB.txt���ɤJ�u�@�Ŷ��A�ܶq�W��MK�A��ڤW������ܥX��

end

