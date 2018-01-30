function [segment_frame] = convertSegment_seconds_to_frames(seg_seconds,featureRate)

start_second = seg_seconds(:,1);
end_second = seg_seconds(:,2);

segment_frame = [start_second*featureRate+1, end_second*featureRate+1];

end
