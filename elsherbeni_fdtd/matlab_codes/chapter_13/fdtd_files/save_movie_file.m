if ~exist('animation','var')
    return;
end

% save animations to movie files
for mi=1:size(animation,2)
    if animation(mi).save_movie && animation(mi).enable
        movie_name = [project_name '\animation' num2str(mi) '.avi'];
        disp(['saving movie file '  movie_name]);       

        writerObj = VideoWriter(movie_name);
        writerObj.FrameRate = 30;
        open(writerObj);
        for mj = 1:animation(mi).frame_number
            writeVideo(writerObj,animation(mi).frames(mj));
        end
        close(writerObj);
        disp([movie_name ' is created']);
    end
end
