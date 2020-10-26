import sys
import os

# https://github.com/kkroening/ffmpeg-python
import ffmpeg

output_folder = sys.argv[1]
movie_fn = sys.argv[2] #"movie.mp4"
movie_fps = int(sys.argv[3]) #int(2)

fn_pattern = os.path.join(output_folder, "output", "frame-*.png")
output_movie = os.path.join(output_folder, "output", movie_fn)
ffmpeg.input(fn_pattern, pattern_type="glob", framerate=movie_fps).output(output_movie).run()