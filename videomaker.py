print("Animator.py loaded...")
import cv2
import os
import time as tm
print("> Imports complete")
frames_folder = "./Results/"

video_name = "video.mp4"
fps = 100 # Adjust fps as needed (e.g., 24, 30)
#print("> Video is 10x slower than real time")
N = 5 * fps # Number of frames in directory 
print(f"> Properties: FPS = {fps} and total number of frames = {N}")

# Get all frame filenames
frame_filenames = [os.path.join(frames_folder, f"Frame_{i}.png") for i in range(N - 1)]  # 2001 for 0-2000 frames

# Check if any frames exist
if not frame_filenames:
    print("> Error: No frames found in the folder!")
    exit()
if N != len(frame_filenames) + 1:
    print(f"> [Error]: The number of frames specified does not match the number of elements in the selected folder({len(frame_filenames)}). Expect errors.")
    tm.sleep(2)
else:
    print("> Frames in folder have been counted and identified.")

# Read the first frame to determine video properties:
first_frame = cv2.imread(frame_filenames[0])
height, width, channels = first_frame.shape
print("> Video properties pulled from Frame 0")

# Define the video writer:
fourcc = cv2.VideoWriter_fourcc(*'XVID') # Used for .mp4
video_writer = cv2.VideoWriter(video_name, fourcc, fps, (width, height))
print("> Defined video writer")

# Process and write frames:
for frame_filename in frame_filenames:
    frame = cv2.imread(frame_filename)
    # You can add pre-processing here (e.g., resizing, color correction)
    video_writer.write(frame)
print("> Frames Processed")

# Release resources:
video_writer.release()
cv2.destroyAllWindows()

print(f"> Video '{video_name}' has been created. Errors (if any) will have been announced by the terminal.")