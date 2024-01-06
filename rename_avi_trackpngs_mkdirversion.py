import glob,shutil,os

avi_files = glob.glob("*Box*/*avi")
track_files = glob.glob("*Box*/*png")

for avi_file in avi_files:
	filename = os.path.basename(avi_file)
	dest_dir = avi_file.split("/")[0] + "/" + avi_file.split("/")[0] + "_avifiles"
	if not os.path.exists(dest_dir):		
		os.mkdir(dest_dir)
	shutil.move(avi_file, dest_dir)
for track_file in track_files:
	filename = os.path.basename(track_file)
	dest_dir = "./trackingpngfiles/" + track_file.split("/")[0] + "_trackedlines"
	if not os.path.exists(dest_dir):		
		os.mkdir(dest_dir)
	#print(dest_dir + "/" + track_file.split("/")[1])
	if not os.path.exists(dest_dir + "/" + track_file.split("/")[1]):
		shutil.move(track_file, dest_dir)
	#print("moving " + track_file + " to " + dest_dir)
