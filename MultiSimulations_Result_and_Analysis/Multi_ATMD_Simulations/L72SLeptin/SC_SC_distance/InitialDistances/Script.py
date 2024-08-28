import os

# Define a function to copy content from one file to another
def copy_content(source_file, target_file):
	with open(source_file, "r") as source:
		content = source.read()
		with open(target_file, "w") as target:
			target.write(content)


# Get a list of all files in the working directory
files = os.listdir()

# Filter the files that start with "file_"
text_files = [file for file in files if file.startswith("file_")]

# Iterate over each text file
for text_file in text_files:
	# Extract the integer part from the file name
	file_number = int(text_file.split("_")[2].split(".")[0])
	# Calculate the new file name
	new_file_name = f"File_72_{file_number + 21}.txt"
	
	# Copy content from the original file to the new file
	copy_content(text_file, new_file_name)
