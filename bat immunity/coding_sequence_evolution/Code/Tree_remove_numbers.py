import os
import re
## Remove any numerical characters appearing in tree files
directory = "/tzachi_storage/yanl/bat_genes/Masked_MSA/Phy"

count = 0
for filename in os.listdir(directory):
    if filename.endswith("_tree.txt"):
        file_path = os.path.join(directory, filename)
        
        with open(file_path, "r") as file:
            contents = file.read()
        
        new_contents = re.sub(r'[\d:.]+', '', contents)  
        
        edited_file_path = os.path.join(directory, filename.replace(".txt", "_edited.txt"))
        with open(edited_file_path, "w") as file:
            file.write(new_contents)
        
        print(f"Processed file: {filename}")
        
        count += 1
        if count == 100:
            break
