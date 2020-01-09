import glob , os

# Obtain list of VTK files in directory
vtk_list = glob.glob("*.vtk")

# Rename the extensions to .txt files to process
for vtk in vtk_list: 
    # Access the last four characters in the path and replace them with .txt
    os.rename(vtk, vtk[:-4] + ".txt")

txt_list = glob.glob("*output.txt")
print(txt_list)

for text in txt_list: 
    with open(text, "w") as foo:
        line = foo.readlines()
        # print (line)
        for number in line:
            # print (number)
            if "UNSTRUCTURED_GRID" in number:
                number = number.replace("UNSTRUCTURED_GRID", "POLYDATA")
            
            if ""