import os

# Get all the filenames in the current directory
filenames = [f for f in os.listdir() if f.startswith("openmm")]

# Loop through the filenames
for filename in filenames:
    # Split the filename by '_' to get the first, second, and third parameters
    drug = filename.split('_')[7]
    excipient = filename.split('_')[8]
    traj = filename
    assert os.path.exists(f"{excipient}.mol2"), f"Error! Missing {excipient}.mol2 file!"
    assert os.path.exists(f"{drug}.mol2"), f"Error! Missing {drug}.mol2 file!"
    # Make the call to make_Analysis.py
    command = f"python make_Analysis.py -d {drug} -e {excipient} --traj {traj}"
    os.system(command)

print("Job's Done")

