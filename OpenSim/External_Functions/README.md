This folder contains the external functions. 

To add a new external function:
- Add a sub-folder with your `.cpp` file and a `CMakeLists.txt`. Look at the examples for inspiration.
- Make sure you set the `TARGET_NAME` in the `CMakeLists.txt` of your new folder to the name of your `.cpp` file.
- Add a new subdirectory to the `CMakeLists.txt` within `OpenSim/External_Functions`.

In visual studio, your new external function will appear under `External_Functions`.
