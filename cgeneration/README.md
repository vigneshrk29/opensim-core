You can run the `generate_dll.ps1` Powershell script (right click - Run with Powershell) to obtain the dynamically linked library `<filename>.dll` corresponding to the external function F described in `OpenSim/External_Functions/<filename>.cpp`.
This step requires MATLAB and CasADi. It also assumes `matlab.exe` is in your user path.

`<filename>.dll` will be saved in the `<path_external_functions>/<filename>/install/bin` folder, where `<path_external_functions>` can be specified in the `generate_dll.ps1` Powershell script.
You can then import `<filename>.dll` in MATLAB or Python as [a CasADi`s external function](https://web.casadi.org/docs/#casadi-s-external-function) when formulating your optimal control problem.
This script assumes certain settings (number of inputs to F and paths) that you may want to edit depending on your target external function and on how things are organized on your machine. To do so, open the Powershell script with Notepad++ and follow the guidelines.
