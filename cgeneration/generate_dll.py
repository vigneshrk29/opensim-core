import casadi as ca
import sys
import shutil
import os 

def test(dim):
    cg = ca.CodeGenerator('foo_jac')
    arg = ca.SX.sym('arg', dim)
    from foo import foo
    y,_,_ = foo(arg)
    F = ca.Function('F',[arg],[y])
    cg.add(F)
    cg.add(F.jacobian())
    cg.generate()

if __name__ == "__main__":
    
# This is a PowerShell script that performs the following steps:
# 1. Copy foo.m, generated when running OpenSim/External_Functions/<filename>.cpp, in the current folder.
# 2. Generate c code (foo_jac.c), containing the function F and its Jacobian, in a format that can be exploited by CasADi.
# 3. Build and install a dynamically linked library (.dll) from the c code.
# 4. Delete foo.m and foo_jac.c.

# There are a few settings you may want to adjust depending on your target external function and on how things are organized on your machine:
# 1. filename is the name of your .cpp file: OpenSim/External_Functions/<filename>.cpp. 
# 2. nInputs_F is the total number of elements in the input vectors of the function F described in OpenSim/External_Functions/<filename>.cpp. 
#       In the example external functions, we have:    
#           PredSim: 87
#           PredSim_pp: 87
#           PredSim_SSCM: 87
#           PredSim_SSCM_pp: 87
#           PredSim_v2: 93
#           PredSim_v2_pp: 93
#           PredSimProsthesis: 87
#           PredSimProsthesis_pp: 87
#           TrackSim_1: 58
#           TrackSim_2: 141
#           PredSim_2D: 30
#           PredSim_2D_pp: 30
#           video_s1_A_tread: 69
#           video_s1_B_tread: 69
#           video_s1_C_tread: 69
#           video_s1_A_tread_KA: 75
#           video_s1_C_tread_KA: 75
#           video_s1_tread_pp: 69
#           video_s1_tread_KA_pp: 75
#           video_s2_C_tread_KA: 75
#           video_s2_tread_KA_pp: 75
#           video_s2_C_KA: 75
#           video_s2_KA_pp: 75
# 			PredSim_mtpPin_cm7: 93
# 			PredSim_mtpPin_pp_cm7: 93
# 3. path_opensim_build_external_functions is the path to the External_Functions folder located in the OpenSim build folder.
# 4. path_external_functions is the path to an existing folder, eg. external_functions, in which a folder <filename> will be created with the build and install folders of each external function. 
# 5. path_cgeneration is the path the folder containing this PowerShell script.
# 6. compiler is by default set to Visual Studio 2017. For Visual Studio 2015, change to "Visual Studio 14 2015 Win64".

    filename = "testLai"
    computername = os.environ.get('COMPUTERNAME', None)
    if computername == "GBW-L-W2003": # Antoine's laptop
        path_opensim_build_external_functions = "C:/Users/u0101727/Documents/Visual Studio 2017/Projects/opensim-recorder-github-py/core/build/OpenSim/External_Functions"
        path_external_functions = "C:/Users/u0101727/Documents/Visual Studio 2017/Projects/opensim-recorder-github-py/external_functions"
        path_cgeneration = "C:/Users/u0101727/Documents/MyRepositories/opensim-recorder-github-py/cgeneration"
    elif computername == "DESKTOP-RV5S4TL": # Antoine's desktop (Stanford)
        path_opensim_build_external_functions = "C:/Users/antoi/Documents/VS2017/opensim-ad/core/build/OpenSim/External_Functions"
        path_external_functions = "C:/Users/antoi/Documents/VS2017/opensim-ad/external_functions"
        path_cgeneration = "C:/Users/antoi/Documents/MyRepositories/opensim-ad/cgeneration"
    compiler = "Visual Studio 15 2017 Win64"
    
    # The remaining setting should not be changed except if you want things to be organized differently.
    fooName = "foo.py"
    
    path_external_filename = os.path.join(path_opensim_build_external_functions, filename)
    path_external_filename_foo = os.path.join(path_external_filename, fooName)
    path_external_functions_filename = os.path.join(path_external_functions, filename)
    path_external_functions_filename_build = os.path.join(path_external_functions_filename, "build")
    path_external_functions_filename_install = os.path.join(path_external_functions_filename, "install")
    shutil.copy2(path_external_filename_foo, path_cgeneration)
    
    test(int(sys.argv[1]))
    
    if not os.path.exists(path_external_functions_filename_build):
        os.makedirs(path_external_functions_filename_build) 
    if not os.path.exists(path_external_functions_filename_install):
        os.makedirs(path_external_functions_filename_install)
        
    os.chdir(path_external_functions_filename_build)
    
    cmd1 = 'cmake "' + path_cgeneration + '" -G "' + compiler + '" -DTARGET_NAME:STRING="' + filename + '"'
    os.system(cmd1)
    cmd2 = "cmake --build . --config RelWithDebInfo --target install"
    os.system(cmd2)
    
    os.remove(os.path.join(path_cgeneration, "foo_jac.c"))
    os.remove(os.path.join(path_cgeneration, fooName))       
