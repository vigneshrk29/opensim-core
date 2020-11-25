# Number of inputs for the different cases. Please adjust $NInputs_F accordingly.
# PredSim: 87
# PredSim_pp: 87
# PredSim_SSCM: 87
# PredSim_SSCM_pp: 87
# PredSim_v2: 93
# PredSim_v2_pp: 93
# PredSimProsthesis: 87
# PredSimProsthesis_pp: 87
# TrackSim_1: 58
# TrackSim_2: 141
# video_subject1_A_tread: 69
# video_subject1_B_tread: 69 
# video_subject1_C_tread: 69
# video_subject1_tread_pp: 69
# video_s1_A_tread_KA: 75
# video_s1_C_tread_KA: 75
# video_s1_tread_KA_pp: 75

# User settings
$filename = "video_s1_C_tread_KA"
$NInputs_F = 75
$path_core_external = "C:\Users\u0101727\Documents\Visual Studio 2017\Projects\opensim-recorder-github-contact\core\build\OpenSim\External_Functions"
$path_external_build = "C:\Users\u0101727\Documents\Visual Studio 2017\Projects\opensim-recorder-github-contact\externalFunctions"
$path_cgeneration = "C:\Users\u0101727\Documents\MyRepositories\opensim-recorder-github-contact\cgeneration"

# Fixed settings (except if you want things to be organized differently)
$foo = "foo.m"
$path_external_filename = Join-Path $path_core_external $filename
$path_external_filename_foo = Join-Path $path_external_filename $foo
$path_external_build_filename = Join-Path $path_external_build $filename
$path_external_build_filename_build = Join-Path $path_external_build_filename "build"
$path_external_build_filename_install = Join-Path $path_external_build_filename "install"

Move-Item -Path $path_external_filename_foo -Destination $foo

matlab -nosplash -nodesktop -r "generate_foo_jac($NInputs_F)" quit
while (!(Test-Path "foo_jac.c")) { Start-Sleep 1 }

mkdir $path_external_build_filename
mkdir $path_external_build_filename_build
mkdir $path_external_build_filename_install

cd $path_external_build_filename_build
cmake $path_cgeneration -G "Visual Studio 15 2017 Win64" -DTARGET_NAME:STRING=$filename 
cmake --build . --config RelWithDebInfo --target install

cd $path_cgeneration
del "foo_jac.c"
del $foo
    