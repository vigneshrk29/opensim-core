# Number of inputs for the different cases. Please adjust $NInputs_F accordingly.
# ShoulderModel_consB: 42

# User settings
$filename = "ShoulderModel_consB"
$NInputs_F = 42
$path_core_external = "C:\Users\u0101727\Documents\Visual Studio 2017\Projects\opensim-recorder-github-scapula\core-build\OpenSim\External_Functions"
$path_external_build = "C:\Users\u0101727\Documents\Visual Studio 2017\Projects\opensim-recorder-github-scapula\External_Functions"
$path_cgeneration = "C:\Users\u0101727\Documents\MyRepositories\opensim-recorder-github-scapula\cgeneration"

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
    