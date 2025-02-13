# #!/bin/bash

# # Find highest supported python version installed
# cmds=('python3.12' 'python3.11' 'python3.10' 'python3.9')

# for cmd in "${cmds[@]}"; do
#     if $cmd --version ; then
#         python=$cmd
#         echo "Using python version" $python
#         break
#     fi
# done

# # Create virtual environment
# $python -m venv PASSAGEvenv
# source "./PASSAGEvenv/bin/activate"
# echo "Virtual environment created at" $VIRTUAL_ENV

# # Create git directory
# gitdir=$VIRTUAL_ENV/GITHUB
# mkdir $gitdir && cd $gitdir
# # echo "$PWD"

# # Pre-install fixed dependencies
# dependencies=("numpy==1.26.4" "scipy==1.14.1" "astropy==6.0.1" "pandas==2.2.3")
# for d in "${dependencies[@]}"; do
#     $python -m pip install $d
# done

# # Clone line-finding repo and install editable package
# git clone https://github.com/jwstwfss/NIRCam-line-finding.git
# cd line-finding
# $python -m pip install -e .

# # Create directory
# codedir=$VIRTUAL_ENV/CODE
# mkdir $codedir && cd $codedir
# cp $gitdir/line-finding/README.md . && cp $gitdir/line-finding/mainPOPPIES.py .
# echo "\n\nIn directory $PWD. To run the line finding code, edit the paths in "\
# "./mainPASSAGE.py, then run 'python ./mainPASSAGE.py'"