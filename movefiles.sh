#!/bin/sh

find ./ -name "b30.042a*" -fprintf moveuptwodir.sh 'mv "%h/%f" "%h/../../../../output_files"\n'
