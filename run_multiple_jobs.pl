#!/usr/bin/perl


#SET THIS ONLY
#the input file names
#(note: not the full path name and without the file extension .nc)
@input_file_list = ("b30.042a.cam2.h0.2000-04",
                    "b30.042a.cam2.h0.2050-04",
                    "b30.042a.cam2.h0.2010-04",
                    "b30.042a.cam2.h0.2020-04",
       	           );
          

#--------------DO NOT TOUCH ANYTHING BELOW THIS LINE------------------------
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------


#this is to allow the batch script to obtain the svn version number
#Note: svn_info.txt will continually get clobbered which is ok
`svn info > svn_info.txt`;


print "\n\n";


#MAKE A RUN DIRECTORY FOR EACH INPUT FILE AND THE JOB FROM IN THERE---------

for($i=0;$i<($#input_file_list + 1);$i++){


    print "RUN NUMBER ",$i+1,"\n";
    print "SETTING UP A RUN FOR INPUT FILE $input_file_list[$i]\n";


    #CREATE NAMES
    #the name of the batch job
    $batch_job_name =  substr $input_file_list[$i], index($input_file_list[$i], 'cam');
    chomp $batch_job_name;
    #the name of the run directory
    $rundir = "rundir_".$batch_job_name;

    
    print "MAKING THE TEMP DIRECTORY FOR INPUT FILE $input_file_list[$i]\n";
    #make the temporary directory
    `mkdir $rundir`;


    #cd into the temporary directory
    print "MOVING TO THE RUN DIRECTORY $rundir\n";
    chdir $rundir;


    #we want to pass a variable to the batch script
    #so we dump the variable to a file then the batch script will source the file
    open (batch_file, '>>batch_envir.txt');
    print batch_file "setenv CLARREO_INPUT_FILE $input_file_list[$i]\n";
    close (batch_file); 


    #link the temp directory to the build directory
    `ln -v --symbolic ../* .`;


    #submit the batch job
    print "SUBMITTING BATCH JOB\n\n";
    `qsub -N $batch_job_name submit-pleiades-multiple`;


    #move back up into the build directory
    chdir "..";


}


print "FINISHED SUBMITTING ALL RUNS";
print "\n\n";

exit;
