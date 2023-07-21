run config_test `0'
which rcr
di "Data from file ${fname}" _newline ///
   "Parameter values for tests: os = ${os}, exe = ${exe}, tol = ${tol}"
*******************************************************************
* Test SAVE option (undocumented)
*******************************************************************
* This option tells stata to create files called in.txt, out.txt, and log.txt
* First delete those files if they already exist
capture erase in.txt
capture erase out.txt
capture erase log.txt
* Second, run the command with the option included
rcr SAT Small_Class ${controls}, save
* Finally, check for existence of the files
confirm file in.txt
confirm file out.txt
confirm file log.txt
