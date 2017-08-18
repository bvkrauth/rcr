!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RCR program
!
! Description:  Performs calculations for RCR model
!
! Author:       Brian Krauth
!               Department of Economics
!               Simon Fraser University
!
! Usage:
!
!       RCR.exe [infile outfile logfile]
!
!       where
!
!           infile      An optional argument giving the name
!                       of the input file.  Default is IN.TXT.
!
!           outfile     An optional argument giving the name of
!                       the output file. Default is OUT.TXT.
!
!           logfile     An optional argument giving the name of
!                       the output file. Default is LOG.TXT.
!
!       The RCR program will read in the INFILE, perform
!       the calculations, and then write the results to OUTFILE.
!       The program may also report information on its status to
!       LOGFILE.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program RCR
    ! The RCRLIB module contains system-specific code
    use rcrlib_gnu, only : SP, DP
    ! The RCRUTIL module contains everything else
    use rcrutil_gnu, only : get_command_arguments, read_data, estimate_model, write_results, &
        infile, outfile, logfile, detail_file, moment_vector, lambda_range, result_matrix
    implicit none

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Begin run code
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Load in arguments from call to program
    call get_command_arguments(infile,outfile,logfile,detail_file)

    ! Read in the data from INFILE (side effect: allocation/creation of moment_vector, lambda_range, result_matrix)
    call read_data(trim(infile))

    ! Perform the calculations and put the results in result_matrix (side effect: allocation of theta_segments, writing to detail_file)
    call estimate_model(moment_vector,lambda_range,result_matrix)

    ! Write out the data to OUTFILE
    call write_results(result_matrix,trim(outfile))

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! End run code
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end program RCR

