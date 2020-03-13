module rcrutil
    use rcrlib, only : SP, DP, sysv, is_finite, sort, stderr, initialize
    implicit none
    public get_command_arguments, read_data, estimate_model, write_results, &
        infile, outfile, logfile, moment_vector, lambda_range, result_matrix
    !---------------------------------------------------------------------------
    ! Declarations for global variables
    !---------------------------------------------------------------------------
    ! Standard_buffer_length is the typical length for a string buffer
    integer, parameter :: standard_buffer_length=2048
    ! Filename_length is the maximum length for a filename
    integer, parameter :: filename_length=1024
    ! Internal_big_number and internal_small_number
    real(kind=DP), parameter :: internal_big_number=huge(1.0_dp), internal_small_number=1.0e-10_dp
    ! External_big_number
    real(kind=DP) :: external_big_number, internal_infinity, internal_nan
    ! It will also be handy to have a string buffer
    character(len=standard_buffer_length) :: buffer
    ! Moment_vector will be the input vector of moments.  See READ_DATA for more information.
    real(kind=DP), dimension(:), allocatable :: moment_vector
    ! Theta_segments will be a vector describing the lambda(theta) function.  See ESTIMATE_THETA_SEGMENTS.
    real(kind=DP), dimension(:), allocatable :: theta_segments
    ! Lambda_range will be the input matrix of (lambda_L,lambda_H) pairs.  See READ_DATA.
    real(kind=DP), dimension(:,:), allocatable :: lambda_range
    ! Result_matrix will be the output matrix, calculated in ESTIMATE_MODEL and output in WRITE_RESULTS
    real(kind=DP), dimension(:,:), allocatable :: result_matrix
    ! Infile will be the input filename
    ! Outfile will be the output filename
    ! Logfile will be the log filename
    character(len=filename_length) :: outfile="out.txt", infile="in.txt", logfile="log.txt", detail_file=""

contains

    !---------------------------------------------------------------------------
    ! SEQ function
    !
    ! Description: Produces a real sequence of length n, between two endpoints.
    !
    ! Usage: seq(first,last,n)
    !
    !           first   A real number
    !           last    A real number
    !           n       An integer
    !
    ! Effect: Returns an evenly-spaced real vector of length N whose first element
    !         is FIRST and whose last element is LAST.
    !
    !---------------------------------------------------------------------------
function seq(first,last,n)
    real(kind=DP), intent(in) :: first,last
    integer, intent(in) :: n
    real(kind=DP), dimension(n) :: seq
    integer :: i
    do i=1,n
        seq(i) = first + (last-first)*real(i-1,kind=DP)/real(n-1,kind=DP)
    end do
end function seq


    !---------------------------------------------------------------------------
    ! BKSOLVE function
    !
    ! Description: Solves a linear system of equations
    !
    ! Usage: bksolve(XX,XY)
    !
    !           XX      A real matrix
    !           XY      A real matrix
    !
    ! Effect: Returns the real matrix corresponding to inv(XX)*XY
    !
    !---------------------------------------------------------------------------
function bksolve(XX,XY)
    real(kind=DP), dimension(:,:), intent(in) :: XX
    real(kind=DP), dimension(size(XX,2),1), intent(in) :: XY
    real(kind=DP), dimension(size(XX,1),size(XX,2)) :: a
    real(kind=DP), dimension(size(XX,2),1) :: bksolve
    a = XX
    bksolve = XY
    call sysv(a,bksolve)
end function bksolve

    function thetastar(moment_vector)
        real(kind=DP), dimension(:), intent(in) :: moment_vector
        real(kind=DP) :: thetastar
        real(kind=DP), dimension(6) :: sm
        sm=simplify_moments(moment_vector)
        ! thetastar is defined as
        ! sm(6) = cov(yhat,zhat)
        ! sm(5) = var(zhat)
        !
        ! The check_moments subroutine should ensure that
        ! var(zhat) >= 0 and that if var(zhat)=0 -> cov(yhat,zhat)=0.
        !
        ! Special values: If var(zhat)=0, then thetastar = 0/0 = NaN.
        !
        thetastar = sm(6)/sm(5)
    end function thetastar

    !---------------------------------------------------------------------------
    ! LAMBDASTAR function
    ! THETASTAR function
    ! LAMBDA0 function
    !
    ! Description: Estimates the scalar model parameter thetastar
    !
    ! Usage: lambdastar(moment_vector)
    !        thetastar(moment_vector)
    !        lambda0(moment_vector)
    !
    ! Arguments:
    !
    !       moment_vector       The vector of cross-moments (input)
    !
    ! Effect: The function returns a scalar equal to the estimated parameter.
    !
    !---------------------------------------------------------------------------
    function lambdastar(moment_vector)
        real(kind=DP), dimension(:), intent(in) :: moment_vector
        real(kind=DP) :: lambdastar
        real(kind=DP), dimension(6) :: sm
        sm=simplify_moments(moment_vector)
        !
        ! lambdastar is defined as sqrt( var(z)/var(zhat) - 1)
        ! sm(2) = var(z)
        ! sm(5) = var(zhat)
        !
        ! The check_moments subroutine should ensure that
        ! var(z) > 0 and that var(z) >= var(zhat) >= 0.
        ! This implies that lambdastar >= 0.
        !
        ! Special values: If var(zhat) = 0, then lambdastar = +Infinity
        !
        lambdastar = sqrt(max(sm(2)/sm(5),1.0_dp) - 1.0_dp)
    end function lambdastar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINES CALLED BY THE MAIN PROGRAM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !---------------------------------------------------------------------------
    ! GET_COMMAND_ARGUMENTS subroutine
    !
    ! Description: Obtains the command arguments, if any, provided when this
    !              program was called, and uses them to set the various
    !              input/output files.
    !              This subroutine uses the Fortran 2003 intrinsic subroutine
    !              GET_COMMAND_ARGUMENT and may thus be system-specific.
    !              Other non-Fortran-2003-compliant compilers may have
    !              a library or intrinsic function called GETARG that does
    !              the same thing.  If there is no such function, this subroutine
    !              can simply be omitted from the program.  In that case, the
    !              filenames will just be the defaults.
    !
    ! Usage: get_command_arguments(infile,outfile,logfile)
    !
    !           infile      The name of the input file.
    !           outfile     The name of the output file
    !           logfile     The name of the log file
    !
    ! Effect: Updates the values of infile, outfile, and logfile to reflect the
    !         user-provided command arguments.
    !---------------------------------------------------------------------------
    subroutine get_command_arguments(infile,outfile,logfile,detail_file)
        character(len=*), intent(inout) :: infile,outfile,logfile,detail_file
        character(len=standard_buffer_length) :: buffer
        ! Get the first argument and store it in BUFFER
        call get_command_argument(1,buffer)
        ! If it is nonempty, then...
        if (len_trim(buffer) > 0) then
            ! Put the result in INFILE
            infile = buffer
            ! Get the second argument
            call get_command_argument(2,buffer)
            ! If it is nonempty, then...
            if (len_trim(buffer) > 0) then
                ! Put the result in OUTFILE
                outfile = buffer
                ! Now get the third argument
                call get_command_argument(3,buffer)
                ! If it is nonempty, then...
                if (len_trim(buffer) > 0) then
                    ! Put the result in LOGFILE
                    logfile=buffer
                    ! Now get the fourth argument
                    call get_command_argument(4,buffer)
                    ! If it is nonempty, then...
                    if (len_trim(buffer) > 0) then
                        ! Put the result in LOGFILE
                        detail_file=buffer
                    end if
                end if
            end if
        end if
    end subroutine get_command_arguments

    !---------------------------------------------------------------------------
    ! READ_DATA subroutine
    !
    ! Description: Reads in the data from a file
    !
    ! Usage: read_data(infile)
    !
    !           infile      The name of the input file to read
    !
    ! Effect: Sets the value of the global real number EXTERNAL_BIG_NUMBER and
    !         creates the following global real arrays:
    !
    !       moment_vector       With values read in from infile
    !       lambda_range        With values read in from infile
    !       result_matrix       With undefined values for now
    !
    !---------------------------------------------------------------------------
    subroutine read_data(infile)
        ! Infile is the name of the input file
        character(len=*), intent(in) :: infile
        ! File_found is a logical we will use to inquire about INFILE's existence
        logical :: file_found
        ! I is a counter
        integer :: i
        ! IOS is used for I/O exception handling
        integer :: ios
        ! N_MOMENTS and  N_LAMBDA will be read in the first line of INFILE, and
        ! will give us the size of MOMENT_VECTOR and LAMBDA_MATRIX respectively
        integer :: n_moments, n_lambda
        ! U will be used to assign unit numbers to files
        integer :: u
        ! K will be calculated below for validating the size of N_MOMENTS.
        integer :: k
        ! CURRENT_DATE and CURRENT_TIME are self-explanatory
        character(len=10) :: current_date, current_time
        ! Initialize floating-point variables.
        call initialize(internal_infinity,internal_nan)
        ! Get the current date and time
        call date_and_time(current_date,current_time)
        ! Start a new log file and write the current date and time to it.
        call write_to_logfile("Log file " // trim(logfile)  // " for RCR version 1.0", new=.TRUE.)
        call write_to_logfile("Run at " // current_time(1:2) // ":" // current_time(3:4) // " on " // current_date(5:6) // "/" // current_date(7:8) // "/" // current_date(1:4))
        call write_to_logfile("Reading data from input file " // trim(infile) // ".")
        ! Check to see if INFILE exists
        inquire (file=infile,exist=file_found)
        call validate(file_found,msg="IO Error: Input file " // trim(infile) // " does not exist.")
        ! Given that INFILE exists, open it
        u = get_next_unit()
        open (unit=u,file=infile,iostat=ios,form="formatted",action="read",position="rewind",status="old")
        ! Check to see if INFILE can be opened
        write (unit=buffer,fmt=*) "IO Error: Input file " , trim(infile) , " exists but cannot be opened.  IOS=" , ios , "."
        call validate(ios==0,msg=trim(buffer))
        !
        ! Line 1 of INFILE should contain three numbers:
        !
        !   N_MOMENTS: an integer giving the size of Moment_vector
        !   N_LAMBDA: an integer giving the size of Lambda_vector
        !   EXTERNAL_BIG_NUMBER: a real number giving the largest finite number in Stata
        !
        ! Read the line in
        read (unit=u,iostat=ios,fmt=*) n_moments, n_lambda, external_big_number
        ! Check to make sure the line got read in correctly
        call validate(ios==0,msg="Error in line 1 of " // infile // ": Unable to read (possibly incorrect format).")
        ! Write the result to the logfile
        write (unit=buffer,fmt=*) "Line 1: n_moments = ", n_moments, ", n_lambda = ",n_lambda, ", external_big_number = ", external_big_number , "."
        call write_to_logfile(trim(buffer(2:)))
        ! Check to make sure n_lambda is a valid (i.e., positive) value
        call validate(n_lambda > 0,"Error in line 1 of " // infile // ": n_lambda is negative or zero." )
        ! Check to make sure n_moments is a valid value
        !   1. It must be at least 9 (i.e., there must be at least one explanatory variable)
        call validate( (n_moments .ge. 9) ,"Error in line 1 of " // infile // ": n_moments < 9.")
        !   2. The number of implied explanatory variables must be an integer
        !   TODO: Add a better description here.
        k = floor((sqrt(9.0 + 8.0*real(n_moments))-1.0)/2.0)
        call validate( (2 * (n_moments + 1) ) == k**2 + k,"Error in line 1 of " // infile // ": Invalid value of n_moments.")
        ! Check to make sure external_big_number is a valid value
        call validate(external_big_number > 0.0_dp,"Error in line 1 of " // infile // ": external_big_number <= 0).")
        ! If external_big_number is bigger than internal_big_number, then issue a warning but don't stop program.
        ! TODO: I'm not satisfied with this.
        if (external_big_number > internal_big_number) then
            call write_to_logfile("Warning: largest real number in executable is less than largest in Stata")
        end if
        ! Now that we know how big they are, we can allocate the arrays
        allocate(moment_vector(n_moments),lambda_range(n_lambda,2),result_matrix(2*n_lambda+3,n_moments+1),stat=ios)
        ! Check to make sure we have allocated the arrays successfully
        call validate(ios == 0,"Error when attempting to allocate memory to variables.")
        !
        ! LINE TWO contains the moment vector.  This is a real vector of length n_moments.
        !
        read (unit=u,iostat=ios,fmt=*) moment_vector
        ! Check to make sure the line got read in correctly
        call validate(ios==0,"Error reading line 2 of " // infile // ": incorrect format.")
        ! Write the result to the logfile
        ! But don't write moment_vector to the buffer. If 10 or more explanatory variables,
        ! the buffer is too small.  That would lead to an overflow and kill the program.
        call write_to_logfile("Line 2: moment_vector = ", arr=moment_vector)
        !
        ! LINES THREE AND ABOVE contain lambda_range, a real matrix with n_lambda rows
        ! and two columns.  Note that the Stata program always has n_lambda=1, though
        ! I've left open the possibility of n_lambda > 1 for future development.
        !
        do i=1,n_lambda
            ! Read in the next line
            read (unit=u,iostat=ios,fmt=*) lambda_range(i,:)
            ! Check to make sure the line got read in correctly
            write (unit=buffer,fmt=*) (2+i)
            call validate(ios==0,"Error reading line " // trim(buffer) //  " of " // infile // ".")
            ! Write the result to the logfile
            write (unit=buffer,fmt=*) "Line ", i+2, ": lambda_range = ", lambda_range(i,:)
            call write_to_logfile(trim(buffer))
            ! Round up to infinity if necessary
            where (abs(lambda_range(i,:)) >=  external_big_number) lambda_range(i,:) = sign(internal_infinity,lambda_range(i,:))
            ! Check to make sure lambda is in the correct order.
            call validate(lambda_range(i,1) <= lambda_range(i,2),"Error in line 3 of " // infile // ": Upper bound on lambda is less than lower bound.")
            ! TODO: If lambda_range is narrow (or a single point) we run into some problems.  The code below will
            ! work around those problems but may not be the best solution.
!           if ((lambda_range(i,2) - lambda_range(i,1)) < 0.00001_dp) then
!               lambda_range(i,2) = lambda_range(i,1) + 0.00001_dp
!           end if
            ! Report changes made to lambda_range, if any.
            write (unit=buffer,fmt=*) "For calculations, lambda_range, row", i, " is set to ", lambda_range(i,:)
            call write_to_logfile(trim(buffer))
        end do
        ! That's it, close the file
        close (unit=u,iostat=ios)
        ! No need to catch an error closing this file, as far as I can think.
        ! Done!
        call write_to_logfile("Data successfully loaded from file " // infile // ".")
    end subroutine read_data


    !---------------------------------------------------------------------------
    ! VALIDATE subroutine
    !
    ! Description: Checks to make sure that a particular condition holds,
    !              and shuts the program down if it doesn't.
    !
    ! Usage: validate(q,msg)
    !
    !           q       A logical condition
    !           msg     A character string
    !
    ! Effect: If q = .TRUE., then nothing.  If q=.FALSE., then
    !         the error message MSG is written to the log file
    !         and to STDERR, and the program is stopped.
    !
    !---------------------------------------------------------------------------
    subroutine validate(q,msg)
        logical, intent(in) :: q
        character(len=*), intent(in) :: msg
        ! If q = .TRUE., do nothing
        if (.not. q) then
            ! Otherwise, write the provided error message to the logfile
            call write_to_logfile(msg)
            ! And stop the program
            call die(msg // " See log file " // trim(logfile) // " for details.")
        end if
    end subroutine validate

    !---------------------------------------------------------------------------
    ! GET_NEXT_UNIT function
    !
    ! Description: Finds the next available unit number for file I/O.
    !
    ! Usage: get_next_unit()
    !
    !
    ! Effect: Returns the smallest unit number that isn't currently attached
    !         to an open file.
    !
    !---------------------------------------------------------------------------
    function get_next_unit()
        integer :: i, get_next_unit
        logical :: opened
        ! I'm assuming here that there are fewer than 100 files open.
        do i=1,100
            ! See if unit i is open
            inquire(unit=i,opened=opened)
            ! If not, then we have the next unit
            if (.not. opened) then
                ! Set get_next_unit to i
                get_next_unit = i
                ! Get out of the loop
                exit
            end if
        end do
    end function get_next_unit



    !---------------------------------------------------------------------------
    ! WRITE_TO_LOGFILE subroutine
    !
    ! Description: Writes a string to the log file
    !
    ! Usage: write_to_logfile(str,new,arr)
    !
    ! Inputs
    !        str              a character string of any length
    !        new              an optional logical scalar
    !        arr              an optional real array
    !
    ! Effect: Writes the string STR and (if supplied) the array ARR to the log
    !         file.  If NEW=.TRUE., then the existing log file is replaced.
    !         Otherwise it is appended.
    !
    !
    !---------------------------------------------------------------------------
    subroutine write_to_logfile(str,new,arr)
        character(len=*), intent(in) :: str
        logical, intent(in), optional :: new
        real(kind=DP), dimension(:), intent(in), optional :: arr
        logical :: file_found,tmpnew
        integer :: ios,i,ialt,u
        ! The default is to append to the logfile, but if the optional argument new=.TRUE., then we replace
        if (present(new)) then
            tmpnew = new
        else
            tmpnew = .FALSE.
        end if
        i=1
        ! See if the file exists already
        inquire(file=logfile,exist=file_found)
        u = get_next_unit()
        if (file_found) then
            if (tmpnew) then
                ! If the file exists, but new=.TRUE., then overwrite it
                open(unit=u,file=logfile,iostat=ios,action="write",position="rewind",status="old")
            else
                ! If the file exists but new = .FALSE. or undefined, then append
                open(unit=u,file=logfile,iostat=ios,action="write",position="append",status="old")
            end if
        else
            ! If the file doesn't exist, then create it
            open(unit=u,file=logfile,iostat=ios,action="write",status="new")
        end if
        i=2
        ! Assuming the file could be opened correctly,...
        if (ios == 0) then
            ! Write STR (and maybe ARR) to the file
            if (present(arr)) then
                write (unit=u,iostat=ios,fmt=*) str, arr
            else
                write (unit=u,iostat=ios,fmt=*) str
            end if
            i=3
            close (unit=u,iostat=ialt)
        end if
        ! If something went wrong, that isn't a fatal error, so we don't need to shut down.
        ! But we can't report the problem in the log file (since the log file is the thing that is having problems).
        ! So we report the problem to standard error.
        if (ios /= 0) then
            if (i==1) then
                write (unit=stderr,fmt=*) "Unable to open log file for writing, ios=", ios
            else
                if (i==2) then
                    write (unit=stderr,fmt=*) "Unable to write to log file, ios=", ios
                else
                    write (unit=stderr,fmt=*) "Unable to close log file, ios=", ios
                end if
            end if
        end if
    end subroutine write_to_logfile



    !---------------------------------------------------------------------------
    ! ESTIMATE_MODEL subroutine
    !
    ! Description: Performs the actual calculations.
    !
    ! Usage: estimate_model(moment_vector,lambda_range,result_matrix)
    !
    ! Arguments:
    !
    !       moment_vector       The vector of cross-moments (input)
    !       lambda_range        The matrix containing the range of lambda (input)
    !       result_matrix       The matrix of results (output)
    !
    ! Effect: The model is estimated and the result returned in RESULT_MATRIX
    !         In addition, the vector THETA_SEGMENTS is created for internal use
    !         by a call to the subroutine ESTIMATE_THETA_SEGMENTS.
    !
    !---------------------------------------------------------------------------
    subroutine estimate_model(moment_vector,lambda_range,result_matrix)
        real(kind=DP), dimension(:), intent(in) :: moment_vector
        real(kind=DP), dimension(:,:), intent(in) :: lambda_range
        real(kind=DP), dimension(:,:), intent(out) :: result_matrix
        integer :: i,j
        call write_to_logfile("Estimating model.")
        ! Check to make sure the moments are consistent
        call check_moments(moment_vector)
        ! We have closed forms for the global parameters lambdastar, thetastar, and lambda(0),
        ! so we just estimate them directly.
        result_matrix(1,:) = estimate_parameter(lambdastar,moment_vector,"lambdastar")
        result_matrix(2,:) = estimate_parameter(thetastar,moment_vector,"thetastar")
        result_matrix(3,:) = estimate_parameter(lambda0,moment_vector,"lambda0")
        ! Here we get to the main estimation problem.  We need to find the range of theta values
        ! consistent with the lambda(theta) function falling in lambda_range.  We have a closed form
        ! solution for lambda(theta), but finding its inverse is an iterative problem.
        !
        ! STEP 1: Estimate THETA_SEGMENTS, which is a global real vector indicating all critical
        !         points (i.e., points where the derivative is zero or nonexistent) of the function
        !         lambda(theta).  The function is continuous and monotonic between these points.
        !         Note that we don't know a priori how many critical points there will be, and so we
        !         don't know how big THETA_SEGMENTS will be.
        call estimate_theta_segments(moment_vector,result_matrix(2,1))
        ! STEP 2: For each row of lambda_range (i.e., each pair of lambda values)...
        do i=1,size(lambda_range,1)
            ! j is the row in result_matrix corresponding to lambda_range(i,:)
            j = 2+2*i
            ! Estimate the corresponding theta range, and put it in result_matrix
            result_matrix(j:j+1,:) = estimate_theta(moment_vector,lambda_range(i,:),result_matrix(1,1),result_matrix(2,1))
        end do
    end subroutine estimate_model


    !---------------------------------------------------------------------------
    ! ESTIMATE_THETA function
    !
    ! Description: Estimates the range of theta consistent with a given lambda_range
    !
    ! Usage: estimate_theta(moment_vector,lambda_range,theta_segments,lambdastar,thetastar)
    !
    ! Arguments:
    !
    !       moment_vector       The vector of cross-moments (input)
    !       lambda_range        The matrix containing the range of lambda (input)
    !       lambdastar          The scalar estimate of lambdastar
    !       thetastar           The scalar estimate of thetastar
    !
    ! Effect: The model is estimated and the result is returned.
    !
    !---------------------------------------------------------------------------
    function estimate_theta(moment_vector,lambda_range,lambdastar,thetastar)
        real(kind=DP), dimension(:), intent(in) :: moment_vector, lambda_range
        real(kind=DP), intent(in) :: lambdastar, thetastar
        integer, parameter :: ntab=10, nmax=10
        real(kind=DP), parameter :: con=1.4_dp, con2=con*con, big=huge(1.0_dp), safe=2.0_dp
        real(kind=DP), dimension(2,size(moment_vector)+1) :: estimate_theta
        real(kind=DP), dimension(:) :: dtheta(1),  simplified_moments(6), current_theta_range(2), current_lambda_range(2)
        real(kind=DP), dimension(:) :: dmoments(size(moment_vector)), deps(size(moment_vector))
        real(kind=DP) :: theta, h=1.0e-1_dp, hh, err, dfridr, errmax=0.0_dp
        real(kind=DP), dimension(ntab-1) :: errt, fac
        real(kind=DP), dimension(ntab,ntab) :: a
        integer :: i,j,k,m,n,ierrmin
        ! These arrays are made allocatable so that we don't have a stack overflow if they turn out to be very big
        real(kind=DP), dimension(:), allocatable :: important_thetas, lambda_segments
        logical, dimension(:), allocatable :: inrange
        allocate(important_thetas(5*size(theta_segments)-4),lambda_segments(5*size(theta_segments)-4),inrange(5*size(theta_segments)-4))
        ! Check to make sure that lambdastar is not in lambda_range.  If so, theta is completely unidentified.
        if ( between(lambdastar,lambda_range)) then
            estimate_theta(1,1) = -internal_infinity
            estimate_theta(2,1) = internal_infinity
            estimate_theta(:,2:) = internal_nan
            return
        end if
        ! IMPORTANT_THETAS is a list of theta values for which lambda(theta) needs to be calculated.
        ! We don't know in advance how many important values there will be, so we make IMPORTANT_THETAS
        ! way too big, and initialize it to all zeros (this choice is arbitrary)
        important_thetas = 0.0_dp
        ! Get simplified moments
        simplified_moments = simplify_moments(moment_vector)
!       estimate_theta=0.0_dp
        ! k is the number of actual important theta values in IMPORTANT_THETAS
        k=1
        ! Go piece by piece through theta_segments
        do i=1,(size(theta_segments)-1)
            ! Get the next pair of thetas.  This represents a range of thetas to check
            current_theta_range = theta_segments(i:i+1)
            ! Skip ahead to the next pair if thetastar is in the current theta range
            if (.not. is_finite(thetastar) .or. (current_theta_range(1) > thetastar) .or. (current_theta_range(2) < thetastar)) then
                ! Otherwise, calculate the range of lambdas associated with that range of thetas
                current_lambda_range = lambdafast(current_theta_range,simplified_moments)
                ! For each of the values in lambda_range
                do j=1,2
                    ! See if that value satisfies lambda(theta)-lambda(j)=0 for some theta in current_theta_range
                    if (( lambda_range(j) > minval(current_lambda_range)) .and. (lambda_range(j)) < maxval(current_lambda_range)) then
                        ! If so, find the theta such that lambda(theta)-lambda(j)=0 and put it in
                        ! our list of IMPORTANT_THETAS.  Of course we can't quite find the exact theta.
                        important_thetas(k) = zbrent(lambda_minus_lambda,current_theta_range(1),current_theta_range(2),1.0e-200_dp, (/ lambda_range(j), simplified_moments /))
                        k = k+1
                    end if
                end do
            end if
        end do
        ! Add THETA_SEGMENTS to the list of IMPORTANT_THETAS
        important_thetas(k:k+size(theta_segments)-1) = theta_segments
        ! Add the OLS theta to the list of IMPORTANT_THETAS
!       important_thetas(size(important_thetas)) = simplified_moments(3)/simplified_moments(2)
        ! Calculate lambda(theta) for every theta in IMPORTANT_THETAS
        ! TODO: Because IMPORTANT_THETAS is big, this next line can induce a stack overflow.
        lambda_segments = lambdafast(important_thetas,simplified_moments)
        ! INRANGE = .TRUE. if a particular value of theta satisfies lambda_range(1) <= lambda(theta) <= lambda_range(2)
        ! Notice that we have put a little error tolerance in here, since zbrent won't find
        ! the exact root.  TODO: Make sure the tolerance here is big enough for the error in zbrent.
        inrange = ((lambda_segments .ge. lambda_range(1)-0.001_dp) .and. (lambda_segments .le. lambda_range(2)+0.001_dp))
        if (k > 1) inrange(1:(k-1)) = .TRUE.
        ! If the lowest value in IMPORTANT_THETAS is in range, then there is no (finite) lower bound
        if (inrange( maxval(minloc(important_thetas)))) then
            estimate_theta(1,1) = -internal_infinity
        else
        ! Otherwise the the lower bound for theta is the minimum value in IMPORTANT_THETAS that is in range
            estimate_theta(1,1) = minval(important_thetas, inrange)
        end if
        ! If the highest value in IMPORTANT_THETAS is in range, then there is no (finite) upper bound
        if (inrange( maxval(maxloc(important_thetas)))) then
            estimate_theta(2,1) = internal_infinity
        else
        ! Otherwise the the upper bound for theta is the maximum value in IMPORTANT_THETAS that is in range
            estimate_theta(2,1) = maxval(important_thetas, inrange)
        end if
        ! Now we find the gradient
        ! Take the gradient at both theta_L and theta_H
        do j=1,2
            theta = estimate_theta(j,1)
            ! The gradient can only be calculated if theta is finite!
            if (is_finite(theta)) then
                ! Gradients are estimated using a simple finite central difference, i.e.
                !           df/dx = (f(x+e)-f(x-e))/2e
                ! where e is some small step size.  The tricky part is getting the right
                ! step size.  The algorithm used here is an adaptation of dfridr in
                ! Numerical Recipes.  However, that algorithm needs an input initial step
                ! size h.
                !
                ! http://www.fizyka.umk.pl/nrbook/c5-7.pdf: "As a function of input h, it is typical for the
                ! accuracy to get better as h is made larger, until a sudden point is reached where nonsensical
                ! extrapolation produces early return with a large error. You should therefore choose
                ! a fairly large value for h, but monitor the returned value err, decreasing h if it is
                ! not small. For functions whose characteristic x scale is of order unity, we typically
                ! take h to be a few tenths."
                !
                ! So we try out starting values (h) until we get one that gives an acceptable
                ! estimated error.
                !
                do n=1,nmax
                    ! Our candidate initial stepsize is 0.1, 0.001, ... 0.0000000001.
                    h = 0.1_dp**n
                    ! Initialize errmax
                    errmax = 0.0_dp
                    ! Initialize the finite-difference vector
                    deps = 0.0_dp
                    ! First, we calculate the scalar-as-vector (dlambda / dtheta)
                    ! hh is the current step size.
                    hh = h
                    ! Calculate an approximate derivative using stepsize hh
                    a(1,1) = maxval((lambdafast((/ (theta+hh) /),simplified_moments)- lambdafast((/ (theta-hh) /),simplified_moments))/(2.0_dp*hh))
                    ! Set the error to very large
                    err = big
                    ! Generate a geometric series
                    fac(1:ntab-1)=geop(con2,con2,ntab-1)
                    ! Now we try progressively smaller stepsizes
                    do k=2,ntab
                        ! The new stepsize hh is the old stepsize divided by 1.4
                        hh=hh/con
                        ! Calculate an approximate derivative with the new stepsize
                        a(1,k) = maxval((lambdafast((/ (theta+hh) /),simplified_moments)- lambdafast((/ (theta-hh) /),simplified_moments))/(2.0_dp*hh))
                        ! Then use Neville's method to estimate the error
                        do m=2,k
                            a(m,k) = (a(m-1,k)*fac(m-1)-a(m-1,k-1))/(fac(m-1)-1.0_dp)
                        end do
                        errt(1:k-1) = max(abs(a(2:k,k)-a(1:k-1,k)),abs(a(2:k,k)-a(1:k-1,k-1)))
                        ierrmin=iminloc(errt(1:k-1))
                        ! If the approximation error is lower than any previous, use that value
                        if (errt(ierrmin) <= err) then
                            err = errt(ierrmin)
                            dfridr = a(1+ierrmin,k)
                        end if
                        ! If the error has increased by a large amount, stop trying new stepsizes
                        if (abs(a(k,k)-a(k-1,k-1)) >= safe*err) exit
                    end do
                    ! errmax is the biggest approximation error so far for the current value of h
                    errmax = max(errmax,err)
                    ! Now we have a candidate derivative dlambda/dtheta
                    dtheta = dfridr
                    ! Second, estimate the vector (dlambda / dmoment_vector)
                    do i=1,size(moment_vector)
                        hh = h
                        deps(i) = hh
                        a(1,1) = (lambdafun(moment_vector+deps,theta) - lambdafun(moment_vector-deps,theta))/(2.0_dp*hh)
                        err = big
                        fac(1:ntab-1)=geop(con2,con2,ntab-1)
                        do k=2,ntab
                            hh=hh/con
                            deps(i) = hh
                            a(1,k) = (lambdafun(moment_vector+deps,theta) - lambdafun(moment_vector-deps,theta))/(2.0_dp*hh)
                            do m=2,k
                                a(m,k) = (a(m-1,k)*fac(m-1)-a(m-1,k-1))/(fac(m-1)-1.0_dp)
                            end do
                            errt(1:k-1) = max(abs(a(2:k,k)-a(1:k-1,k)),abs(a(2:k,k)-a(1:k-1,k-1)))
                            ierrmin=iminloc(errt(1:k-1))
                            if (errt(ierrmin) <= err) then
                                err = errt(ierrmin)
                                dfridr = a(1+ierrmin,k)
                            end if
                            if (abs(a(k,k)-a(k-1,k-1)) >= safe*err) exit
                        end do
                        ! errmax is the biggest approximation error so far for the current value of h
                        errmax = max(errmax,err)
                        dmoments(i) = dfridr
                        deps(i) = 0.0_dp
                    end do
                    ! At this point we have estimates of the derivatives stored in dtheta and dmoments
                    ! We also have the maximum approximation error for the current h stored in errmax
                    ! If that approximation error is "good enough" we are done and can exit the loop
                    if (errmax < 0.01_dp) exit
                    ! Otherwise we will try again with a smaller h
                    if (n .eq. nmax) then
                        call write_to_logfile("Warning: Inaccurate SE for thetaL/H. Try normalizing variables to mean zero.")
                    end if
                end do
                ! Finally, we apply the implicit function theorem to calculate the gradient
                ! that we actually need:   dtheta/dmoments = -(dlambda/dmoments)/(dlambda/dtheta)
                estimate_theta(j,2:) = -dmoments/dtheta(1)
            else
                ! If theta is infinite, then the gradient is NaN.
                estimate_theta(j,2:) = internal_nan
            end if
        end do
    contains
        !---------------------------------------------------------------------------
        ! Functions and subroutines used by ESTIMATE_THETA
        !---------------------------------------------------------------------------
    function lambda_minus_lambda(theta,simplified_moments_and_lambda)
        real(kind=DP), intent(in) :: theta
        real(kind=DP), dimension(:), intent(in) :: simplified_moments_and_lambda
        real(kind=DP) :: lambda_minus_lambda
        ! Potential FPE
        lambda_minus_lambda = minval(lambdafast( (/ theta /),simplified_moments_and_lambda(2:7))) - simplified_moments_and_lambda(1)
    end function lambda_minus_lambda
end function estimate_theta


    !---------------------------------------------------------------------------
    ! ESTIMATE_THETA_SEGMENTS subroutine
    !
    ! Description: Estimates the critical points of the lambda(theta) function.
    !              This information is needed by ESTIMATE_THETA in order to
    !              invert lambda(theta).
    !
    ! Usage: estimate_theta_segments(moment_vector,thetastar)
    !
    ! Arguments:
    !
    !       moment_vector       The vector of cross-moments (input)
    !       thetastar           The value of the parameter thetastar
    !
    ! Effect: The global real array THETA_SEGMENTS is created.
    !
    !---------------------------------------------------------------------------
    subroutine estimate_theta_segments(moment_vector,thetastar)
        real(kind=DP), dimension(:), intent(in) :: moment_vector
        real(kind=DP), intent(in) :: thetastar
        integer, parameter :: k=30000 ! For some reason a bigger number produces an FP overflow
        real(kind=DP), dimension(k) :: thetavec, lambdavec
        logical, dimension(k) :: localmin, localmax
        real(kind=DP) :: thetamax, tmptheta
        real(kind=DP), dimension(6) :: sm
        integer i,j,u, ios, tmp(2)
        logical :: file_found
        sm = simplify_moments(moment_vector)
        !
        ! THETAMAX is the largest value of theta for which we can calculate both lambda(theta) and lambda(-theta)
        ! without generating a floating point exception.
        !
        thetamax = sqrt(huge(1.0_dp)/max(1.0_dp,sm(5),sm(2)-sm(5)))
        ! The calculation above seems clever, but it turns out not to always work.  So I've put in a hard limit as well
        thetamax = minval( (/1.0e100_dp , thetamax /) )
        ! Create a starting set of theta values at which to calculate lambda
        thetavec  = (/ -thetamax , seq(-50.0_dp,50.0_dp,k-2), thetamax /)
        if (is_finite(thetastar)) then
            ! Figure out where thetastar lies in thetavec
            i = count(thetavec < thetastar)
            ! If i=0 or i=k, then thetastar is finite but outside of [-thetamax,thetamax].
            ! This is unlikely, but we should check anyway.
            if ( ( (i > 0) .and. (i < k))) then
                ! Adjust i to ensure that -thetamax and thetamax are still included in thetavec
                i = min(max(i,2),k-2)
                ! Replace the two elements of thetavec that bracket thetastar with
                ! two more carefully-chosen numbers.  See BRACKET_THETA_STAR for details
                thetavec(i:i+1) = bracket_theta_star(moment_vector)
                ! TODO: There is a potential bug here.  The bracket_theta_star function is
                ! used to take the two values in thetavec that are closest to thetastar
                ! and replace them with values that are guaranteed to give finite and
                ! nonzero lambda.  But there's nothing to guarantee that these are still
                ! the two values in thetavec that are the closest to thetastar.
            else
                write (unit=buffer,fmt=*) "Warning: absolute value of thetastar (", thetastar, ") is larger than thetamax (", thetamax, ")."
                call write_to_logfile(trim(buffer(2:)))
            end if
        end if
        ! Re-sort thetavec
        thetavec = sort(thetavec)
        ! Calculate lambda for every theta in thetavec
        lambdavec = lambdafast(thetavec,simplify_moments(moment_vector))
        ! If a detail_file has been specified, output thetavec and lambdavec to that file
        if (len_trim(detail_file) > 0) then
            ! Get a new unit number
            u = get_next_unit()
            ! See if the file exists already
            inquire(file=detail_file,exist=file_found)
            ! Open the file
            if (file_found) then
                open(unit=u,file=detail_file,iostat=ios,action="write",position="rewind",status="old")
            else
                open(unit=u,file=detail_file,iostat=ios,action="write",status="new")
            end if
            ! Assuming the file could be opened correctly,...
            if (ios == 0) then
                write (unit=u,fmt=*) "theta , lambda"
                ! Write out thetavec
                do i=1,size(thetavec)
                    write (unit=u,fmt=*) thetavec(i) , "," , lambdavec(i)
                end do
                ! Close the file
                close (unit=u,iostat=ios)
            else
                ! If the file can't be opened correctly, then put a warning message in the logfile
                call write_to_logfile("WARNING: Unable to write detailed lambda(theta) function to detail file " // trim(detail_file))
            end if
        end if
        ! LOCALMIN=.TRUE if the corresponding element of THETAVEC appears to be a local minimum
        localmin = (/ .false. , ((lambdavec(2:(k-1)) < lambdavec(1:(k-2))) .and. (lambdavec(2:(k-1)) < lambdavec(3:k))) , .false. /)
        ! LOCALMAX=.TRUE if the corresponding element of THETAVEC appears to be a local maximum
        localmax = (/ .false. , ((lambdavec(2:(k-1)) > lambdavec(1:(k-2))) .and. (lambdavec(2:(k-1)) > lambdavec(3:k))) , .false. /)
        if (is_finite(thetastar)) then
            ! Figure out where THETASTAR lies in THETAVEC.  We need to do this calculation again because we sorted THETAVEC
            i = count(thetavec < thetastar)
            if ( ( (i > 0) .and. (i < k))) then
                ! The two values bracketing THETASTAR are never local optima
                localmin(i:i+1) = .false.
                localmax(i:i+1) = .false.
            end if
        end if
        ! Right now, we only have approximate local optima.  We need to apply
        ! an iterative optimization algorithm to improve the precision.
        do j=1,size(localmin)
            if (localmin(j) == .TRUE.) then
                thetavec(j) = brent(thetavec(j-1),thetavec(j),thetavec(j+1),lambda_for_brent,1.0e-10_dp,simplify_moments(moment_vector))
            else if (localmax(j) == .TRUE.) then
                thetavec(j) = brent(thetavec(j-1),thetavec(j),thetavec(j+1),negative_lambda_for_brent,1.0e-10_dp,simplify_moments(moment_vector))
            end if
        end do
        ! Now we are ready to create THETA_SEGMENTS.  We didn't know how big it was, so we needed to
        ! allocate it dynamically
        if (is_finite(thetastar) .and. (i > 0) .and. (i < k) ) then
            allocate(theta_segments(count(localmin)+count(localmax)+4))
            ! THETA_SEGMENTS contains the two limits (-Inf,+Inf), the pair of values that bracket thetastar, and any local optima
            theta_segments = (/ -thetamax , pack(thetavec,localmin) , pack(thetavec,localmax) , thetavec(i:i+1), thetamax /)
        else
            ! If thetastar is not finite, then we have two less elements in THETA_SEGMENTS
            allocate(theta_segments(count(localmin)+count(localmax)+2))
                theta_segments = (/ -thetamax , pack(thetavec,localmin) , pack(thetavec,localmax) , thetamax /)
        end if
        ! Sort the result (definitely necessary)
        theta_segments = sort(theta_segments)
    end subroutine estimate_theta_segments



    !---------------------------------------------------------------------------
    ! ESTIMATE_PARAMETER function
    !
    ! Description: Calculates a parameter estimate and its gradient
    !
    ! Usage: estimate_parameter(func,moment_vector,xopt)
    !
    ! Arguments:
    !
    !       func                A scalar-valued function of the form func(x,xopt)
    !                           where x is a vector, and xopt is
    !                           an optional vector.  Note that the function
    !                           must accept xopt, and it must be optional
    !                           in that function.
    !
    !       moment_vector       The vector of cross-moments (input)
    !
    !
    !
    ! Effect: The function returns a vector of length = length(moment_vector)+1.
    !
    !       estimate_parameter(1) is the actual parameter estimate, i.e.
    !           func(moment_vector)
    !
    !       estimate_parameter(2:) is the gradient, i.e.
    !
    !           dfunc(moment_vector) / dmoment_vector
    !
    !---------------------------------------------------------------------------
    function estimate_parameter(func,moment_vector,fname)
        real(kind=DP), dimension(:), intent(in) :: moment_vector
        character(len=*), intent(in) :: fname
        real(kind=DP), dimension(size(moment_vector)+1) :: estimate_parameter
        real(kind=DP) :: h=1.0e-4_dp, err, dfridr, hh, errmax
        real(kind=DP), dimension(size(moment_vector)) :: deps
        real(kind=DP), parameter :: con=1.4_dp, con2=con*con, big=huge(1.0_dp), safe=2.0_dp
        integer, parameter :: ntab=10, nmax=10
        integer :: i,ierrmin, j, k, n
        real(kind=DP), dimension(ntab-1) :: errt, fac
        real(kind=DP), dimension(ntab,ntab) :: a
        ! The function takes this form
        interface
            function func(x)
                use rcrlib, only : SP, DP
                real(kind=DP), dimension(:), intent(in) :: x
                real(kind=DP)  :: func
            end function func
        end interface
        ! Estimate the parameter
        estimate_parameter(1) = func(moment_vector)
        ! Then estimate its derivative using a finite difference approximation
        ! See the comments in ESTIMATE_THETA for a detailed description.
        if (is_finite(estimate_parameter(1))) then
            do n=1,nmax
                h = 0.1_dp**n
                errmax = 0.0_dp
                ! We are estimating the gradient, i.e., a vector of derivatives the same size as moment_vector
                do i=1,size(moment_vector)
                    ! Initialize DEPS
                    deps = 0.0_dp
                    ! HH is the step size.  It is chosen by an algorithm borrowed from the dfridr function
                    ! in Numerical Recipes.  We start with HH set to a predetermined value H.  After that, each
                    ! successive value of HH is the previous value divided by CON (which is set to 1.4)
                    hh = h
                    ! Set element i of DEPS to HH.
                    deps(i)=hh
                    ! Calculate the first approximation
                    a(1,1) = (func(moment_vector+deps)-func(moment_vector-deps))/(2.0_dp*hh)
                    ! The error is assumed to be a big number
                    err = big
                    ! Generate a geometric series
                    fac(1:ntab-1)=geop(con2,con2,ntab-1)
                    ! Try a total of NTAB different step sizes
                    do j=2,ntab
                        ! Generate the next step size
                        hh=hh/con
                        ! Set DEPS based on that step size
                        deps(i)=hh
                        ! Calculate the approximate derivative for that step size
                        a(1,j) = (func(moment_vector+deps)-func(moment_vector-deps))/(2.0_dp*hh)
                        ! Next we estimate the approximation error for the current step size
                        do k=2,j
                            a(k,j) = (a(k-1,j)*fac(k-1)-a(k-1,j-1))/(fac(k-1)-1.0_dp)
                        end do
                        errt(1:j-1) = max(abs(a(2:j,j)-a(1:j-1,j)),abs(a(2:j,j)-a(1:j-1,j-1)))
                        ierrmin=iminloc(errt(1:j-1))
                        ! If the error is smaller than the lowest previous error, use that hh
                        if (errt(ierrmin) <= err) then
                            err = errt(ierrmin)
                            dfridr = a(1+ierrmin,j)
                        end if
                        ! If the error is much larger than the lowest previous error, stop
                        if (abs(a(j,j)-a(j-1,j-1)) >= safe*err) exit
                    end do
                    errmax = max(errmax,err)
                    estimate_parameter(i+1) = dfridr
                end do
                if (errmax < 0.01_dp) exit
                if (n .eq. nmax) then
                    call write_to_logfile("Warning: Inaccurate SE for " // fname // ". Try normalizing variables to mean zero.")
                end if
            end do
        else
            estimate_parameter(2:) = internal_nan
        end if
    end function estimate_parameter

    !---------------------------------------------------------------------------
    ! WRITE_RESULTS subroutine
    !
    ! Description: Writes the results out to a text file
    !
    ! Usage: write_results(outfile)
    !
    !           outfile     The name of the file to write to
    !
    !---------------------------------------------------------------------------
    subroutine write_results(result_matrix,outfile)
        ! result_matrix is the name of the matrix containing our results
        real(kind=DP), dimension(:,:), intent(inout) :: result_matrix
        ! Outfile is the name of the output file
        character(len=*), intent(in) :: outfile
        ! We will translate result_matrix into stata format, and put the result in OUTPUT_MATRIX
        real(kind=DP), dimension(size(result_matrix,1),size(result_matrix,2)) :: output_matrix
        ! In writing the results out, we need a string buffer.  It needs to be big enough to
        ! accomodate all of the results.
        character(len=size(result_matrix,2)*40) :: buffer
        ! File_found is a logical we will use to inquire about OUTFILE's existence
        logical :: file_found
        ! IOS is used for I/O exception handling, I is a counter, and U is used to number file units
        integer :: ios,i,u
        ! Start by translating our results into stata format - Stata can't handle Infinity or NaN.
        output_matrix = translate_to_stata_form(result_matrix)
        ! Put some information in the log file
        call write_to_logfile("Writing results to output file " // outfile // ".")
        write (buffer,iostat=ios,fmt=*) "Actual results = ", result_matrix
        call write_to_logfile(trim(buffer))
        write (buffer,iostat=ios,fmt=*) "Results sent to Stata = ", output_matrix
        call write_to_logfile(trim(buffer))
        ! Now we are ready to write the results to OUTFILE. Start by getting a fresh unit number
        u = get_next_unit()
        ! See if OUTFILE already exists
        inquire(file=outfile,exist=file_found)
        if (file_found) then
            ! If it exists, replace it
            open(unit=u,file=outfile,iostat=ios,action="write",position="rewind",status="old",recl=100000)
        else
            ! If it doesn't exist, create a new one
            open(unit=u,file=outfile,iostat=ios,action="write",status="new",recl=100000) ! if it doesn't exist, create a new one
        end if
        ! Assuming the file could be opened correctly,...
        call validate(ios==0,msg="Error: Unable to open output file " // outfile // ".")
        ! Write out the data one row at a time
        do i=1,size(output_matrix,1)
            ! Start by writing the data to the string buffer
            write (buffer,iostat=ios,fmt=*) output_matrix(i,:)
            ! Then write the buffer to the file, taking out extra whitespace
            write (unit=u,iostat=ios,fmt='(a)') trim(remove_duplicates(trim(buffer)," "))
            call validate(ios==0,msg="Error: Unable to write to output file " // outfile // ".")
        end do
        ! Close the file
        close(unit=u,iostat=ios)
        call write_to_logfile("RCR successfully concluded.")
    contains
        !---------------------------------------------------------------------------
        ! The TRANSLATE_TO_STATA_FORM function takes a real matrix and translates
        ! it to a form that can be handled by Stata
        !---------------------------------------------------------------------------
        function translate_to_stata_form(result_matrix) result(rm)
            real(kind=DP), dimension(:,:), intent(in) :: result_matrix
            real(kind=DP), dimension(size(result_matrix,1),size(result_matrix,2)) :: rm
            rm = result_matrix
            ! Replace infinity (or any other number larger than external_big_number) with external_big_number
            where (abs(rm) .ge. external_big_number) rm = sign(external_big_number,rm)
            ! Replace NaN with zero.
            where (isnan(rm))  rm = 0.0000_dp
        end function translate_to_stata_form
        !---------------------------------------------------------------------------
        ! The REMOVE_DUPLICATES function removes duplicate characters (i.e., spaces)
        ! from a string.
        !---------------------------------------------------------------------------
        function remove_duplicates(str,chr,padchr)
            character(len=*), intent(in) :: str
            character(len=len(str)) :: remove_duplicates,tmpstr
            character, intent(in) :: chr
            character, intent(in), optional :: padchr
            character :: z
            integer :: i,j,k,l
            tmpstr=str
            if (present(padchr)) then
                z = padchr
            else
                z = "Z"
            end if
            ! Check to make sure we don't have any Z's in the string
            j = index(tmpstr,z)
            if (j > 0) then
                call die("Error in remove_duplicates.  String should not have '" // z // "' in it")
            end if
            ! Find out the longest sequence of duplicates
            do i=1,len(str)
                j = index(tmpstr,repeat(chr,i))
                if (j == 0) then
                    j = i-1
                    exit
                end if
            end do
            ! Working backwards from there, replace with a single copy
            do i=j,2,-1
                do l=1,len(str)
                    k = index(tmpstr,repeat(chr,i))
                    if (k==0) then
                        exit
                    else
                        ! Note that we can't make tmpstr shorter, so we pad the end with an arbitrary character (Z)
                        tmpstr=tmpstr(:k) // tmpstr(k+i:) // repeat(z,i-1)
                    end if
                end do
            end do
            ! Now we look for the first "Z"
            k = index(tmpstr,z)
            ! And we clear out everything from that point
            remove_duplicates=tmpstr(:k-1)
        end function remove_duplicates
    end subroutine write_results

    !---------------------------------------------------------------------------
    ! SIMPLIFY_MOMENTS function
    !
    ! Description: Converts the large vector moment_vector into the 6 statistics
    !               we need to estimate the model
    !
    ! Usage: simplify_moments(moment_vector)
    !
    !           moment_vector       The vector of moments
    !
    ! Effect: Outputs a vector of length 6 containing the simplified moments
    !
    ! Note: This is a very inefficient function.  But making it more efficient
    !       would take a lot of effort, and might make it harder to read.
    !---------------------------------------------------------------------------
    function simplify_moments(moment_vector)
        real(kind=DP), dimension(:), intent(in) :: moment_vector
        real(kind=DP), dimension(size(moment_vector)+1) :: mvtmp
        real(kind=DP), dimension(6) :: simplify_moments
        real(kind=DP), dimension(:,:), allocatable :: XX, xtmp,XY,XZ
        integer :: i,j,k,l,m
        ! Get sizes
        m = size(moment_vector)
        k = 1 + floor((sqrt(real(1+8*m,kind=DP))-1.0_dp)/2.0_dp)
        ! Allocate the necessary arrays
        allocate(xtmp(k,k),XX(k-2,k-2),XY(k-2,1),XZ(k-2,1))
        mvtmp = (/ 1.0_dp , moment_vector /)
        ! The array XTMP will contain the full cross-product matrix E(WW') where W = [1 X Y Z]
        l=1
        do i=1,k
            do j = i,k
                xtmp(i,j) = mvtmp(l)
                xtmp(j,i) = mvtmp(l)
                l=l+1
            end do
        end do
        ! The array XX will contain the symmetric matrix E(XX')
        XX = xtmp(1:size(XX,1),1:size(XX,2))
        ! The array XY will contain the vector E(XY)
        XY(:,1) = xtmp(size(xtmp,1)-1,1:size(xtmp,2)-2)
        ! The array XZ will contain the vector E(XZ)
        XZ(:,1) = xtmp(size(xtmp,1),1:size(xtmp,2)-2)
        ! Now we fill in simplify_moments with the various moments.
        ! varY
        simplify_moments(1) = moment_vector(m-2) - (moment_vector(k-2))**(2.0_dp)
        ! varZ
        simplify_moments(2) = moment_vector(m) - (moment_vector(k-1))**(2.0_dp)
        ! covYZ
        simplify_moments(3) = moment_vector(m-1) - moment_vector(k-1)*moment_vector(k-2)
        ! varYhat
        simplify_moments(4) = maxval(matmul(transpose(XY),bksolve(XX,XY)))  - (moment_vector(k-2))**(2.0_dp)
        ! varZhat
        simplify_moments(5) = maxval(matmul(transpose(XZ),bksolve(XX,XZ)))  - (moment_vector(k-1))**(2.0_dp)
        ! covYZhat
        simplify_moments(6) = maxval(matmul(transpose(XY),bksolve(XX,XZ)))  - moment_vector(k-1)*moment_vector(k-2)
        ! When there is only one control variable yhat and zhat are perfectly correlated (positively or negatively)
        ! With rounding error, this can lead to a correlation that is > 1 in absolute value.  This can
        ! create problems, so we force the correlation to be exactly 1.
        ! TODO: This also could happen if there is more than one control variable but only one happens
        ! to have a nonzero coefficient.  I don't know how to handle that case.
        if (k .eq. 4) then
            simplify_moments(6) = sign(sqrt(simplify_moments(4)*simplify_moments(5)),simplify_moments(6))
        end if
        ! One thing to notice is that matmul produces an array (even if it is a 1-by-1 array)
        ! So we took maxval of each array, which converts the array to a scalar.
        ! That seems safe, but maybe it can be improved later.
        !
        ! Just to be sure, deallocate the local variables.  The compiler probably does this automatically.
        deallocate(xtmp,XX,XY,XZ)
    end function simplify_moments

    !---------------------------------------------------------------------------
    ! LAMBDAFUN function
    ! LAMBDAFAST function
    !
    ! Description: Calculates the lambda(theta) function. LAMBDAFAST is the main
    !              function here.  LAMBDAFUN is just a wrapper that takes inputs
    !              and produces outputs in a slightly different format that is
    !              useful for calculating gradients.
    !
    ! Usage: lambdafun(moment_vector,theta)
    !        lambdafast(theta,simplified_moments)
    !
    !           moment_vector       The vector of moments
    !           theta               A scalar (for lambdafun) or vector (for lambdafast)
    !                               values of theta at which to calculate lambda(theta)
    !           simplified_moments  The vector of simplified moments
    !
    ! Effect: A scalar (for lambdafun) or vector (for lambdafast) corresponding to
    !         lambda(theta) for the input theta.
    !
    !---------------------------------------------------------------------------
    function lambdafun(moment_vector,theta)
        real(kind=DP), dimension(:), intent(in) :: moment_vector
        real(kind=DP), intent(in) :: theta
        real(kind=DP) :: lambdafun
        ! Potential FPE
        lambdafun = maxval(lambdafast( (/ theta /),simplify_moments(moment_vector)))
    end function lambdafun

    function lambda0(moment_vector)
        real(kind=DP), dimension(:), intent(in) :: moment_vector
        real(kind=DP) :: lambda0
        ! lambda0 is defined as:
        ! (cov(y,z)/cov(yhat,zhat)-1) / sqrt(var(y)/var(yhat)-1)
        !
        ! The check_moments subroutine should ensure that
        !  var(y) >= var(yhat) > 0, so the denominator is
        ! always positive and finite.
        !
        ! Special values: If cov(yhat,zhat)=0, then lambda0 can
        !   be +Infinity, -Infinity, or NaN depending on the sign
        !   of cov(y,z).
        !
        lambda0 = lambdafun(moment_vector,0.0_dp)
    end function lambda0


    function lambdafast(theta,simplifiedMoments)
        real(kind=DP), dimension(:), intent(in) :: theta
        real(kind=DP), dimension(6), intent(in) :: simplifiedMoments
        real(kind=DP), dimension(size(theta)) :: lambdafast
        real(kind=DP) :: y,z,yz,yhat,zhat,yzhat
        y = simplifiedMoments(1)
        z = simplifiedMoments(2)
        yz = simplifiedMoments(3)
        yhat = simplifiedMoments(4)
        zhat = simplifiedMoments(5)
        yzhat = simplifiedMoments(6)
        ! Potential FPE
                lambdafast = (yhat - 2.0_dp*theta*yzhat + theta**2.0_dp * zhat)
                lambdafast = lambdafast/(y - yhat - (2.0_dp)*theta*(yz-yzhat) + theta**(2.0_dp)*(z-zhat))
        lambdafast = (yz - yzhat - theta*(z-zhat))/(yzhat - theta*zhat)*sqrt(lambdafast)
    end function lambdafast


    function lambda_for_brent(theta,simplifiedMoments)
        real(kind=DP), intent(in) :: theta
        real(kind=DP), dimension(6), intent(in) :: simplifiedMoments
        real(kind=DP) :: lambda_for_brent
        lambda_for_brent = maxval(lambdafast( (/ theta /), simplifiedMoments))
    end function lambda_for_brent

    function negative_lambda_for_brent(theta,simplifiedMoments)
        real(kind=DP), intent(in) :: theta
        real(kind=DP), dimension(6), intent(in) :: simplifiedMoments
        real(kind=DP) :: negative_lambda_for_brent
        negative_lambda_for_brent = -lambda_for_brent(theta,simplifiedMoments)
    end function negative_lambda_for_brent

        !---------------------------------------------------------------------------
        ! BRACKET_THETA_STAR returns a pair of numbers, one of which is slightly
        !                    less than thetastar and one that is slightly more.
        !                    The numbers are chosen to be as small as possible while
        !                    avoiding major floating-point approximation error.
        !---------------------------------------------------------------------------
    function bracket_theta_star(moment_vector) result(bracket)
        real(kind=DP), dimension(:), intent(in) :: moment_vector
        real(kind=DP), dimension(:) :: bracket(2),true_limit(2),candidate(2),sm(6),tmp2(2)
        real(kind=DP) :: theta_star
        integer :: i,j
        ! Get the value of THETASTAR.  If we are in this function it should be finite.
        theta_star = thetastar(moment_vector)
        ! Get the limit of lambda(theta) as theta approaches THETASTAR,from below and from above
        ! These limits are generally not finite.
        sm = simplify_moments(moment_vector)
        true_limit = (/1.0_dp, -1.0_dp /)* internal_big_number * sign(1.0_dp,sm(3) - sm(6) * sm(2)/sm(5))
        ! Pick a default value
        bracket = theta_star + (/-1.0_dp, 1.0_dp/)*max(abs(theta_star),1.0_dp)*0.1_dp
        j = 0
        do i=1,100
            ! For the candidate bracket, consider THETASTAR plus or minus some small number epsilon
            ! (epsilon gets smaller each iteration)
            candidate = theta_star + (/-1.0_dp, 1.0_dp/)*max(abs(theta_star),1.0_dp)*0.1_dp**real(i,kind=DP)
            ! In order to be a good bracket, candidate must satisfy some conditions.
            !    1. The bracket must be wide enough that the system can tell that CANDIDATE(1) < THETASTAR < CANDIDATE(2)
            !    2. The bracket must be narrow enough that lambda(candidate) is the same sign as true_limit.
            !    3. The bracket must be wide enough that lambda(candidate) is finite and nonzero.
            !       If candidate is very close to thetastar, then the calculated lambda(candidate)
            !       can be *either* NaN or zero.  The reason for this is that lambda(candidate) is
            !       a ratio of two things that are going to zero.  Approximation error will eventually
            !       make both the numerator and denominator indistingushable from zero (NaN), but sometimes
            !       the numerator will reach indistinguishable-from-zero faster (giving zero for the ratio).
            !
            if ((candidate(1) < theta_star) .and. (candidate(2) > theta_star)) then
                tmp2 = lambdafast(candidate,sm)
                if (all(is_finite(tmp2)) .and. (tmp2(1)*sign(1.0_dp,true_limit(1)) > 0.0_dp) .and. (tmp2(2)*sign(1.0_dp,true_limit(2)) > 0.0_dp) ) then
                    j = i
                    bracket = candidate
                end if
            else
                exit
            end if
        end do
        if (j == 0) then
            call write_to_logfile("Warning: Unable to find a good bracket for thetastar")
        end if
    end function bracket_theta_star

    !---------------------------------------------------------------------------
    ! DIE subroutine
    !
    ! Description: Closes the program with an error message
    !
    ! Usage: die(msg)
    !
    ! Inputs
    !        msg              a character string of any length
    !
    ! Effect: Writes the string MSG to standard error and closes the program.
    !         In addition, if possible, this subroutine will pass a message
    !         to Stata that the program exited with an error.  This takes
    !         the form of just writing the word "ERROR" to the OUTFILE.
    !
    !---------------------------------------------------------------------------
    subroutine die(msg)
        character(len=*), intent(in), optional :: msg
        integer :: i,u,ios
        logical :: file_found
        ! We need to tell Stata there's been a problem.  We will do this
        ! by writing the word "ERROR" to the output file.
        ! Get a new unit number
        u = get_next_unit()
        ! See if OUTFILE already exists
        inquire(file=outfile,exist=file_found)
        if (file_found) then
            ! If it exists, replace it
            open(unit=u,file=outfile,iostat=ios,action="write",position="rewind",status="old")
        else
            ! If it doesn't exist, create a new one
            open(unit=u,file=outfile,iostat=ios,action="write",status="new") ! if it doesn't exist, create a new one
        end if
        ! Assuming the file could be opened correctly, add the error notation to the file
        write (unit=u,iostat=ios,fmt=*) "ERROR"
        ! Close the file
        close(unit=u,iostat=ios)
        ! Write the error message to standard error.
        if (present(msg)) then
            write (unit=stderr,iostat=ios,fmt=*) msg
        else
            write (unit=stderr,iostat=ios,fmt=*) "Program has terminated with an unknown error. "
        end if
        ! Exit the program
        stop
    end subroutine die

    !---------------------------------------------------------------------------
    ! CHECK_MOMENTS subroutine
    !
    ! Description: Checks to see if the moments in MOMENT_MATRIX are logically
    !              consistent
    !
    ! Usage: check_moments(moment_vector)
    !
    ! Arguments:
    !
    !       moment_vector       The vector of cross-moments (input)
    !
    ! Effect: Stops the program with an error message if any of the tested
    !         conditions fail.
    !---------------------------------------------------------------------------
    subroutine check_moments(moment_vector)
        real(kind=DP), dimension(:), intent(in) :: moment_vector
        real(kind=DP), dimension(6) :: sm
        sm = simplify_moments(moment_vector)
!       write (*,*) "var(y) =  " , sm(1)
!       write (*,*) "var(z) = " , sm(2)
!       write (*,*) "cov(y,z) = " , sm(3)
!       write (*,*) "var(yhat) = " , sm(4)
!       write (*,*) "var(zhat) = " , sm(5)
!       write (*,*) "cov(yhat,zhat) = " , sm(6)
        ! First make sure that moment_vector describes a valid covariance matrix
        call validate(sm(1) .ge. 0.0_dp,"Error - invalid data in moment_vector: var(y) < 0")
        call validate(sm(2) .ge. 0.0_dp,"Error - invalid data in moment_vector: var(z) < 0")
        call validate( abs(sm(3)) .le. sqrt(sm(1)*sm(2)) ,"Error - invalid data in moment_vector: cov(y,z) > sqrt(var(y)*var(z))")
        call validate(sm(4) .ge. 0.0_dp,"Error - invalid data in moment_vector: var(yhat) < 0")
        call validate(sm(4) .le. sm(1),"Error - invalid data in moment_vector: var(yhat) > var(y)")
        call validate(sm(5) .ge. 0.0_dp,"Error - invalid data in moment_vector: var(zhat) < 0")
        call validate(abs(sm(6)) .le. sqrt(sm(4)*sm(5)),"Error - invalid moment_vector: cov(yhat,zhat) > sqrt(var(yhat)*var(zhat))")
        ! Next make sure that the identifying conditions are satisfied.
        ! TODO: Maybe these could be addressed with warnings rather than error messages?
        call validate(sm(1) > 0.0_dp,"Error - model not identified: var(y)=0")
        call validate(sm(2) > 0.0_dp,"Error - model not identified: var(z)=0")
        call validate(sm(4) > 0, "Error - model not identified: var(yhat)=0")
        call validate(sm(4) < sm(1),"Error - model not identified: y is an exact linear function of X")
        ! TODO: We may also want to check for var(zhat)=0.  The model is identified in this case, but
        !       we may need to take special steps to get the calculations right.
!       end if
    end subroutine check_moments

    !---------------------------------------------------------------------------
    ! FIND_ROOT function (no longer used)
    !
    ! Description: Finds the root of a function.
    !
    ! Usage: find_root(func,xrange,xopt)
    !
    !           func    A real-valued function of the form FUNC(X,XOPT)
    !           xrange  A length-2 vector indicating the range of
    !                   X over which to look for the root.
    !           xopt    An additional vector argument to FUNC
    !
    ! Effect: Returns the real number ROOT, where:
    !
    !           1. XRANGE(1) <= ROOT <= XRANGE(2)
    !           2. FUNC(ROOT,XOPT) = 0 (approximately)
    !
    !---------------------------------------------------------------------------
    function find_root(func,xrange,xopt) result(root)
        real(kind=DP), dimension(2), intent(in) :: xrange
        real(kind=DP), dimension(:), intent(in) :: xopt
        real(kind=DP) :: xtmp,ytmp,yrange(2), root(2)
        integer :: i,j
        interface
            function func(x,xopt)
                use rcrlib, only : SP,DP
                real(kind=DP), intent(in) :: x
                real(kind=DP), dimension(:), intent(in) :: xopt
                real(kind=DP)  :: func
            end function func
        end interface
        root=xrange
        do i=1,10000
            do j=1,2
                yrange(j) = func(root(j),xopt)
            end do
            ! TODO: This next line could be improved upon.  Suppose root(1)=0 and root(2)=1.8e308.  Then this
            !       line will produce xtmp=9.0e307, then xtmp=4.5e307, etc.  As a result, it will take many iterations to
            !       get to some reasonable range.  That's why I've set the number of iterations at 10,000.
            !       Even that might not be sufficient to get convergence.
            xtmp = sum(root)/2.0_dp
            ytmp = func(xtmp,xopt)
            if (yrange(1)*ytmp > 0) then
                root(1) = xtmp
            else
                root(2) = xtmp
            end if
            if ((root(2) - root(1)) < 100.0_dp*epsilon(1.0_dp)) then
                exit
            end if
        end do
    end function find_root




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! I/O UTILITY FUNCTIONS AND SUBROUTINES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MATH UTILITY FUNCTIONS AND SUBROUTINES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    !---------------------------------------------------------------------------
    ! BETWEEN function
    !
    ! Description: Determines if one real number lies between two others.
    !
    ! Usage: between(x,xrange)
    !
    !           x       A real number
    !           xrange  A length-2 real vector
    !
    ! Effect: Returns .TRUE. if XRANGE(1) <= x <= XRANGE(2).
    !
    !---------------------------------------------------------------------------
    function between(x,xrange)
        real(kind=DP), intent(in) :: x
        real(kind=DP), dimension(2), intent(in) :: xrange
        logical :: between
        between = ( ( x >= xrange(1) ) .and. ( x <= xrange(2) ) )
    end function between



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ESTIMATION FUNCTIONS AND SUBROUTINES (used for ESTIMATE_MODEL)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !---------------------------------------------------------------------------
    ! BRENT function
    !
    ! Description: Finds the minimum of a function.  The code is adapted from
    !              Numerical Recipes, with changes made to fit the needs here.
    !
    ! Usage: brent(ax,bx,cx,func,tol,xopt)
    !
    !           ax      The lower bound of the search range
    !           bx      The upper bound of the search range.
    !           cx      A point between ax and bx such that func(cx) < func(ax) and
    !                   func(cx) < func(bx)
    !           func    The function to optimize. It should have form func(x,xopt)
    !                   where x is a DP scalar and xopt is a length-6 DP vector.
    !           tol     A scalar measuring the tolerance.
    !           xopt    The vector of simplified moments.
    !
    ! Result: The value of x between ax and bx that minimizes func(x,xopt).
    !
    !---------------------------------------------------------------------------
        function brent(ax,bx,cx,func,tol,xopt)
        ! use nrtype; use nrutil, only : nrerror
        real(kind=DP), intent(in) :: ax,bx,cx,tol
        real(kind=DP), dimension(:), intent(in) :: xopt
        real(kind=DP) :: brent
        interface
            function func(x,xopt)
                use rcrlib, only : SP,DP
                real(kind=DP), intent(in) :: x
                real(kind=DP), dimension(6), intent(in) :: xopt
                real(kind=DP) :: func
            end function func
        end interface
        integer, parameter :: itmax=1000
        real(kind=DP), parameter :: cgold=0.3819660_dp, zeps=1.0e-3_dp*epsilon(ax)
        integer :: iter
        real(kind=DP) :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
        a = min(ax,cx)
        b = max(ax,cx)
        v=bx
        w=v
        x=v
        e=0.0_dp
        fx = func(x,xopt)
        fv = fx
        fw = fx
        do iter=1,itmax
            xm = 0.5_dp*(a+b)
            tol1 = tol*abs(x)+zeps
            tol2 = 2.0_dp*tol1
            if ( abs(x-xm) <= (tol2 - 0.5_dp*(b-a))) then
                brent=x
                return
            end if
            if (abs(e) > tol1) then
                r = (x-w)*(fx-fv)
                q = (x-v)*(fx-fw)
                p = (x-v)*q - (x-w)*r
                q = 2.0_dp*(q-r)
                if (q > 0.0_dp) then
                    p = -p
                end if
                q = abs(q)
                etemp = e
                e = d
                if (abs(p) >= abs(0.5_dp*q*etemp) .or. p <= q*(a-x) .or. p >= q*(b-x) ) then
                    e = merge(a-x,b-x,x >= xm)
                    d = cgold*e
                else
                    d = p/q
                    u = x+d
                    if (u-a < tol2 .or. b-u < tol2) then
                        d = sign(tol1,xm-x)
                    end if
                end if
            else
                e = merge(a-x,b-x, x >= xm)
                d = cgold*e
            end if
            u = merge(x+d,x+sign(tol1,d),abs(d) >= tol1)
            fu = func(u,xopt)
            if (fu <= fx) then
                if (u >= x) then
                    a = x
                else
                    b = x
                end if
                call shft(v,w,x,u)
                call shft(fv,fw,fx,fu)
            else
                if (u < x) then
                    a = u
                else
                    b = u
                end if
                if ( fu <= fw .or. w == x) then
                    v = w
                    fv = fw
                    w = u
                    fw = fu
                else if (fu <= fv .or. v == x .or. v == w) then
                    v = u
                    fv = fu
                end if
            end if
        end do
        brent = x
        call write_to_logfile("brent: exceed maximum iterations")
        contains
        subroutine shft(a,b,c,d)
            real(kind=DP), intent(out) :: a
            real(kind=DP), intent(inout) :: b,c
            real(kind=DP), intent(in) :: d
            a = b
            b = c
            c = d
        end subroutine shft
    end function brent



    !---------------------------------------------------------------------------
    ! ZBRENT function
    !
    ! Description: Finds the root of a function.  The code is adapted from
    !              Numerical Recipes, with changes made to fit the needs here.
    !
    ! Usage: zbrent(func,x1,x2,tol,xopt)
    !
    !           func    The function to optimize. It should have form func(x,xopt)
    !                   where x is a DP scalar and xopt is a length-6 DP vector.
    !           x1      The lower bound of the search range
    !           x2      The upper bound of the search range.  It is assumed
    !                   that there is exactly one root of func between x1 and x2.
    !           tol     A scalar measuring the tolerance.
    !           xopt    The vector of simplified moments.
    !
    ! Result: A value of x between x1 and x2 such that func(x,xopt) = 0.
    !
    !---------------------------------------------------------------------------
    function zbrent(func,x1,x2,tol,xopt)
        real(kind=DP), intent(in) :: x1,x2,tol
        real(kind=DP), dimension(:), intent(in) :: xopt
        real(kind=DP) :: zbrent
        interface
            function func(x,xopt)
                use rcrlib, only : sp,dp
                real(kind=DP), intent(in) :: x
                real(kind=DP), dimension(:), intent(in) :: xopt
                real(kind=DP) :: func
            end function func
        end interface
        integer, parameter :: itmax=1000
        real(kind=DP), parameter :: eps=epsilon(x1)
        integer :: iter
        real(kind=DP) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
        a = x1
        b = x2
        fa = func(a,xopt)
        fb = func(b,xopt)
        if (( fa > 0.0_dp .and. fb > 0.0_dp) .or. ( fa < 0.0_dp .and. fb < 0.0_dp) ) then
            call write_to_logfile("Error in zbrent: Root is not bracketed")
            call die("root must be bracketed for zbrent")
        end if
        c = b
        fc = fb
        do iter=1,itmax
            if (( fb > 0.0_dp .and. fc > 0.0_dp) .or. ( fb < 0.0_dp .and. fc < 0.0_dp) ) then
                c = a
                fc = fa
                d = b-a
                e = d
            end if
            if (abs(fc) < abs(fb)) then
                a = b
                b = c
                c = a
                fa = fb
                fb = fc
                fc = fa
            end if
            ! check for convergence
            tol1=2.0_dp*eps*abs(b)+0.5_dp*tol
            xm = 0.5_dp*(c-b)
            if (abs(xm) <= tol1 .or. fb == 0.0_dp) then
                zbrent = b
                return
            end if
            if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
                s = fb/fa
                if (a == c) then
                    p = 2.0_dp*xm*s
                    q = 1.0_dp-s
                else
                    q = fa/fc
                    r = fb/fc
                    p = s*(2.0_dp*xm*q*(q-r)-(b-a)*(r-1.0_dp))
                    q = (q-1.0_dp)*(r-1.0_dp)*(s-1.0_dp)
                end if
                if (p > 0.0_dp) then
                    q = -q
                end if
                p = abs(p)
                if (2.0_dp*p < min(3.0_dp*xm*q-abs(tol1*q),abs(e*q))) then
                    e = d
                    d = p/q
                else
                    d = xm
                    e = d
                end if
            else
                d = xm
                e = d
            end if
            a = b
            fa = fb
            b = b+merge(d,sign(tol1,xm),abs(d) > tol1)
            fb = func(b,xopt)
        end do
        call write_to_logfile("zbrent: exceeded maximum iterations")
        zbrent=b
    end function zbrent

    !---------------------------------------------------------------------------
    ! GEOP function
    !
    ! Description: Creates a geometric progression.  The code is adapted from
    !              Numerical Recipes, with changes made to fit the needs here.
    !
    ! Usage: geop(first,factor,n)
    !
    !           first   Starting value for the series
    !           factor  Multiplication factor
    !           n       Length of desired series
    !
    ! Result: The n-length sequence (first, first*factor,first*factor^2,...)
    !
    !---------------------------------------------------------------------------
    function geop(first,factor,n)
        real(kind=DP), intent(in) :: first, factor
        integer, intent(in) :: n
        real(kind=DP), dimension(n) :: geop
        integer :: k, k2
        real(kind=DP) :: temp
        if (n > 0) geop(1)=first
        do k=2,n
            geop(k) = geop(k-1)*factor
        end do
    end function geop

    !---------------------------------------------------------------------------
    ! IMINLOC function
    !
    ! Description: Returns the location in an array of the minimum.  The code is
    !              adapted from Numerical Recipes for use in numerical differentiation.
    !
    ! Usage: iminloc(arr)
    !
    !           arr     An array of DP reals
    !
    ! Result: An integer indicating the location in ARR of its minimum
    !
    !---------------------------------------------------------------------------
    function iminloc(arr)
        real(kind=DP), dimension(:), intent(in) :: arr
        integer, dimension(1) :: imin
        integer :: iminloc
        imin = minloc(arr(:))
        iminloc=imin(1)
    end function iminloc


end module rcrutil
