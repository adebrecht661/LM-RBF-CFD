!       Solves for the best-fit parameters of cubic equation using MINPACK 
!       routines lmdif, lmder, and lmstr (using provided driver routines)
        
        program lm_example
!       use lmdif1 as module for LM solver
        use LM
        implicit none
        
!       declare variables        
        integer i,seed(12),sleep_time_int
        integer a_actual,b_actual,c_actual
        integer m,n,info,lwa
        double precision sleep_time, tol
        logical ex
!       declare matrices
        integer, allocatable :: iwa(:),ipvt(:)
        double precision, allocatable :: x(:),y(:),error(:),soln(:),fvec(:)
        double precision, allocatable :: soln_guess(:),wa(:),fjac(:,:)  
!       initialize variables
        seed = 0
        a_actual = 3            ! cubic coefficient
        b_actual = 10           ! quadratic coefficient
        c_actual = -5           ! linear coefficient
        n = 3
        tol = 1e-10             ! tolerance for solver

!       Ask user for # of points        
        print *, 'How many data points would you like to generate?'
        read(*,*) m

!       allocate arrays according to number of points m (# of equations) and free parameters n        
        lwa = m*n+5*n+m
        allocate(x(m),y(m),error(m),fvec(m))
        allocate(soln_guess(n),iwa(n),soln(n),wa(lwa))
        allocate(fjac(m,n),ipvt(n))

!       Ask user for starting guess for parameters        
        print *, 'Input your initial guess for the free parameters a,b,c:'
        read(*,*) soln_guess
        
!       Generate some random points and random error

        inquire(file='/dev/urandom', exist=ex)          ! determine if entropy file exists
        if (ex) then                                    ! if entropy file exists, use it to seed RNG
           open(unit=20, file='/dev/urandom', access="stream", &
            form="UNFORMATTED", status="old", action="read")
          read(20) seed
        else                                            ! if entropy file does not exist, seed RNG from system clock
          print *, '/dev/urandom not present. Generating random seed from &
           system clock. This may take up to 60 seconds.'
          do i=1,12
            seed(i) = time()
            call random_seed(put=seed)          ! reseed with more-full seed each time
            call random_number(sleep_time)
            sleep_time_int = nint(sleep_time*5)
            call sleep(sleep_time_int)          ! sleep for random amount between each seed entry
          end do
          print *, 'Random seed generated.'
        end if

!       fill random seed, then randomly fill x values from 0 to 100 (mostly just because - can change to smaller or larger with little effect)
        call random_seed(put=seed)
        call random_number(x)
        x = x*100
                
!       calculate some data points and add random error to create data set for best fit
        do i=1,m
          y(i) = target(x(i))
        end do
        call random_number(error)
        error = error*0.1 - 0.05
        y = y + error
        
!       solve with system f(x) = 0 only
        soln = soln_guess
        call lmdif1(fcn,m,n,soln,fvec,tol,info,iwa,wa,lwa)

!       Print conditions on result or error if unsuccessful        
        print *, 'lmdif result:'
        select case (info)
          case (0)
            print *, 'Improper input parameters.'
          case (1)
            print *, 'The relative error in the sum of squares is at &
             &most ', tol, '.'
          case (2)
            print *, 'The relative error between the estimated &
             &solution and the solution is at most ', tol, '.'
          case (3)
            print *, 'The relative errors in the sum of squares and &
             &between the estimated and actual solution are at most ', tol, '.'
          case (4)
            print *, 'Function is orthogonal to the columns of the Jacobian.'
          case (5)
            print *, 'Subroutine calls have reached or exceeded ', &
             200*(n+1), '.'
          case (6)
            print *, 'Tolerance is too small. No further reduction in &
             &the sum of squares is possible.'
          case (7)
            print *, 'Tolerance is too small. No further improvement in &
             &the approximate solution is possible.'
        end select   

!       if solution is successful, print result                    
        if (info == 1 .or. info == 2 .or. info == 3) then
          print *, 'a = ', soln(1)
          print *, 'b = ', soln(2)
          print *, 'c = ', soln(3)
        end if
        
!       change length of work array
        deallocate(wa)
        lwa = 5*n+m
        allocate(wa(lwa))
        
!       solve with functions and jacobian
        soln = soln_guess       
        call lmder1(jac_fcn,m,n,soln,fvec,fjac,m,tol,info,ipvt,wa,lwa)
        
!       Print conditions on result or error if unsuccessful        
        print *, 'lmder result:'
        select case (info)
          case (0)
            print *, 'Improper input parameters.'
          case (1)
            print *, 'The relative error in the sum of squares is at &
             &most ', tol, '.'
          case (2)
            print *, 'The relative error between the estimated &
             &solution and the solution is at most ', tol, '.'
          case (3)
            print *, 'The relative errors in the sum of squares and &
             &between the estimated and actual solution are at most ', tol, '.'
          case (4)
            print *, 'Function is orthogonal to the columns of the Jacobian.'
          case (5)
            print *, 'Subroutine calls have reached or exceeded ', &
             100*(n+1), '.'
          case (6)
            print *, 'Tolerance is too small. No further reduction in &
             &the sum of squares is possible.'
          case (7)
            print *, 'Tolerance is too small. No further improvement in &
             &the approximate solution is possible.'
        end select   

!       if solution is successful, print result                    
        if (info == 1 .or. info == 2 .or. info == 3) then
          print *, 'a = ', soln(1)
          print *, 'b = ', soln(2)
          print *, 'c = ', soln(3)
        end if

!       change size of fjac
        deallocate(fjac)
        allocate(fjac(n,n))

!       solve with functions and jacobian in less space
        soln = soln_guess
        call lmstr1(jac_row_fcn,m,n,soln,fvec,fjac,n,tol,info,ipvt,wa,lwa)
        
!       Print conditions on result or error if unsuccessful        
        print *, 'lmstr result:'
        select case (info)
          case (0)
            print *, 'Improper input parameters.'
          case (1)
            print *, 'The relative error in the sum of squares is at &
             &most ', tol, '.'
          case (2)
            print *, 'The relative error between the estimated &
             &solution and the solution is at most ', tol, '.'
          case (3)
            print *, 'The relative errors in the sum of squares and &
             &between the estimated and actual solution are at most ', tol, '.'
          case (4)
            print *, 'Function is orthogonal to the columns of the Jacobian.'
          case (5)
            print *, 'Subroutine calls have reached or exceeded ', &
             100*(n+1), '.'
          case (6)
            print *, 'Tolerance is too small. No further reduction in &
             &the sum of squares is possible.'
          case (7)
            print *, 'Tolerance is too small. No further improvement in &
             &the approximate solution is possible.'
        end select   

!       if solution is successful, print result                    
        if (info == 1 .or. info == 2 .or. info == 3) then
          print *, 'a = ', soln(1)
          print *, 'b = ', soln(2)
          print *, 'c = ', soln(3)
        end if
                  
        stop
        
        contains

!       calculate the actual values of the starting function at a point        
        double precision function target(x)
        double precision x
        
        target = a_actual*x**3 + b_actual*x**2 + c_actual*x
        
        return
        
        end function

!       calculate function of free parameters at each data point        
        subroutine fcn(m,n,soln,fvec,iflag)        
        integer m,n,iflag
        double precision soln(n),fvec(m)
        
        fvec = y - soln(1)*x**3 - soln(2)*x**2 - soln(3)*x
        
        return
        
        end subroutine

!       calculate function of free parameters and jacobian at each data point        
        subroutine jac_fcn(m,n,soln,fvec,fjac,ldfjac,iflag)
        integer m,n,ldfjac,iflag
        double precision soln(n),fvec(m),fjac(ldfjac,n)

        if (iflag == 1) then
          fvec = y - soln(1)*x**3 - soln(2)*x**2 - soln(3)*x
        else if (iflag == 2) then
          fjac(1:m,1) = -x**3
          fjac(1:m,2) = -x**2
          fjac(1:m,3) = -x
        end if

        return
        
        end
     
!       calculate function of free parameters at each data point and row of jacobian at one data point        
        subroutine jac_row_fcn(m,n,soln,fvec,fjrow,iflag)
        integer m,n,iflag
        double precision soln(n),fvec(m),fjrow(n)

        if (iflag == 1) then
          fvec = y - soln(1)*x**3 - soln(2)*x**2 - soln(3)*x
        else if (iflag /= 1 .and. iflag > 0) then
          fjrow(1) = -x(iflag-1)**3
          fjrow(2) = -x(iflag-1)**2
          fjrow(3) = -x(iflag-1)
        end if

        return
        
        end
        
        end program
