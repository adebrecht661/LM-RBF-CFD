        module output
        use grid
        use declarations
        use fluid_vars
        implicit none
        
        private LMCode
        
        contains

        subroutine ResultToTxt()
        
        integer stat,stat1,cur_pos

!       if solution is successful, print result to file for processing
        if (info  == 1 .or. info == 2 .or. info == 3) then
         open(unit=40,file='out/result.txt',iostat=stat,status='replace')
         open(unit=50,file='out/info.txt',iostat=stat1,status='replace')
         call LMCode
         if (stat /= 0) then
          print *, 'The file could not be opened for writing.'
         else if (stat1 /= 0) then
          print *, 'Output information could not be written to file. Continuing.'
          write (40,*) 'x y u v p'
          do cur_pos = 1, nodes
           write (40,*) dble(node_array(cur_pos)%x)/res, dble(node_array(cur_pos)%y)/res, u(cur_pos), v(cur_pos), p(cur_pos)
          end do
         else
          ! Print out resolution and boundary conditions
          select case (BC)
           case (channel)
            write (50,*) 'The resolution of this run was ', res, ' points/CU. Channel flow boundary conditions were used.'
            write (50,*) 'The Reynolds number of this run was ', Re, '.'
           case (couette)
            write (50,*) 'The resolution of this run was ', res, ' points/CU. Couette flow boundary conditions were used.'
            write (50,*) 'The Reynolds number of this run was ', Re, '.'
           case (channelCent)
            write (50,*) 'The resolution of this run was ', res, ' points/CU. Channel flow boundary conditions with a central obstruction were used.' 
            write (50,*) 'The Reynolds number of this run was ', Re, '.'
           case (channelObs)
            write (50,*) 'The resolution of this run was ', res, ' points/CU. Channel flow boundary conditions with an obstruction were used.' 
            write (50,*) 'The Reynolds number of this run was ', Re, '.'
           case (channelCircle)
            write (50,*) 'The resolution of this run was ', res, ' points/CU. Channel flow boundary conditions with a circular obstruction were used.' 
            write (50,*) 'The Reynolds number of this run was ', Re, '.'
          end select
          write (40,*) 'x y u v p'
          do cur_pos = 1, nodes
           write (40,*) dble(node_array(cur_pos)%x)/res, dble(node_array(cur_pos)%y)/res, u(cur_pos), v(cur_pos), p(cur_pos)
          end do
         end if
        else
         open(unit=50,file='out/info.txt',iostat=stat1,status='replace')
         call LMCode
         print *, 'An error occurred. There may be more information in info.txt.'
        end if
        
        end subroutine
        
        subroutine LMCode()
        
!       Print conditions on result or error if unsuccessful (write to info.txt?)
        select case (info)
         case (0)
          write (50,*) 'Improper input parameters to lmdif.'
         case (1)
          write (50,*) 'The relative error in the sum of squares is at most ', tol, '.'
         case (2)
          write (50,*) 'The relative error between the estimated solution and the solution is at most ', tol, '.'
         case (3)
          write(50,*) 'The relative errors in the sum of squares and between the estimated and actual solution are at most ', tol, '.'
         case (4)
          write(50,*) 'Function is orthogonal to the columns of the Jacobian.'
         case (5)
          write(50,*) 'Subroutine calls have reached or exceeded ', 200*(nvar+1), '.'
         case (6)
          write(50,*) 'Tolerance is too small. No further reduction in the sum of squares is possible.'
         case (7)
          write(50,*) 'Tolerance is too small. No further improvement in the approximate solution is possible.'
        end select
        
        return
        
        end subroutine
        
        end module
