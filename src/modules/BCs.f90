        module boundaries
        use declarations
        use fluid_vars
        implicit none
        
        contains
        
!       calculate function of free parameters at each data point
        subroutine fcn(nfunc,nvar,alpha,fvec,iflag)        
        integer, intent(in) :: nfunc, nvar
        double precision, intent(in) :: alpha(nvar)
        integer, intent(out) :: iflag
        double precision, intent(out) :: fvec(nfunc)
        integer cur_pos
        
        fvec = 0

        do cur_pos = 1, nodes
         ! Boundary condition for horizontal walls
         if (node_array(cur_pos)%y == y_dim-1 .or. node_array(cur_pos)%y == 0 .or. (BC == channelCent .and. ((node_array(cur_pos)%y == floor(y_dim/3.) .or. node_array(cur_pos)%y == floor(2*y_dim/3.)) .and. node_array(cur_pos)%x >= floor(x_dim/3.) .and. node_array(cur_pos)%x <= floor(2*x_dim/3.))) .or. (BC == channelObs .and. (node_array(cur_pos)%y == floor(y_dim/2.) .and. node_array(cur_pos)%x >= floor(x_dim/3.) .and. node_array(cur_pos)%x <= floor(2*x_dim/3.)))) then
          fvec(cur_pos*3-2:cur_pos*3) = horizontalBC(cur_pos)
          ! Condition for top wall in Couette flow
          if (BC == couette .and. node_array(cur_pos)%y == y_dim-1) fvec(cur_pos*3-2) = fvec(cur_pos*3-2) - uwall
         ! Boundary conditions for vertical walls
         else if ((BC == channelCent .and. ((node_array(cur_pos)%x == floor(x_dim/3.) .or. node_array(cur_pos)%x == floor(2*x_dim/3.)) .and. node_array(cur_pos)%y > floor(y_dim/3.) .and. node_array(cur_pos)%y < floor(2*y_dim/3.))) .or. (BC == channelObs .and. ((node_array(cur_pos)%x == floor(x_dim/3.) .or. node_array(cur_pos)%x == floor(2*x_dim/3.)) .and. node_array(cur_pos)%y < floor(y_dim/2.)))) then
          fvec(cur_pos*3-2:cur_pos*3) = verticalBC(cur_pos)
         else if (BC == channelCircle .and. ((((node_array(cur_pos)%y-cent(2))**2 + (node_array(cur_pos)%x-cent(1))**2) >= r2) .and. (((node_array(cur_pos)%y+1-cent(2))**2 + (node_array(cur_pos)%x-cent(1))**2) < r2 .or. ((node_array(cur_pos)%y-cent(2))**2 + (node_array(cur_pos)%x+1-cent(1))**2) < r2 .or. ((node_array(cur_pos)%y-cent(2))**2 + (node_array(cur_pos)%x-1-cent(1))**2) < r2 .or. ((node_array(cur_pos)%y-1-cent(2))**2 + (node_array(cur_pos)%x-cent(1))**2) < r2))) then
         ! BC x=0
         else if (node_array(cur_pos)%x == 0) then
          fvec(cur_pos*3-2:cur_pos*3) = leftBC(cur_pos)
         ! BC x=right side
         else if (node_array(cur_pos)%x == x_dim-1) then
          fvec(cur_pos*3-2:cur_pos*3) = rightBC(cur_pos)
         ! If point is not on any boundary, apply interior conditions
         else
          fvec(cur_pos*3-2:cur_pos*3) = interiorConditions(cur_pos)
         end if
        end do

        return
        
        end subroutine
        
        ! Boundary conditions
        function horizontalBC(cur_pos)
        integer, intent(in) :: cur_pos
        double precision horizontalBC(3)
        
        horizontalBC = (/ u(cur_pos), v(cur_pos), (1./Re)*(d2vdx2(cur_pos) + d2vdy2(cur_pos)) - dpdy(cur_pos) /)
        
        end function
        
        function verticalBC(cur_pos)
        integer, intent(in) :: cur_pos
        double precision verticalBC(3)
        
        verticalBC = (/ u(cur_pos), v(cur_pos), (1./Re)*(d2udx2(cur_pos) + d2udy2(cur_pos)) - dpdx(cur_pos) /)
        
        end function
        
        function rightBC(cur_pos)
        integer, intent(in) :: cur_pos
        double precision rightBC(3)
        
        if (right_outflow) then
         rightBC = (/ dudx(cur_pos) + dvdy(cur_pos), dvdx(cur_pos), (1./Re)*(d2udx2(cur_pos) + d2udy2(cur_pos)) - u(cur_pos)*dudx(cur_pos) - v(cur_pos)*dudy(cur_pos) - dpdx(cur_pos) /)
        else if (right_set) then
         rightBC = (/ u(cur_pos) - uright, v(cur_pos) - vright, p(cur_pos) - pright /)
        end if
        
        end function
        
        function leftBC(cur_pos)
        integer, intent(in) :: cur_pos
        double precision leftBC(3)
        
        leftBC = (/ u(cur_pos) - uleft, v(cur_pos) - vleft, p(cur_pos) - pleft /)
        
        end function
        
        function interiorConditions(cur_pos)
        integer, intent(in) :: cur_pos
        double precision interiorConditions(3)
        
        interiorConditions = (/ dudx(cur_pos) + dvdy(cur_pos), u(cur_pos)*dudx(cur_pos) + v(cur_pos)*dudy(cur_pos) + dpdx(cur_pos) - (1./Re)*(d2udx2(cur_pos) + d2udy2(cur_pos)), u(cur_pos)*dvdx(cur_pos) + v(cur_pos)*dvdy(cur_pos) + dpdy(cur_pos) - (1./Re)*(d2vdx2(cur_pos) + d2vdy2(cur_pos)) /)
        
        end function
        
        end module
