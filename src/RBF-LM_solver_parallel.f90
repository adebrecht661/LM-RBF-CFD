        program fluid_solver
        use declarations
        use LM
        use grid
        use Gaussian_RBF
!        use TPS_RBF            !!!! Not working yet
        use input
        use output
!	use system
        implicit none

!       get grid size and resolution for input to lmdif
        call GetProblemInfo
        call CalcShapeParam

!       Total number of nodes in the system
        nodes = y_dim*x_dim
        call SetUpGrid
        nvar = 3*nodes          ! Three variables (u,v,p) at each point
        nfunc = nvar            ! nfunc >= nvar, need at least three functions per point

!       Calculate length of working array and allocate matrices
        lwa = nfunc*nvar+5*nvar+nfunc
        allocate(alpha(nvar),iwa(nvar),fvec(nfunc),wa(lwa))

        call GetAlpha(alpha)

!       solve with system f(x) = 0 only
        call lmdif1(fcn,nfunc,nvar,alpha,fvec,tol,info,iwa,wa,lwa)

        print *, 'Completed calculation. Writing out.'
        call ResultToTxt

        deallocate(alpha,iwa,fvec,wa)

        stop
        
        contains
        
!       calculate function of free parameters at each data point
        subroutine fcn(nfunc,nvar,alpha,fvec,iflag)        
        integer, intent(in) :: nfunc, nvar
        double precision, intent(inout) :: alpha(nvar),fvec(nfunc)
        integer, intent(out) :: iflag
        double precision c
        integer cur_pos
        
        fvec = 0

!!!!!!! Less expensive now - but it may be worth making the RBF functions into array accesses instead at some point
        select case (BC)
         case (1)               ! Channel flow
          !$OMP parallel do
          do cur_pos = 1, nodes
           ! Boundary condition for top and bottom walls
           if (node_array(cur_pos)%y == y_dim-1 .or. node_array(cur_pos)%y == 0) then
            fvec(cur_pos*3-2) = sum(alpha(:nodes)*psi(cur_pos,node_array,nodes,c))
            fvec(cur_pos*3-1) = sum(alpha(nodes:nodes*2)*psi(cur_pos,node_array,nodes,c))
            fvec(cur_pos*3) = (-1.*sum(alpha(nodes*2:nodes*3)*dy(cur_pos,node_array,nodes,c))) + (1./Re)*(sum(alpha(nodes:nodes*2)*d2x(cur_pos,node_array,nodes,c)) + sum(alpha(nodes:nodes*2)*d2y(cur_pos,node_array,nodes,c)))
           ! BC x=0 - Above conditions satisfy need for .and. .not. (y=0 .or. y=1), because this will be skipped 
           else if (node_array(cur_pos)%x == 0) then
            fvec(cur_pos*3-2) = sum(alpha(:nodes)*psi(cur_pos,node_array,nodes,c)) - uleft
            fvec(cur_pos*3-1) = sum(alpha(nodes:nodes*2)*psi(cur_pos,node_array,nodes,c)) - vleft
            fvec(cur_pos*3) = sum(alpha(nodes*2:nodes*3)*psi(cur_pos,node_array,nodes,c)) - pleft
           ! BC x=right side
           else if (node_array(cur_pos)%x == x_dim-1) then
            fvec(cur_pos*3-2) = sum(alpha(:nodes)*dx(cur_pos,node_array,nodes,c)) + sum(alpha(nodes:nodes*2)*dy(cur_pos,node_array,nodes,c))
            fvec(cur_pos*3-1) = sum(alpha(nodes:nodes*2)*dx(cur_pos,node_array,nodes,c))
            fvec(cur_pos*3) = sum(-alpha(nodes*2:nodes*3)*dx(cur_pos,node_array,nodes,c)) + (1./Re)*(sum(alpha(:nodes)*d2x(cur_pos,node_array,nodes,c)) + sum(alpha(:nodes)*d2y(cur_pos,node_array,nodes,c))) - sum(alpha(:nodes)*psi(cur_pos,node_array,nodes,c))*sum(alpha(:nodes)*dx(cur_pos,node_array,nodes,c)) - sum(alpha(nodes:nodes*2)*psi(cur_pos,node_array,nodes,c))*sum(alpha(:nodes)*dy(cur_pos,node_array,nodes,c))
           ! If point is not on any boundary, apply interior conditions
           else
            fvec(cur_pos*3-2) = sum(alpha(:nodes)*dx(cur_pos,node_array,nodes,c)) + sum(alpha(nodes:nodes*2)*dy(cur_pos,node_array,nodes,c))
            fvec(cur_pos*3-1) = sum(alpha(:nodes)*psi(cur_pos,node_array,nodes,c))*sum(alpha(:nodes)*dx(cur_pos,node_array,nodes,c)) + sum(alpha(nodes:nodes*2)*psi(cur_pos,node_array,nodes,c))*sum(alpha(:nodes)*dy(cur_pos,node_array,nodes,c)) + sum(alpha(nodes*2:nodes*3)*dx(cur_pos,node_array,nodes,c)) - (1./Re)*(sum(alpha(:nodes)*d2x(cur_pos,node_array,nodes,c)) + sum(alpha(:nodes)*d2y(cur_pos,node_array,nodes,c)))
            fvec(cur_pos*3) = sum(alpha(:nodes)*psi(cur_pos,node_array,nodes,c))*sum(alpha(nodes:nodes*2)*dx(cur_pos,node_array,nodes,c)) + sum(alpha(nodes:nodes*2)*psi(cur_pos,node_array,nodes,c))*sum(alpha(nodes:nodes*2)*dy(cur_pos,node_array,nodes,c)) + sum(alpha(nodes*2:nodes*3)*dy(cur_pos,node_array,nodes,c)) - (1./Re)*(sum(alpha(nodes:nodes*2)*d2x(cur_pos,node_array,nodes,c)) + sum(alpha(nodes:nodes*2)*d2y(cur_pos,node_array,nodes,c)))
           end if
          end do
          !$OMP end parallel do
         case (2)
          !$OMP parallel do
          do cur_pos = 1, nodes
           ! Boundary condition for top wall
           if (node_array(cur_pos)%y == y_dim-1) then
            fvec(cur_pos*3-2) = sum(alpha(:nodes)*psi(cur_pos,node_array,nodes,c)) - uwall
            fvec(cur_pos*3-1) = sum(alpha(nodes:nodes*2)*psi(cur_pos,node_array,nodes,c))
            fvec(cur_pos*3) = sum(alpha(nodes*2:nodes*3)*psi(cur_pos,node_array,nodes,c))
           ! Boundary condition for bottom wall
           else if (node_array(cur_pos)%y == 0) then
            fvec(cur_pos*3-2) = sum(alpha(:nodes)*psi(cur_pos,node_array,nodes,c))
            fvec(cur_pos*3-1) = sum(alpha(nodes:nodes*2)*psi(cur_pos,node_array,nodes,c))
            fvec(cur_pos*3) = sum(alpha(nodes*2:nodes*3)*psi(cur_pos,node_array,nodes,c))
           ! If point is not on any boundary, apply interior conditions
           else
            fvec(cur_pos*3-2) = sum(alpha(:nodes)*dx(cur_pos,node_array,nodes,c)) + sum(alpha(nodes:nodes*2)*dy(cur_pos,node_array,nodes,c))
            fvec(cur_pos*3-1) = sum(alpha(:nodes)*psi(cur_pos,node_array,nodes,c))*sum(alpha(:nodes)*dx(cur_pos,node_array,nodes,c)) + sum(alpha(nodes:nodes*2)*psi(cur_pos,node_array,nodes,c))*sum(alpha(:nodes)*dy(cur_pos,node_array,nodes,c)) + sum(alpha(nodes*2:nodes*3)*dx(cur_pos,node_array,nodes,c)) - (1./Re)*(sum(alpha(:nodes)*d2x(cur_pos,node_array,nodes,c)) + sum(alpha(:nodes)*d2y(cur_pos,node_array,nodes,c)))
            fvec(cur_pos*3) = sum(alpha(:nodes)*psi(cur_pos,node_array,nodes,c))*sum(alpha(nodes:nodes*2)*dx(cur_pos,node_array,nodes,c)) + sum(alpha(nodes:nodes*2)*psi(cur_pos,node_array,nodes,c))*sum(alpha(nodes:nodes*2)*dy(cur_pos,node_array,nodes,c)) + sum(alpha(nodes*2:nodes*3)*dy(cur_pos,node_array,nodes,c)) - (1./Re)*(sum(alpha(nodes:nodes*2)*d2x(cur_pos,node_array,nodes,c)) + sum(alpha(nodes:nodes*2)*d2y(cur_pos,node_array,nodes,c)))
           end if
          end do
          !$OMP end parallel do
         case (3)
          !$OMP parallel do
          do cur_pos = 1, nodes
           ! Boundary condition for horizontal walls
           if (node_array(cur_pos)%y == y_dim-1 .or. node_array(cur_pos)%y == 0 .or. ((node_array(cur_pos)%y == floor(y_dim/3.) .or. node_array(cur_pos)%y == floor(2*y_dim/3.)) .and. node_array(cur_pos)%x >= floor(x_dim/3.) .and. node_array(cur_pos)%x <= floor(2*x_dim/3.))) then
            fvec(cur_pos*3-2) = sum(alpha(:nodes)*psi(cur_pos,node_array,nodes,c))
            fvec(cur_pos*3-1) = sum(alpha(nodes:nodes*2)*psi(cur_pos,node_array,nodes,c))
            fvec(cur_pos*3) = (-1.*sum(alpha(nodes*2:nodes*3)*dy(cur_pos,node_array,nodes,c))) + (1./Re)*(sum(alpha(nodes:nodes*2)*d2x(cur_pos,node_array,nodes,c)) + sum(alpha(nodes:nodes*2)*d2y(cur_pos,node_array,nodes,c)))
           ! Boundary condition for vertical walls
           else if ((node_array(cur_pos)%x == floor(x_dim/3.) .or. node_array(cur_pos)%x == floor(2*x_dim/3.)) .and. node_array(cur_pos)%y > floor(y_dim/3.) .and. node_array(cur_pos)%y < floor(2*y_dim/3.)) then
            fvec(cur_pos*3-2) = sum(alpha(:nodes)*psi(cur_pos,node_array,nodes,c))
            fvec(cur_pos*3-1) = sum(alpha(nodes:nodes*2)*psi(cur_pos,node_array,nodes,c))
            fvec(cur_pos*3) = (-1.*sum(alpha(nodes*2:nodes*3)*dx(cur_pos,node_array,nodes,c))) + (1./Re)*(sum(alpha(:nodes)*d2x(cur_pos,node_array,nodes,c)) + sum(alpha(:nodes)*d2y(cur_pos,node_array,nodes,c)))
           ! BC x=0 - Above conditions satisfy need for .and. .not. (y=0 .or. y=1), because this will be skipped 
           else if (node_array(cur_pos)%x == 0) then
            fvec(cur_pos*3-2) = sum(alpha(:nodes)*psi(cur_pos,node_array,nodes,c)) - uleft
            fvec(cur_pos*3-1) = sum(alpha(nodes:nodes*2)*psi(cur_pos,node_array,nodes,c)) - vleft
            fvec(cur_pos*3) = sum(alpha(nodes*2:nodes*3)*psi(cur_pos,node_array,nodes,c)) - pleft
           ! BC x=right side
           else if (node_array(cur_pos)%x == x_dim-1) then
            fvec(cur_pos*3-2) = sum(alpha(:nodes)*dx(cur_pos,node_array,nodes,c)) + sum(alpha(nodes:nodes*2)*dy(cur_pos,node_array,nodes,c))
            fvec(cur_pos*3-1) = sum(alpha(nodes:nodes*2)*dx(cur_pos,node_array,nodes,c))
            fvec(cur_pos*3) = sum(-alpha(nodes*2:nodes*3)*dx(cur_pos,node_array,nodes,c)) + (1./Re)*(sum(alpha(:nodes)*d2x(cur_pos,node_array,nodes,c)) + sum(alpha(:nodes)*d2y(cur_pos,node_array,nodes,c))) - sum(alpha(:nodes)*psi(cur_pos,node_array,nodes,c))*sum(alpha(:nodes)*dx(cur_pos,node_array,nodes,c)) - sum(alpha(nodes:nodes*2)*psi(cur_pos,node_array,nodes,c))*sum(alpha(:nodes)*dy(cur_pos,node_array,nodes,c))
           ! If point is not on any boundary, apply interior conditions
           else
            fvec(cur_pos*3-2) = sum(alpha(:nodes)*dx(cur_pos,node_array,nodes,c)) + sum(alpha(nodes:nodes*2)*dy(cur_pos,node_array,nodes,c))
            fvec(cur_pos*3-1) = sum(alpha(:nodes)*psi(cur_pos,node_array,nodes,c))*sum(alpha(:nodes)*dx(cur_pos,node_array,nodes,c)) + sum(alpha(nodes:nodes*2)*psi(cur_pos,node_array,nodes,c))*sum(alpha(:nodes)*dy(cur_pos,node_array,nodes,c)) + sum(alpha(nodes*2:nodes*3)*dx(cur_pos,node_array,nodes,c)) - (1./Re)*(sum(alpha(:nodes)*d2x(cur_pos,node_array,nodes,c)) + sum(alpha(:nodes)*d2y(cur_pos,node_array,nodes,c)))
            fvec(cur_pos*3) = sum(alpha(:nodes)*psi(cur_pos,node_array,nodes,c))*sum(alpha(nodes:nodes*2)*dx(cur_pos,node_array,nodes,c)) + sum(alpha(nodes:nodes*2)*psi(cur_pos,node_array,nodes,c))*sum(alpha(nodes:nodes*2)*dy(cur_pos,node_array,nodes,c)) + sum(alpha(nodes*2:nodes*3)*dy(cur_pos,node_array,nodes,c)) - (1./Re)*(sum(alpha(nodes:nodes*2)*d2x(cur_pos,node_array,nodes,c)) + sum(alpha(nodes:nodes*2)*d2y(cur_pos,node_array,nodes,c)))
           end if
          end do
         case (4)
          !$OMP parallel do
          do cur_pos = 1, nodes
           ! Boundary condition for horizontal walls
           if (node_array(cur_pos)%y == y_dim-1 .or. node_array(cur_pos)%y == 0 .or. (node_array(cur_pos)%y == floor(y_dim/2.) .and. node_array(cur_pos)%x >= floor(x_dim/3.) .and. node_array(cur_pos)%x <= floor(2*x_dim/3.))) then
            fvec(cur_pos*3-2) = sum(alpha(:nodes)*psi(cur_pos,node_array,nodes,c))
            fvec(cur_pos*3-1) = sum(alpha(nodes:nodes*2)*psi(cur_pos,node_array,nodes,c))
            fvec(cur_pos*3) = (-1.*sum(alpha(nodes*2:nodes*3)*dy(cur_pos,node_array,nodes,c))) + (1./Re)*(sum(alpha(nodes:nodes*2)*d2x(cur_pos,node_array,nodes,c)) + sum(alpha(nodes:nodes*2)*d2y(cur_pos,node_array,nodes,c)))
           ! Boundary condition for vertical walls
           else if ((node_array(cur_pos)%x == floor(x_dim/3.) .or. node_array(cur_pos)%x == floor(2*x_dim/3.)) .and. node_array(cur_pos)%y < floor(y_dim/2.)) then
            fvec(cur_pos*3-2) = sum(alpha(:nodes)*psi(cur_pos,node_array,nodes,c))
            fvec(cur_pos*3-1) = sum(alpha(nodes:nodes*2)*psi(cur_pos,node_array,nodes,c))
            fvec(cur_pos*3) = (-1.*sum(alpha(nodes*2:nodes*3)*dx(cur_pos,node_array,nodes,c))) + (1./Re)*(sum(alpha(:nodes)*d2x(cur_pos,node_array,nodes,c)) + sum(alpha(:nodes)*d2y(cur_pos,node_array,nodes,c)))
           ! BC x=0 - Above conditions satisfy need for .and. .not. (y=0 .or. y=1), because this will be skipped 
           else if (node_array(cur_pos)%x == 0) then
            fvec(cur_pos*3-2) = sum(alpha(:nodes)*psi(cur_pos,node_array,nodes,c)) - uleft
            fvec(cur_pos*3-1) = sum(alpha(nodes:nodes*2)*psi(cur_pos,node_array,nodes,c)) - vleft
            fvec(cur_pos*3) = sum(alpha(nodes*2:nodes*3)*psi(cur_pos,node_array,nodes,c)) - pleft
           ! BC x=right side
           else if (node_array(cur_pos)%x == x_dim-1) then
            fvec(cur_pos*3-2) = sum(alpha(:nodes)*dx(cur_pos,node_array,nodes,c)) + sum(alpha(nodes:nodes*2)*dy(cur_pos,node_array,nodes,c))
            fvec(cur_pos*3-1) = sum(alpha(nodes:nodes*2)*dx(cur_pos,node_array,nodes,c))
            fvec(cur_pos*3) = sum(-alpha(nodes*2:nodes*3)*dx(cur_pos,node_array,nodes,c)) + (1./Re)*(sum(alpha(:nodes)*d2x(cur_pos,node_array,nodes,c)) + sum(alpha(:nodes)*d2y(cur_pos,node_array,nodes,c))) - sum(alpha(:nodes)*psi(cur_pos,node_array,nodes,c))*sum(alpha(:nodes)*dx(cur_pos,node_array,nodes,c)) - sum(alpha(nodes:nodes*2)*psi(cur_pos,node_array,nodes,c))*sum(alpha(:nodes)*dy(cur_pos,node_array,nodes,c))
           ! If point is not on any boundary, apply interior conditions
           else
            fvec(cur_pos*3-2) = sum(alpha(:nodes)*dx(cur_pos,node_array,nodes,c)) + sum(alpha(nodes:nodes*2)*dy(cur_pos,node_array,nodes,c))
            fvec(cur_pos*3-1) = sum(alpha(:nodes)*psi(cur_pos,node_array,nodes,c))*sum(alpha(:nodes)*dx(cur_pos,node_array,nodes,c)) + sum(alpha(nodes:nodes*2)*psi(cur_pos,node_array,nodes,c))*sum(alpha(:nodes)*dy(cur_pos,node_array,nodes,c)) + sum(alpha(nodes*2:nodes*3)*dx(cur_pos,node_array,nodes,c)) - (1./Re)*(sum(alpha(:nodes)*d2x(cur_pos,node_array,nodes,c)) + sum(alpha(:nodes)*d2y(cur_pos,node_array,nodes,c)))
            fvec(cur_pos*3) = sum(alpha(:nodes)*psi(cur_pos,node_array,nodes,c))*sum(alpha(nodes:nodes*2)*dx(cur_pos,node_array,nodes,c)) + sum(alpha(nodes:nodes*2)*psi(cur_pos,node_array,nodes,c))*sum(alpha(nodes:nodes*2)*dy(cur_pos,node_array,nodes,c)) + sum(alpha(nodes*2:nodes*3)*dy(cur_pos,node_array,nodes,c)) - (1./Re)*(sum(alpha(nodes:nodes*2)*d2x(cur_pos,node_array,nodes,c)) + sum(alpha(nodes:nodes*2)*d2y(cur_pos,node_array,nodes,c)))
           end if
          end do
         case (5)
          !$OMP parallel do
          do cur_pos = 1, nodes
           ! Boundary condition for horizontal walls
           if (node_array(cur_pos)%y == y_dim-1 .or. node_array(cur_pos)%y == 0) then
            fvec(cur_pos*3-2) = sum(alpha(:nodes)*psi(cur_pos,node_array,nodes,c))
            fvec(cur_pos*3-1) = sum(alpha(nodes:nodes*2)*psi(cur_pos,node_array,nodes,c))
            fvec(cur_pos*3) = (-1.*sum(alpha(nodes*2:nodes*3)*dy(cur_pos,node_array,nodes,c))) + (1./Re)*(sum(alpha(nodes:nodes*2)*d2x(cur_pos,node_array,nodes,c)) + sum(alpha(nodes:nodes*2)*d2y(cur_pos,node_array,nodes,c)))
           ! BC for circle
           else if ((((node_array(cur_pos)%y-cent(2))**2 + (node_array(cur_pos)%x-cent(1))**2) >= r2) .and. (((node_array(cur_pos)%y+1-cent(2))**2 + (node_array(cur_pos)%x-cent(1))**2) < r2 .or. ((node_array(cur_pos)%y-cent(2))**2 + (node_array(cur_pos)%x+1-cent(1))**2) < r2 .or. ((node_array(cur_pos)%y-cent(2))**2 + (node_array(cur_pos)%x-1-cent(1))**2) < r2 .or. ((node_array(cur_pos)%y-1-cent(2))**2 + (node_array(cur_pos)%x-cent(1))**2) < r2)) then
            fvec(cur_pos*3-2) = sum(alpha(:nodes)*psi(cur_pos,node_array,nodes,c))
            fvec(cur_pos*3-1) = sum(alpha(nodes:nodes*2)*psi(cur_pos,node_array,nodes,c))
            !!!!!! Definitely not right - but need to conserve momentum in all directions (or perhaps just r [normal] direction??)
            fvec(cur_pos*3) = (-1.*sum(alpha(nodes*2:nodes*3)*dx(cur_pos,node_array,nodes,c)) + (1./Re)*(sum(alpha(:nodes)*d2x(cur_pos,node_array,nodes,c)) + sum(alpha(:nodes)*d2y(cur_pos,node_array,nodes,c))))**2 + (-1.*sum(alpha(nodes*2:nodes*3)*dy(cur_pos,node_array,nodes,c))) + (1./Re)*(sum(alpha(nodes:nodes*2)*d2x(cur_pos,node_array,nodes,c)) + sum(alpha(nodes:nodes*2)*d2y(cur_pos,node_array,nodes,c)))**2
           ! BC x=0 - Above conditions satisfy need for .and. .not. (y=0 .or. y=1), because this will be skipped 
           else if (node_array(cur_pos)%x == 0) then
            fvec(cur_pos*3-2) = sum(alpha(:nodes)*psi(cur_pos,node_array,nodes,c)) - uleft
            fvec(cur_pos*3-1) = sum(alpha(nodes:nodes*2)*psi(cur_pos,node_array,nodes,c)) - vleft
            fvec(cur_pos*3) = sum(alpha(nodes*2:nodes*3)*psi(cur_pos,node_array,nodes,c)) - pleft
           ! BC x=right side
           else if (node_array(cur_pos)%x == x_dim-1) then
            fvec(cur_pos*3-2) = sum(alpha(:nodes)*dx(cur_pos,node_array,nodes,c)) + sum(alpha(nodes:nodes*2)*dy(cur_pos,node_array,nodes,c))
            fvec(cur_pos*3-1) = sum(alpha(nodes:nodes*2)*dx(cur_pos,node_array,nodes,c))
            fvec(cur_pos*3) = sum(-alpha(nodes*2:nodes*3)*dx(cur_pos,node_array,nodes,c)) + (1./Re)*(sum(alpha(:nodes)*d2x(cur_pos,node_array,nodes,c)) + sum(alpha(:nodes)*d2y(cur_pos,node_array,nodes,c))) - sum(alpha(:nodes)*psi(cur_pos,node_array,nodes,c))*sum(alpha(:nodes)*dx(cur_pos,node_array,nodes,c)) - sum(alpha(nodes:nodes*2)*psi(cur_pos,node_array,nodes,c))*sum(alpha(:nodes)*dy(cur_pos,node_array,nodes,c))
           ! If point is not on any boundary, apply interior conditions
           else
            fvec(cur_pos*3-2) = sum(alpha(:nodes)*dx(cur_pos,node_array,nodes,c)) + sum(alpha(nodes:nodes*2)*dy(cur_pos,node_array,nodes,c))
            fvec(cur_pos*3-1) = sum(alpha(:nodes)*psi(cur_pos,node_array,nodes,c))*sum(alpha(:nodes)*dx(cur_pos,node_array,nodes,c)) + sum(alpha(nodes:nodes*2)*psi(cur_pos,node_array,nodes,c))*sum(alpha(:nodes)*dy(cur_pos,node_array,nodes,c)) + sum(alpha(nodes*2:nodes*3)*dx(cur_pos,node_array,nodes,c)) - (1./Re)*(sum(alpha(:nodes)*d2x(cur_pos,node_array,nodes,c)) + sum(alpha(:nodes)*d2y(cur_pos,node_array,nodes,c)))
            fvec(cur_pos*3) = sum(alpha(:nodes)*psi(cur_pos,node_array,nodes,c))*sum(alpha(nodes:nodes*2)*dx(cur_pos,node_array,nodes,c)) + sum(alpha(nodes:nodes*2)*psi(cur_pos,node_array,nodes,c))*sum(alpha(nodes:nodes*2)*dy(cur_pos,node_array,nodes,c)) + sum(alpha(nodes*2:nodes*3)*dy(cur_pos,node_array,nodes,c)) - (1./Re)*(sum(alpha(nodes:nodes*2)*d2x(cur_pos,node_array,nodes,c)) + sum(alpha(nodes:nodes*2)*d2y(cur_pos,node_array,nodes,c)))
           end if
          end do
        end select
        
        return
        
        end subroutine

        end program
