        !!! Current status - modified, compiles, not debugged
        program fluid_solver
        use declarations
        use LM
        use grid
!        use TPS_RBF            !!!! Not working yet
        use input
        use output
        use boundaries
        implicit none

!       get grid size and resolution for input to lmdif
        call GetProblemInfo

!       Calculate length of working array and allocate matrices
        allocate(alpha(nvar),iwa(nvar),fvec(nfunc),wa(lwa))

        call GetAlpha(alpha)

!       solve with system f(x) = 0 only
        call lmdif1(fcn,nfunc,nvar,alpha,fvec,tol,info,iwa,wa,lwa)

        print *, 'Completed calculation. Writing out.'
        call ResultToTxt

        deallocate(alpha,iwa,fvec,wa)

        stop

        end program
