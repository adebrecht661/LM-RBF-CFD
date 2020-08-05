        module input
        use declarations
        use Gaussian_RBF
        use grid
        implicit none
        
        contains
        
        subroutine GetProblemInfo()
        double precision cent_real(2),r2_real
        logical ex
        
        NAMELIST/ProblemData/ x_dim,y_dim,res,tol,Re,BC,r2_real,cent_real,uleft,vleft,pleft,uwall,alpha_max,alpha_min,width,right_outflow,right_set,uright,vright,pright
        
        inquire(file='problem.data', exist=ex)
        if (ex) then
         OPEN(UNIT=30, FILE='problem.data', STATUS="OLD")
         READ(30,NML=ProblemData)
         CLOSE(30)
         select case (BC)
          case (channel)
           print *, 'Using channel flow boundary conditions.'
          case (couette)
           print *, 'Using Couette flow boundary conditions.'
          case (channelCent)
           print *, 'Using channel flow boundary conditions with central obstruction (center square third of grid).'
          case (channelObs)
           print *, 'Using channel flow boundary conditions with obstruction.'
          case (channelCircle)
           print *, 'Using channel flow boundary conditions with spherical obstruction.'
         end select
        else
         print *, 'Missing problem.data. See README for more information.'
        end if
        
        r2 = floor(r2_real*res**2)
        cent = floor(cent_real*res)       
        
        x_dim = x_dim*res
        y_dim = y_dim*res
        
        nodes = x_dim*y_dim
        nvar = 3*nodes          ! Three variables (u,v,p) at each point
        nfunc = nvar            ! nfunc >= nvar, need at least three functions per point
        
        lwa = nfunc*nvar+5*nvar+nfunc
        
        call CalcShapeParam
        call SetUpGrid
        call CalcRBF
        
        return
        
        end subroutine
        
        subroutine GetAlpha(alpha)
        double precision, intent(out) :: alpha(nvar)
        integer i
        logical ex
        
        inquire(file='alpha.dat',exist=ex)
        if (ex) then
         open(unit=30, file='alpha.dat',status='old')
         do i = 1, nvar
          read (30,*) alpha(i)
         end do
        else
         call random_seed
         call random_number(alpha)
         alpha = alpha*(alpha_max - alpha_min) + alpha_min
        end if
        
        end subroutine
        
        end module
