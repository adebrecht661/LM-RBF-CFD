        module fluid_vars
        use declarations
        use Gaussian_RBF
        implicit none
        
        contains
        
        ! Fluid variables u, v, and p
        double precision function u(cur_pos)
        integer, intent(in) :: cur_pos
        
        u = sum(alpha(:nodes)*rbf(cur_pos,:))
        
        end function
        
        double precision function v(cur_pos)
        integer, intent(in) :: cur_pos
        
        v = sum(alpha(nodes+1:nodes*2)*rbf(cur_pos,:))
        
        end function
        
        double precision function p(cur_pos)
        integer, intent(in) :: cur_pos
        
        p = sum(alpha(nodes*2+1:nodes*3)*rbf(cur_pos,:))
        
        end function
        
        ! Derivatives in x direction
        double precision function dudx(cur_pos)
        integer, intent(in) :: cur_pos
        
        dudx = sum(alpha(:nodes)*drbfdx(cur_pos,:))
        
        end function
        
        double precision function dvdx(cur_pos)
        integer, intent(in) :: cur_pos
        
        dvdx = sum(alpha(nodes+1:nodes*2)*drbfdx(cur_pos,:))
        
        end function
        
        double precision function dpdx(cur_pos)
        integer, intent(in) :: cur_pos
        
        dpdx = sum(alpha(nodes*2+1:nodes*3)*drbfdx(cur_pos,:))
        
        end function
        
        ! Derivatives in y direction
        double precision function dudy(cur_pos)
        integer, intent(in) :: cur_pos
        
        dudy = sum(alpha(:nodes)*drbfdy(cur_pos,:))
        
        end function
        
        double precision function dvdy(cur_pos)
        integer, intent(in) :: cur_pos
        
        dvdy = sum(alpha(nodes+1:nodes*2)*drbfdy(cur_pos,:))
        
        end function
        
        double precision function dpdy(cur_pos)
        integer, intent(in) :: cur_pos
        
        dpdy = sum(alpha(nodes*2+1:nodes*3)*drbfdy(cur_pos,:))
        
        end function
        
        ! Second derivatives in x direction
        double precision function d2udx2(cur_pos)
        integer, intent(in) :: cur_pos
        
        d2udx2 = sum(alpha(:nodes)*d2rbfdx2(cur_pos,:))
        
        end function

        double precision function d2vdx2(cur_pos)
        integer, intent(in) :: cur_pos
        
        d2vdx2 = sum(alpha(nodes+1:nodes*2)*d2rbfdx2(cur_pos,:))
        
        end function
        
        double precision function d2pdx2(cur_pos)
        integer, intent(in) :: cur_pos
        
        d2pdx2 = sum(alpha(nodes*2+1:nodes*3)*d2rbfdx2(cur_pos,:))
        
        end function
        
        ! Second derivatives in y direction
        double precision function d2udy2(cur_pos)
        integer, intent(in) :: cur_pos
        
        d2udy2 = sum(alpha(:nodes)*d2rbfdy2(cur_pos,:))
        
        end function
        
        double precision function d2vdy2(cur_pos)
        integer, intent(in) :: cur_pos
        
        d2vdy2 = sum(alpha(nodes+1:nodes*2)*d2rbfdy2(cur_pos,:))
        
        end function
        
        double precision function d2pdy2(cur_pos)
        integer, intent(in) :: cur_pos
        
        d2pdy2 = sum(alpha(nodes*2+1:nodes*3)*d2rbfdy2(cur_pos,:))
        
        end function
        
        end module
