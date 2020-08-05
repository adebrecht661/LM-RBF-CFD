        module Gaussian_RBF
        use declarations
        use grid
        implicit none
       
        contains
        
        subroutine CalcRBF()
        integer i,j
        
        allocate(rbf(nodes,nodes),drbfdy(nodes,nodes),drbfdx(nodes,nodes),d2rbfdx2(nodes,nodes),d2rbfdy2(nodes,nodes))
        
        do j = 1, nodes
         do i = 1, nodes
          r2 = (node_array(j)%x-node_array(i)%x)**2 + (node_array(j)%y-node_array(i)%y)**2
          rbf(i,j) = exp((-r2)/c)
          drbfdy(i,j) = -2.*(node_array(j)%y - node_array(i)%y)*rbf(i,j)/c
          drbfdx(i,j) = -2.*(node_array(j)%x - node_array(i)%x)*rbf(i,j)/c
          d2rbfdy2(i,j) = -2.*rbf(i,j)/c + rbf(i,j)*4.*((node_array(j)%y - node_array(i)%y)/c)**2
          d2rbfdx2(i,j) = -2.*rbf(i,j)/c + rbf(i,j)*4.*((node_array(j)%x - node_array(i)%x)/c)**2
         end do
        end do
        
        return
        
        end subroutine
        
        subroutine CalcShapeParam
        
        c = width*sqrt(2.)*(1./res)

        return
                
        end subroutine
        
        end module
