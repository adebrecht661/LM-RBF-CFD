        module TPS_RBF
        use grid
        implicit none
        private
!       Public functions return 
        public psi, dx, dy, d2x, d2y, CalcShapeParam
       
        contains
        
        function d2y(cur_pos,node_array,nodes,c) result(deriv)
        
        integer nodes,cur_pos,i
        double precision c
        double precision deriv(nodes)
        type(node) node_array(nodes)
        
        do i = 1,nodes
         deriv(i) = d2rbfdy2(node_array(cur_pos),node_array(i),c)
        end do
        
        end function
        
        function d2x(cur_pos,node_array,nodes,c) result(deriv)
        
        integer nodes,cur_pos,i
        double precision c
        double precision deriv(nodes)
        type(node) node_array(nodes)
        
        do i = 1,nodes
         deriv(i) = d2rbfdx2(node_array(cur_pos),node_array(i),c)
        end do
        
        end function
        
        function dy(cur_pos,node_array,nodes,c) result(deriv)
        
        integer nodes,cur_pos,i
        double precision c
        double precision deriv(nodes)
        type(node) node_array(nodes)
        
        do i = 1,nodes
         deriv(i) = drbfdy(node_array(cur_pos),node_array(i),c)
        end do
        
        end function
        
        function dx(cur_pos,node_array,nodes,c) result(deriv)
        
        integer nodes,cur_pos,i
        double precision c
        double precision deriv(nodes)
        type(node) node_array(nodes)
        
        do i = 1,nodes
         deriv(i) = drbfdx(node_array(cur_pos),node_array(i),c)
        end do
        
        end function
        
        function psi(cur_pos,node_array,nodes,c) result(psis)
        
        integer nodes,cur_pos,i
        double precision c
        double precision psis(nodes)
        type(node) node_array(nodes)
        
        do i = 1,nodes
         psis(i) = rbf(node_array(cur_pos),node_array(i),c)
        end do
        
        end function
        
        subroutine CalcShapeParam(res,c)
        
        integer res
        double precision c
        
        c = 2000*sqrt(2.)*(1./res)

        return
                
        end subroutine

!       Functions for calculating value of radial basis function and its derivatives
        double precision function rbf(cur_node,node_i,c)
        type(node) :: cur_node,node_i
        double precision c
        integer r2
        
        r2 = (node_i%x-cur_node%x)**2 + (node_i%y-cur_node%y)**2

        if (r2 == 0) then
         rbf = 0
        else
         rbf = log(sqrt(dble(r2)))*r2**2
        end if
        
        end function
        
        double precision function drbfdx(cur_node,node_i,c)
        type(node) :: cur_node,node_i
        double precision c
        integer r2
        
        r2 = (node_i%x-cur_node%x)**2 + (node_i%y-cur_node%y)**2
        
        if (r2 == 0) then
         drbfdx = 0
        else
         drbfdx = -(node_i%x-cur_node%x)*r2*(4*log(sqrt(dble(r2)))+1)
        end if
        
        end function
        
        double precision function drbfdy(cur_node,node_i,c)
        type(node) :: cur_node,node_i
        double precision c
        integer r2
        
        r2 = (node_i%x-cur_node%x)**2 + (node_i%y-cur_node%y)**2
        
        if (r2 == 0) then
         drbfdy = 0
        else
         drbfdy = -(node_i%y-cur_node%y)*r2*(4*log(sqrt(dble(r2)))+1)
        end if
        
        end function        

        double precision function d2rbfdx2(cur_node,node_i,c)
        type(node) :: cur_node,node_i
        double precision c
        integer r2
        
        r2 = (node_i%x-cur_node%x)**2 + (node_i%y-cur_node%y)**2
        
        if (r2 == 0) then
         d2rbfdx2 = 0
        else
         d2rbfdx2 = 8*(node_i%x-cur_node%x)**2*log(sqrt(dble(r2))) - 4*r2*log(sqrt(dble(r2))) + 6*(node_i%x-cur_node%x)**2 - r2
        end if
        
        end function

        double precision function d2rbfdy2(cur_node,node_i,c)
        type(node) :: cur_node,node_i
        double precision c
        integer r2
        
        r2 = (node_i%x-cur_node%x)**2 + (node_i%y-cur_node%y)**2
        
        if (r2 == 0) then
         d2rbfdy2 = 0
        else
         d2rbfdy2 = 8*(node_i%y-cur_node%y)**2*log(sqrt(dble(r2))) - 4*r2*log(sqrt(dble(r2))) + 6*(node_i%y-cur_node%y)**2 - r2
        end if
        
        end function
        
        end module
