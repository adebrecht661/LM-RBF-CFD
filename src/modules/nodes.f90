        module grid
        use declarations
        implicit none
        
        type node
         sequence
         integer x,y
        end type node
        
        type(node), allocatable :: node_array(:)
        
        contains
        
        subroutine SetUpGrid()
        integer i,j,n
        type(node), allocatable :: temp(:)
        
        i = 0
        j = 0
        allocate(node_array(nodes))
        
!       define node locations : x and y go from 0 to dim-1 [(0,0) to (x_dim-1,y_dim-1)]  
        do n = 1, nodes
         node_array(n)%x = mod(i,x_dim)
         node_array(n)%y = mod(j,y_dim)
         i = i + 1
         if (mod(i,x_dim) == 0 .and. i > 0) j = j + 1
        end do
        
        select case (BC)
         case (channelCent)
          temp = pack(node_array, .not. (node_array%y > floor(y_dim/3.) .and. node_array%y < floor(2*y_dim/3.) .and. node_array%x > floor(x_dim/3.) .and. node_array%x < floor(2*x_dim/3.)))
          nodes = size(temp)
          deallocate(node_array)
          allocate(node_array(nodes))
          node_array = temp
          deallocate(temp)
         case (channelObs)
          temp = pack(node_array, .not. (node_array%y < floor(y_dim/2.) .and. node_array%x > floor(x_dim/3.) .and. node_array%x < floor(2*x_dim/3.)))
          nodes = size(temp)
          deallocate(node_array)
          allocate(node_array(nodes))
          node_array = temp
          deallocate(temp)
         case (channelCircle)
          temp = pack(node_array, .not. (((node_array%y-cent(2))**2 + (node_array%x-cent(1))**2) < r2))
          nodes = size(temp)
          deallocate(node_array)
          allocate(node_array(nodes))
          node_array = temp
          deallocate(temp)
        end select
        
        return
        
        end subroutine
        
        end module
