        module declarations
        implicit none
        
        integer x_dim,y_dim,nodes,nvar,nfunc,info,lwa,BC,cent(2),r2,res
        double precision tol,uleft,vleft,pleft,uwall,uright,vright,pright,Re,c,alpha_max,alpha_min,width
        logical right_outflow, right_set
        
        integer, allocatable :: iwa(:)
        double precision, allocatable :: alpha(:),fvec(:),wa(:),rbf(:,:),drbfdy(:,:),drbfdx(:,:),d2rbfdx2(:,:),d2rbfdy2(:,:)
        
        integer, parameter :: couette = 2, channel = 1, channelCent = 3, channelObs = 4, channelCircle = 5
        
        end module
