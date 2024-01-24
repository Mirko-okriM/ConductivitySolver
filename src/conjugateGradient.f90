!preconditioned conjugate gradient with Jacobi-Diagonal-Preconditioner
r=b                                                                     !original logic would be r=b-A*x, since x0=zeroVector is chosen r=b follows
deallocate(b)                                                           !b while not be used anymore    
call JacobiDiagonalPreconditioner(MD, r, zVec)                          !z0=M^-1*r0
p=zVec                                                                  !p0=z0
residualFactor=sumOfAbsComponents(r)
currResidual=sumOfAbsComponents(r)/residualFactor
write(*,*) 'Iteration 0 residual (L1-Norm): ', currResidual
i=1
do while (i<=nIteration .AND. currResidual>targetResidual)
    call symHeptaMatrixTimesVector(MD, ODp1, ODp2, ODp3, p, posODp1, posODp2, posODp3, Apk)         ! calc Apk=A*pk
    rDotProductOld=dotProduct(r,zVec)                                                               ! calc rDotProductOld=r^T*z
    alpha=rDotProductOld/dotProduct(p,Apk)                                                          ! calc alpha=(r^T*z)/(p^T*A*pk)
    call vectorPlusScalarTimesVector(xField,alpha,p,xField)                                         ! calc xNew=x+alpha*p
    call vectorPlusScalarTimesVector(r,-alpha,Apk,r)                                                ! calc rNew=r-alpha*A*p
    call JacobiDiagonalPreconditioner(MD, r, zVec)                                                  ! calc zNew=M^-1*rNew
    beta=dotProduct(r,zVec)/rDotProductOld                                                          ! calc beta=(rNew^T*zNew)/(r^T*z)
    call vectorPlusScalarTimesVector(zVec,beta,p,p)                                                 ! calc pNew=zNew+beta*p
    currResidual=sumOfAbsComponents(r)/residualFactor
    !calc current conductivity (code should be improved ... if-condition in every loop is not very smart)
    if (evalDirection=='X') then
        totalFluxBack=calcFluxX(1,BC_B,xField)
        totalFluxFront=calcFluxX(nx,BC_F,xField)
        aveFluxX=(abs(totalFluxBack)+abs(totalFluxFront))/2
        kEffX=aveFluxX/abs(BC_F-BC_B)*(nx*dx)/(ny*dy*nz*dz)
        write(*,*) 'Iteration ', i,' residual (L1-Norm): ', currResidual, ', kEffX= ',kEffX
    else if (evalDirection=='Y') then
        totalFluxWest=calcFluxY(1,BC_W,xField)
        totalFluxEast=calcFluxY(ny,BC_E,xField)
        aveFluxY=(abs(totalFluxWest)+abs(totalFluxEast))/2
        kEffY=aveFluxY/abs(BC_E-BC_W)*(ny*dy)/(nx*dx*nz*dz)
        write(*,*) 'Iteration ', i,' residual (L1-Norm): ', currResidual, ', kEffY= ',kEffY
    else if (evalDirection=='Z') then
        totalFluxSouth=calcFluxZ(1,BC_S,xField)
        totalFluxNorth=calcFluxZ(nz,BC_N,xField)
        aveFluxZ=(abs(totalFluxSouth)+abs(totalFluxNorth))/2
        kEffZ=aveFluxZ/abs(BC_N-BC_S)*(nz*dz)/(nx*dx*ny*dy)
        write(*,*) 'Iteration ', i,' residual (L1-Norm): ', currResidual, ', kEffZ= ',kEffZ
    end if
    i=i+1
end do  