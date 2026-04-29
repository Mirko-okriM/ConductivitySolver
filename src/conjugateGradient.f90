!initialize variables for conjugate gradient method
Apk=0.0d0
r=0.0d0
p=0.0d0
xField=0.0d0
!set G-vector up
Gx=0
Gy=0
Gz=0
if (evalDirection=='X') then !logic explained in handwritten notes
    Gx=1
else if (evalDirection=='Y') then
    Gy=1
else if (evalDirection=='Z') then
    Gz=1
end if 
outputInfoCounter=0
!preconditioned conjugate gradient with Jacobi-Diagonal-Preconditioner
r=b                                                                     !original logic would be r=b-A*x, since x0=zeroVector is chosen r=b follows
!deallocate(b)                                                          !b while not be used anymore    
call JacobiDiagonalPreconditioner(MD, r, zVec)                          !z0=M^-1*r0
p=zVec                                                                  !p0=z0
residualFactor=sumOfAbsComponents(r)
currResidual=sumOfAbsComponents(r)/residualFactor
write(*,*) 'Iteration 0 residual (L1-Norm): ', currResidual
i=1
do while (i<=nIteration .AND. currResidual>targetResidual)
    call CSRMatrixTimesVector(A_mat, colInd, rowPtr, p, Apk, nCellsRock)         ! calc Apk=A*pk
    rDotProductOld=dotProduct(r,zVec)                                                               ! calc rDotProductOld=r^T*z
    alpha=rDotProductOld/dotProduct(p,Apk)                                                          ! calc alpha=(r^T*z)/(p^T*A*pk)
    call vectorPlusScalarTimesVector(xField,alpha,p,xField)                                         ! calc xNew=x+alpha*p
    call vectorPlusScalarTimesVector(r,-alpha,Apk,r)                                                ! calc rNew=r-alpha*A*p
    call JacobiDiagonalPreconditioner(MD, r, zVec)                                                  ! calc zNew=M^-1*rNew
    beta=dotProduct(r,zVec)/rDotProductOld                                                          ! calc beta=(rNew^T*zNew)/(r^T*z)
    call vectorPlusScalarTimesVector(zVec,beta,p,p)                                                 ! calc pNew=zNew+beta*p
    currResidual=sumOfAbsComponents(r)/residualFactor
    xAve=sumOfComponents(xField)/nCells
    call forceAveZeroFieldConstraint(xField, xAve, xField) !necessary? ChatGPT says yes ...? 
    if  (outputInfoCounter == calcFluxStep) then
        outputInfoCounter=0
        if (searchDirection=='forwardSearch') then
            include './src/fluxCalculation_forwardSearch.f90' 
        else if (searchDirection=='backwardSearch') then   
            include './src/fluxCalculation_backwardSearch.f90' 
        end if  
        !Solver feedback
        if (evalDirection=='X') then
            write(*,*) 'Iteration ', i,' residual (L1-Norm): ', currResidual, ', k_xx= ', qAveX 
        else if (evalDirection=='Y') then
            write(*,*) 'Iteration ', i,' residual (L1-Norm): ', currResidual, ', k_yy= ', qAveY
        else if (evalDirection=='Z') then
            write(*,*) 'Iteration ', i,' residual (L1-Norm): ', currResidual, ', k_zz= ', qAveZ
        end if        
    else
    if (evalDirection=='X') then
            write(*,*) 'Iteration ', i,' residual (L1-Norm): ', currResidual
        else if (evalDirection=='Y') then
            write(*,*) 'Iteration ', i,' residual (L1-Norm): ', currResidual
        else if (evalDirection=='Z') then
            write(*,*) 'Iteration ', i,' residual (L1-Norm): ', currResidual
        end if
    end if    
    i=i+1
    outputInfoCounter=outputInfoCounter+1
end do  

if (searchDirection=='forwardSearch') then
        include './src/fluxCalculation_forwardSearch.f90' 
else if (searchDirection=='backwardSearch') then   
        include './src/fluxCalculation_backwardSearch.f90' 
end if     
