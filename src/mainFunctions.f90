    integer(kind=bitLength/8) function unsignedToSigned(greyValue) !should work for input integers kind<4 (i.e. 32 bit)
        integer :: greyValue
        if (0 <= greyValue .AND. greyValue <= (2**bitLength)/2-1) then
            unsignedToSigned=greyValue
        else
            unsignedToSigned=greyValue-(2**bitLength)
        end if
    end function unsignedToSigned

    integer function signedToUnsigned(greyValue) !should work for input integers kind<4 (i.e. 32 bit)
        integer(kind=bitLength/8) :: greyValue
        if (0 <= greyValue .AND. greyValue <= (2**bitLength)/2-1) then
            signedToUnsigned=greyValue
        else
            signedToUnsigned=greyValue+(2**bitLength)
        end if
    end function signedToUnsigned

    integer function cellPos(x,y,z)
        integer :: x,y,z
        cellPos=nx*ny*(z-1)+nx*(y-1)+x
    end function

    logical function isInnerCell(x,y,z)
        integer :: x,y,z
        isInnerCell=.TRUE.
        if (x == 1 .OR. x == nx .OR. y == 1 .OR. y == ny .OR. z == 1 .OR. z == nz) then
            isInnerCell=.FALSE.
        end if
    end function

    real(kind=4*realPrecision) function kLookUp(x,y,z)
        integer :: x,y,z,i
        do i = 1, size(rawValues)
            if (currSample(cellPos(x, y, z)) == rawValues(i)) then
              kLookUp = correspConductivity(i)
              exit ! Exit the loop when the condition is met
            end if
        end do
    end function

    real(kind=4*realPrecision) function kEval(x1,y1,z1,x2,y2,z2)
        integer :: x1,y1,z1,x2,y2,z2
        real(kind=4*realPrecision) :: kP, kNB
        kP=kLookUp(x1,y1,z1)
        kNB=kLookUp(x2,y2,z2)
        kEval=2*kP*kNB/(kP+kNB)  !harmonic
        !kEval=(kP+kNB)/2         !arithmetic
    end function
    
    subroutine countZerosInVector(v)
        real(kind=4*realPrecision), dimension(:) :: v
        integer :: i, nZeros
        nZeros=0
        !$omp parallel do reduction(+:nZeros)
        do i=1,size(v)
            if (v(i)==0) then
                nZeros=nZeros+1
            end if
        end do
        !$omp end parallel do
        write(*,*) 'nSize: ',size(v), 'nZeros: ', nZeros
    end subroutine

    !functions for the implementation of the conjugate gradient method 
    real(kind=4*realPrecision) function dotProduct(v1, v2)
        real(kind=4*realPrecision), dimension(:) :: v1, v2
        integer :: i
        dotProduct=0
        !$omp parallel do reduction(+:dotProduct) !schedule(dynamic, chunkSize)
            do i=1,size(v1)
                dotProduct=dotProduct+v1(i)*v2(i)
            end do
        !$omp end parallel do
    end function
    
    subroutine JacobiDiagonalPreconditioner(D, r, zVec) !Jacobi-Diagonal-Preconditioner
        real(kind=4*realPrecision), dimension(:) :: D, r, zVec
        integer :: i
        !$omp parallel do !schedule(dynamic, chunkSize)
            do i=1,size(zVec)
                zVec(i)=r(i)/D(i)   !z^(k+1)=b_i/A_ii
            end do
        !$omp end parallel do
    end subroutine JacobiDiagonalPreconditioner
    
    subroutine JacobiPreconditioner(MD, OD1, OD2, OD3, b, zOld, posOD1, posOD2, posOD3, zVec) !only valid for symmetric heptadiagonal matrices
        real(kind=4*realPrecision), dimension(:) :: MD, OD1, OD2, OD3, b, zOld, zVec
        integer :: posOD1, posOD2, posOD3, i
        integer, parameter :: nCells=nx*ny*nz
        !$omp parallel do !schedule(dynamic, chunkSize)
            do i=1,nx*ny*nz !logic might be enhanced with less if statements, but maybe code readability might suffer
                zVec(i)=b(i)
                if (i<=nCells-posOD1) then               !first positiv offdiagonal is active
                    zVec(i)=zVec(i)-OD1(i)*zOld(i+posOD1)
                end if
                if (i<=nCells-posOD2) then               !second positiv offdiagonal is active
                    zVec(i)=zVec(i)-OD2(i)*zOld(i+posOD2)
                end if
                if (i<=nCells-posOD3) then               !third positiv offdiagonal is active
                    zVec(i)=zVec(i)-OD3(i)*zOld(i+posOD3)
                end if
                if (posOD1<i) then                       !first negativ offdiagonal is active
                    zVec(i)=zVec(i)-OD1(i-posOD1)*zOld(i-posOD1)
                end if
                if (posOD2<i) then                       !second negativ offdiagonal is active
                    zVec(i)=zVec(i)-OD2(i-posOD2)*zOld(i-posOD2)
                end if
                if (posOD3<i) then                       !third negativ offdiagonal is active
                    zVec(i)=zVec(i)-OD3(i-posOD3)*zOld(i-posOD3)
                end if
                zVec(i)=zVec(i)/MD(i)
            end do
        !$omp end parallel do
    end subroutine      

    subroutine symHeptaMatrixTimesVector(MD, OD1, OD2, OD3, v, posOD1, posOD2, posOD3, resVec) !only valid for symmetric heptadiagonal matrices
        real(kind=4*realPrecision), dimension(:) :: MD, OD1, OD2, OD3, v, resVec
        integer :: posOD1, posOD2, posOD3, i
        integer, parameter :: nCells=nx*ny*nz
        !$omp parallel do !schedule(dynamic, chunkSize)
            do i=1,nx*ny*nz !logic might be enhanced with less if statements, but maybe code readability might suffer
                resVec(i)=MD(i)*v(i)
                if (i<=nCells-posOD1) then               !first positiv offdiagonal is active
                    resVec(i)=resVec(i)+OD1(i)*v(i+posOD1)
                end if
                if (i<=nCells-posOD2) then               !second positiv offdiagonal is active
                    resVec(i)=resVec(i)+OD2(i)*v(i+posOD2)
                end if
                if (i<=nCells-posOD3) then               !third positiv offdiagonal is active
                    resVec(i)=resVec(i)+OD3(i)*v(i+posOD3)
                end if
                if (posOD1<i) then                       !first negativ offdiagonal is active
                    resVec(i)=resVec(i)+OD1(i-posOD1)*v(i-posOD1)
                end if
                if (posOD2<i) then                       !second negativ offdiagonal is active
                    resVec(i)=resVec(i)+OD2(i-posOD2)*v(i-posOD2)
                end if
                if (posOD3<i) then                       !third negativ offdiagonal is active
                    resVec(i)=resVec(i)+OD3(i-posOD3)*v(i-posOD3)
                end if
            end do
        !$omp end parallel do
    end subroutine symHeptaMatrixTimesVector

    subroutine vectorPlusScalarTimesVector(v1, s, v2, resVec)
        real(kind=4*realPrecision), dimension(:) :: v1, v2, resVec
        real(kind=4*realPrecision) :: s
        integer :: i
        !$omp parallel do !schedule(dynamic, chunkSize)
            do i=1,size(v1)
                resVec(i)=v1(i)+s*v2(i)
            end do
        !$omp end parallel do
    end subroutine vectorPlusScalarTimesVector

    real(kind=4*realPrecision) function sumOfAbsComponents(v1)
        real(kind=4*realPrecision), dimension(:) :: v1
        integer :: i
        sumOfAbsComponents=0
        !$omp parallel do reduction(+:sumOfAbsComponents) !schedule(dynamic, chunkSize)
            do i=1,size(v1)
                sumOfAbsComponents=sumOfAbsComponents+abs(v1(i))
            end do
        !$omp end parallel do
    end function

    !functions for postprocessing (determination of heatFlux)
    real(kind=4*realPrecision) function calcFluxX(xTarget,BC,xField)
        integer :: y, z, xTarget, currCell
        integer(kind=1) :: BC
        real(kind=4*realPrecision) :: currK
        real(kind=4*realPrecision), dimension(:) :: xField
        calcFluxX=0.0
        !$omp parallel do collapse(2) private(currK, currCell) reduction(+:calcFluxX)
        do z=1,nz
            do y=1,ny
                currK=kLookUp(xTarget,y,z)
                currCell=cellPos(xTarget,y,z)
                calcFluxX=calcFluxX-currK*(BC-xField(currCell))
            end do
        end do              
        !$omp end parallel do
        calcFluxX=-calcFluxX*(2*dy*dz/dx);
    end function

    real(kind=4*realPrecision) function calcFluxY(yTarget,BC,xField)
        integer :: x, z, yTarget, currCell
        integer(kind=1) :: BC
        real(kind=4*realPrecision) :: currK
        real(kind=4*realPrecision), dimension(:) :: xField
        calcFluxY=0.0
        !$omp parallel do collapse(2) private(currK, currCell) reduction(+:calcFluxY)
        do z=1,nz
            do x=1,nx
                currK=kLookUp(x,yTarget,z)
                currCell=cellPos(x,yTarget,z)
                calcFluxY=calcFluxY-currK*(BC-xField(currCell))
            end do
        end do              
        !$omp end parallel do
        calcFluxY=-calcFluxY*(2*dx*dz/dy);
    end function
    
    real(kind=4*realPrecision) function calcFluxZ(zTarget,BC,xField)
        integer :: x, y, zTarget, currCell
        integer(kind=1) :: BC
        real(kind=4*realPrecision) :: currK
        real(kind=4*realPrecision), dimension(:) :: xField
        calcFluxZ=0.0
        !$omp parallel do collapse(2) private(currK, currCell) reduction(+:calcFluxZ)
        do y=1,ny
            do x=1,nx
                currK=kLookUp(x,y,zTarget)
                currCell=cellPos(x,y,zTarget)
                calcFluxZ=calcFluxZ-currK*(BC-xField(currCell))
            end do
        end do      
        !$omp end parallel do
        calcFluxZ=-calcFluxZ*(2*dx*dy/dz);
    end function