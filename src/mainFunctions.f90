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
    
    ! integer(kind=4) function rowPtr(currRow) !since all rows have 7 elements, rowPtr can be calculated on the fly and memory is saved
        ! integer (kind=4) :: currRow
        ! rowPtr=7*(currRow-1)+1
    ! end function
    
    subroutine bubbleSort(keys, vals) !bubble sort for two arrays, the keys array is sorted and the position of the vals-elements is set accordingly
        integer, intent(inout) :: keys(:)
        real(kind=4*realPrecision), intent(inout) :: vals(:)
        integer :: i, j, n, tk
        real(kind=4*realPrecision) :: tv
        n = size(keys)
        do i = 1, n-1
            do j = 1, n-i
                if (keys(j) > keys(j+1)) then
                    tk = keys(j); keys(j) = keys(j+1); keys(j+1) = tk
                    tv = vals(j); vals(j) = vals(j+1); vals(j+1) = tv
                end if
            end do
        end do
    end subroutine
    
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
    
    subroutine CSRMatrixTimesVector(A_mat, colInd, rowPtr, vIn, resVec, nCells)
        real(kind=4*realPrecision), dimension(:) :: A_mat,  vIn, resVec
        integer(kind=4), dimension(:) :: colInd, rowPtr
        integer :: i, nCells, currRow
        real(kind=4*realPrecision) :: sum_tmp
        resVec=0.0d0
        !$omp parallel do default(shared) private(sum_tmp) !schedule(dynamic, chunkSize)
            do currRow=1,nCells
            sum_tmp = 0.0d0
            !$omp simd reduction(+:sum_tmp) 
            !according to chatgpt "simd reduction(+:sum)" improves chache management and 
            !calculation is performed in one step (using vectorregisters to speed up calculation)
                do i=rowPtr(currRow),rowPtr(currRow+1)-1
                    sum_tmp=sum_tmp+A_mat(i)*vIn(colInd(i))
                end do
                resVec(currRow)=sum_tmp
            end do
        !$omp end parallel do
    end subroutine CSRMatrixTimesVector

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
    
    real(kind=4*realPrecision) function sumOfComponents(v1)
        real(kind=4*realPrecision), dimension(:) :: v1
        integer :: i
        sumOfComponents=0
        !$omp parallel do reduction(+:sumOfComponents) !schedule(dynamic, chunkSize)
            do i=1,size(v1)
                sumOfComponents=sumOfComponents+v1(i)
            end do
        !$omp end parallel do
    end function
    
    subroutine forceAveZeroFieldConstraint(vecIn, s, resVec)
        real(kind=4*realPrecision), dimension(:) :: vecIn, resVec
        real(kind=4*realPrecision) :: s
        integer :: i
        !$omp parallel do !schedule(dynamic, chunkSize)
            do i=1,size(vecIn)
                resVec(i)=vecIn(i)-s
            end do
        !$omp end parallel do
    end subroutine forceAveZeroFieldConstraint

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