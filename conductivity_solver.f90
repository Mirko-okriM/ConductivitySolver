!Numerical solver for the calculation of effective conductivity kEff in porous media
! - diffusion type problem div(k*grad(T))=0 with k=k(x,y,z) is solved utilizing the finite volume method
! - two opposing sides are prescribed with a fixedValue Dirichlet boundary condition (of dT=1), remaining sides are zeroGradient Neumann boundary conditions
! - calculation of conducitivity at cell faces is performed with i) harmonic mean or ii) arithmetic mean (adjust in 'function kEval')
! - linear system is solved with the preconditioned conjugate gradient method using a JacobiDiagonalPreconditioner 
! - code is parallelized with openmp (number of processors has to be adjusted with environmental variable 'OMP_NUM_THREADS', e.g. powershell '$env:OMP_NUM_THREADS=8')
! - compile with 'gfortran.exe -fopenmp -O3 .\conductivity_solver.f90 -o conductivity_solver.exe'
! - output-file can be opened in paraview with 'loadData'-file (dimensions of sample have to be adjusted in loadData-file)
! - contact: mirko.siegert@web.de
! - Version 1.0 (28.11.23)
! - Changelog: [...]

module GlobalModule
    implicit none
    !integer, parameter :: chunkSize=100000	!activate if schedule(dynamic, chunkSize) is activated during parallel-looping
    !integer, parameter :: nProcessors=8
    integer, parameter :: nx = 100  										!Image width
    integer, parameter :: ny = 100  										!Image height
    integer, parameter :: nz = 100   										!Image depth
    integer, parameter :: voxelRes=1 										!VoxelResolution of input-raw data --> SHOULD STAY 1, VALUE DOES NOT REALLY HAVE INFLUENCE ON RESULT, MIGHT BE CHANGED TO TYPE REAL IN FUTURE VERSION
    integer, parameter :: realPrecision=2 									!floating-point precision: 1 (single) or 2 (double)
    integer, parameter :: nIteration=10000									!max. number of iterations
    real(kind=4*realPrecision), parameter :: targetResidual=1e-6						!target residual 	
    integer(kind=1), parameter :: bitLength = 8  								!Bites per pixel, typically grey-values in raw images are stored as either 8-bit or 16-bit integers
    integer, dimension(2) :: rawValues=[0,255] 									!bit-grayvalues of each phase (unsigned format)
    real(kind=4*realPrecision), dimension(2) :: correspConductivity=[0.0001,1.0] 				!conductivity of each phase (only double value, e.g. 1.0)
    character(512) :: file_path = 'C:\Users\Mirko\Desktop\solver_conductivity\veldsteen_poreus_100cube.raw' 	!input file (veldsteen_poreus_250cube, Berea_2d25um_binary_400cube)
    character(1) :: evalDirection = 'Z'										!direction of calculation (X or Y or Z)
    logical :: writeResults=.FALSE.										!write calculated temperature-field (.TRUE. or .FALSE.)
    character(512) :: output_path ='C:\Users\Mirko\Desktop\solver_conductivity\output.raw'			!name and path of the calculated temperature-field
    integer, parameter :: dx=voxelRes
    integer, parameter :: dy=voxelRes
    integer, parameter :: dz=voxelRes
    integer(kind=bitLength/8), allocatable :: currSample(:)

    !-Functions-------------------------------------------------------------------------------
    contains

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

    !functions for the implementation of the conjugate gradient method (ALS SUBROUTINE SCHREIBEN; SODASS DIE BERECHNUNG DIREKT IN DEN... )
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
				zVec(i)=r(i)/D(i)	!z^(k+1)=b_i/A_ii
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

end module GlobalModule



!-------------------------------------------------------------------------------------
!-----------MAIN----------------------------------------------------------------------
!-------------------------------------------------------------------------------------
program hello
    use omp_lib
    use GlobalModule
    implicit none
    integer :: x, y, z, i, j, currCell
    integer, parameter :: nCells=nx*ny*nz
    integer, parameter :: posODp1=1                 !cellPos(3,2,2)-cellPos(2,2,2)
    integer, parameter :: posODm1=-1                !cellPos(1,2,2)-cellPos(2,2,2)
    integer, parameter :: posODp2=nx                !cellPos(2,3,2)-cellPos(2,2,2)
    integer, parameter :: posODm2=-nx               !cellPos(2,1,2)-cellPos(2,2,2)
    integer, parameter :: posODp3=nx*ny             !cellPos(2,2,3)-cellPos(2,2,2)
    integer, parameter :: posODm3=-nx*ny            !cellPos(2,2,1)-cellPos(2,2,2)
    real(kind=4*realPrecision), allocatable :: MD(:), b(:), ODp1(:), ODp2(:), ODp3(:)!, ODm1(:), ODm2(:), ODm3(:) not necessary since Matrix is symmetric
    real(kind=4*realPrecision) :: kF, kB, kE, kW, kN, kS
    real(kind=4*realPrecision) :: aP, aF, aB, aE, aW, aN, aS, bCell
    real(kind=4*realPrecision) :: currResidual
    real(kind=4*realPrecision) :: totalFluxFront, totalFluxBack, totalFluxEast, totalFluxWest, totalFluxNorth, totalFluxSouth
    real(kind=4*realPrecision) :: aveFluxX, aveFluxY, aveFluxZ, kEffX, kEffY, kEffZ
    integer(kind=1) BC_F, BC_B, BC_E, BC_W, BC_N, BC_S
    character(8) :: patchTypeFront, patchTypeBack, patchTypeEast, patchTypeWest, patchTypeNorth, patchTypeSouth
    real(kind=4*realPrecision), allocatable :: r(:), p(:), xField(:), Apk(:), zVec(:) !variables for conjugate gradient method
    real(kind=4*realPrecision) :: alpha, beta, rDotProductOld, residualFactor

    !openmp might not work if giant arrays are not defined as allocatable variables
    allocate(currSample(nCells))
    allocate(MD(nCells))
    allocate(b(nCells))
    allocate(ODp1(nCells-posODp1))
    allocate(ODp2(nCells-posODp2))
    allocate(ODp3(nCells-posODp3))
    allocate(r(nCells))
    allocate(p(nCells))
    allocate(xField(nCells))
    allocate(Apk(nCells))
	allocate(zVec(nCells))

	!call omp_set_num_threads(nProcessors)

    !initialize variables for conjugate gradient method
    Apk=0.0
    r=0.0
    p=0.0
    xField=0.0

    !Setup boundary logic
    patchTypeFront='wall'   !x+
    patchTypeBack='wall'    !x-
    patchTypeEast='wall'    !y+
    patchTypeWest='wall'    !y-
    patchTypeNorth='wall'   !z+
    patchTypeSouth='wall'   !z-
    if (evalDirection=='Z') then
        patchTypeNorth='patch'  !z+
        patchTypeSouth='patch'  !z-
    elseif (evalDirection=='X') then
        patchTypeFront='patch'  !x+
        patchTypeBack='patch'   !x-
    elseif (evalDirection=='Y') then
        patchTypeEast='patch'   !y+
        patchTypeWest='patch'   !y-
    endif
    BC_F=0  !fixed front bc value if face is 'patch', if value is changed, than fluxCalculation might not be correct anymore
    BC_B=1  !fixed back bc value if face is 'patch', if value is changed, than fluxCalculation might not be correct anymore
    BC_E=0  !fixed east bc value if face is 'patch', if value is changed, than fluxCalculation might not be correct anymore
    BC_W=1  !fixed west bc value if face is 'patch', if value is changed, than fluxCalculation might not be correct anymore
    BC_N=0  !fixed north bc value if face is 'patch', if value is changed, than fluxCalculation might not be correct anymore
    BC_S=1  !fixed south bc value if face is 'patch', if value is changed, than fluxCalculation might not be correct anymore

    ! Open the file for binary reading
	write(*,*) '1) Reading raw-data ...'
    open(10, file=file_path, form='unformatted', access='stream', status='old')
	read(10) currSample
	close(10)
	write(*,*) '... done'
	
    !write(*,*) 'Wert von Position(20,20,20):',signedToUnsigned(currSample(cellPos(20,20,20)))

    !Convert rawValues to signed values
    do i = 1, size(rawValues)
        rawValues(i)=unsignedToSigned(rawValues(i))
        write(*,*) 'Wert ',i,': ',rawValues(i)
    end do

    !determination of matrix elements for each cell
	write(*,*) '2) Setting up linear system ...'
	!$omp parallel do collapse(3) private(kF, kB, kE, kW, kN, kS, aP, aF, aB, aE, aW, aN, aS, bCell, currCell)
		do z=1,nz
			do y=1,ny
				do x=1,nx
					currCell=cellPos(x,y,z);
					if (isInnerCell(x,y,z)) then !hardcoded equation result
						kF=kEval(x,y,z,x+1,y,z);
						kB=kEval(x,y,z,x-1,y,z);
						kE=kEval(x,y,z,x,y+1,z);
						kW=kEval(x,y,z,x,y-1,z);
						kN=kEval(x,y,z,x,y,z+1);
						kS=kEval(x,y,z,x,y,z-1);
						aP=-(kB+kW+kS+kF+kE+kN);
						MD(currCell)=aP;
						b(currCell)=0
						ODp1(currCell)=kF;
						!ODm1(currCell+posODm1)=kB; !not necessary since matrix is symmetric
						ODp2(currCell)=kE;
						!ODm2(currCell+posODm2)=kW; !not necessary since matrix is symmetric
						ODp3(currCell)=kN;
						!ODm3(currCell+posODm3)=kS; !not necessary since matrix is symmetric
					else !boundaryCell (code of this part might no be very elegant)
						aP=0; aF=0; aB=0; aE=0; aW=0; aN=0; aS=0;
						bCell=0;
						!evalFront
						if (x<nx) then !NB is cell
							kF=kEval(x,y,z,x+1,y,z);
							aP=aP-kF;
							aF=kF;
						else !NB is boundary
							!aP=aP+0; !default case (not necessary to evaluate)
							!aF=0;    !default case (not necessary to evaluate)
							if (patchTypeFront=='patch') then
								kF=kLookUp(x,y,z); !kF is equal to kP since cell is at boundary
								aP=aP-2*kF;
								!aF=0;                 !default case (not necessary to evaluate)
								bCell=-2*kF*BC_F;
							end if
						end if
						!evalBack
						if (x>1) then !NB is cell
							kB=kEval(x,y,z,x-1,y,z);
							aP=aP-kB;
							aB=kB;
						else !NB is boundary
							!aP=aP+0;
							!aB=0;
							if (patchTypeBack=="patch") then
								kB=kLookUp(x,y,z); !kB is equal to kP since cell is at boundary
								aP=aP-2*kB;
								!aB=0;                 !default case (not necessary to evaluate)
								bCell=-2*kB*BC_B;
							end if
						end if
						!evalEast
						if (y<ny) then !NB is cell
							kE=kEval(x,y,z,x,y+1,z);
							aP=aP-kE;
							aE=kE;
						else !NB is boundary
							!aP=aP+0;
							!aE=0;
							if (patchTypeEast=="patch") then
								kE=kLookUp(x,y,z); !kE is equal to kP since cell is at boundary
								aP=aP-2*kE;
								!aE=0;                 !default case (not necessary to evaluate)
								bCell=-2*kE*BC_E;
							end if
						end if
						!evalWest
						if (y>1) then !NB is cell
							kW=kEval(x,y,z,x,y-1,z);
							aP=aP-kW;
							aW=kW;
						else !NB is boundary
							!aP=aP+0;
							!aW=0;
							if (patchTypeWest=="patch") then
								kW=kLookUp(x,y,z); !kN is equal to kP since cell is at boundary
								aP=aP-2*kW;
								!aW=0;                 !default case (not necessary to evaluate)
								bCell=-2*kW*BC_W;
							end if
						end if
						!evalNorth
						if (z<nz) then !NB is cell
							kN=kEval(x,y,z,x,y,z+1);
							aP=aP-kN;
							aN=kN;
						else !NB is boundary
							!aP=aP+0;
							!aN=0;
							if (patchTypeNorth=="patch") then
								kN=kLookUp(x,y,z); !kN is equal to kP since cell is at boundary
								aP=aP-2*kN;
								!aN=0;                 !default case (not necessary to evaluate)
								bCell=-2*kN*BC_N;
							end if
						end if
						!evalSouth
						if (z>1) then !NB is cell
							kS=kEval(x,y,z,x,y,z-1);
							aP=aP-kS;
							aS=kS;
						else !NB is boundary
							!aP=aP+0;
							!aS=0;
							if (patchTypeSouth=="patch") then
								kS=kLookUp(x,y,z); !kS is equal to kP since cell is at boundary
								aP=aP-2*kS;
								!aS=0;                 !default case (not necessary to evaluate)
								bCell=-2*kS*BC_S;
							end if
						end if
						!store coefficients to matrix
						MD(currCell)=aP;
						b(currCell)=bCell;
						if (currCell<=nCells-posODp1) then   !first positiv offdiagonal is active
							ODp1(currCell)=aF;
						end if
	!                    if (posODp1<currCell) then          !first negativ offdiagonal is active, old condition (currCell>-posODm1)
	!                        ODm1(currCell+posODm1)=aB;      !not necessary since matrix is symmetric
	!                    end if
						if (currCell<=nCells-posODp2) then   !second positiv offdiagonal is active
							ODp2(currCell)=aE;
						end if
	!                    if (posODp2<currCell) then          !second negativ offdiagonal is active, old condition (currCell>-posODm2)
	!                        ODm2(currCell+posODm2)=aW;      !not necessary since matrix is symmetric
	!                    end if
						if (currCell<=nCells-posODp3) then   !third positiv offdiagonal is active
							ODp3(currCell)=aN;
						end if
	!                    if (posODp3<currCell) then          !third negativ offdiagonal is active, old condition (currCell>-posODm3)
	!                        ODm3(currCell+posODm3)=aS;      !not necessary since matrix is symmetric
	!                    end if
					end if
				end do
			end do
		end do
	!$omp end parallel do
	write(*,*) '... done'
	
	! call countZerosInVector(b)
	! call countZerosInVector(MD)
	! call countZerosInVector(ODp1)
	! call countZerosInVector(ODp2)
	! call countZerosInVector(ODp3)
	
	!preconditioned conjugate gradient with Jacobi-Diagonal-Preconditioner
	write(*,*) '3) Solving linear system ...'
	r=b																		!original logic would be r=b-A*x, since x0=zeroVector is chosen r=b follows
	deallocate(b)															!b while not be used anymore	
	call JacobiDiagonalPreconditioner(MD, r, zVec)							!z0=M^-1*r0
	p=zVec																	!p0=z0
	residualFactor=sumOfAbsComponents(r)
    currResidual=sumOfAbsComponents(r)/residualFactor
    write(*,*) 'Iteration 0 residual (L1-Norm): ', currResidual
    i=1
	do while (i<=nIteration .AND. currResidual>targetResidual)
		call symHeptaMatrixTimesVector(MD, ODp1, ODp2, ODp3, p, posODp1, posODp2, posODp3, Apk)			! calc Apk=A*pk
		rDotProductOld=dotProduct(r,zVec)																! calc rDotProductOld=r^T*z
		alpha=rDotProductOld/dotProduct(p,Apk) 															! calc alpha=(r^T*z)/(p^T*A*pk)
		call vectorPlusScalarTimesVector(xField,alpha,p,xField) 										! calc xNew=x+alpha*p
		call vectorPlusScalarTimesVector(r,-alpha,Apk,r)  												! calc rNew=r-alpha*A*p
		call JacobiDiagonalPreconditioner(MD, r, zVec)													! calc zNew=M^-1*rNew
		beta=dotProduct(r,zVec)/rDotProductOld  														! calc beta=(rNew^T*zNew)/(r^T*z)
		call vectorPlusScalarTimesVector(zVec,beta,p,p)													! calc pNew=zNew+beta*p
		currResidual=sumOfAbsComponents(r)/residualFactor
		!calc current conductivity
		if (evalDirection=='X') then
			totalFluxBack=calcFluxX(1,BC_B,xField)
			totalFluxFront=calcFluxX(nx,BC_F,xField)
			aveFluxX=(abs(totalFluxBack)+abs(totalFluxFront))/2;
			kEffX=aveFluxX/abs(BC_F-BC_B)*(nx*dx)/(ny*dy*nz*dz)
			write(*,*) 'Iteration ',i,' residual (L1-Norm): ', currResidual, ', kEffX= ',kEffX
		else if (evalDirection=='Y') then
			totalFluxWest=calcFluxY(1,BC_W,xField)
			totalFluxEast=calcFluxY(ny,BC_E,xField)
			aveFluxY=(abs(totalFluxWest)+abs(totalFluxEast))/2;
			kEffY=aveFluxY/abs(BC_E-BC_W)*(ny*dy)/(nx*dx*nz*dz)
			write(*,*) 'Iteration ',i,' residual (L1-Norm): ', currResidual, ', kEffY= ',kEffY
		else if (evalDirection=='Z') then
			totalFluxSouth=calcFluxZ(1,BC_S,xField)
			totalFluxNorth=calcFluxZ(nz,BC_N,xField)
			aveFluxZ=(abs(totalFluxSouth)+abs(totalFluxNorth))/2;
			kEffZ=aveFluxZ/abs(BC_N-BC_S)*(nz*dz)/(nx*dx*ny*dy)
			write(*,*) 'Iteration ',i,' residual (L1-Norm): ', currResidual, ', kEffZ= ',kEffZ
		end if
		i=i+1
	end do	
	write(*,*) '... done'

	
	write(*,*) 'Final result:'
	if (evalDirection=='X') then
		write(*,*) 'totalFluxFront= ', totalFluxFront, ', totalFluxBack= ', totalFluxBack, ', kEffX= ', kEffX
	else if (evalDirection=='Y') then
		write(*,*) 'totalFluxEast= ', totalFluxEast, ', totalFluxWest= ', totalFluxWest, ', kEffY= ', kEffY
	else if (evalDirection=='Z') then
		write(*,*) 'totalFluxNorth= ', totalFluxNorth, ', totalFluxSouth= ', totalFluxSouth, ', kEffZ= ', kEffZ
	end if

	if (writeResults) then
		write(*,*) 'Writing result file ...'
		open(1, file=output_path, access='stream', status='replace', form='unformatted')
		write(1) real(xField,4)
		close(1)
		write(*,*) '... done'
	end if

    !write(*,*) "Press Enter to exit..."
    !read(*,*)
end program
