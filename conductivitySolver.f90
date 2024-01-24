!Numerical solver for the calculation of effective conductivity kEff in porous media
! - diffusion type problem div(k*grad(T))=0 with k=k(x,y,z) is solved utilizing the finite volume method
! - two opposing sides are prescribed with a fixedValue Dirichlet boundary condition (of dT=1), remaining sides are zeroGradient Neumann boundary conditions
! - calculation of conducitivity at cell faces is performed with i) harmonic mean or ii) arithmetic mean (adjust in 'function kEval')
! - linear system is solved with the preconditioned conjugate gradient method using a JacobiDiagonalPreconditioner 
! - code is parallelized with openmp (number of processors has to be adjusted with environmental variable 'OMP_NUM_THREADS', e.g. powershell '$env:OMP_NUM_THREADS=8')
! - compile with                 
!                gfortran.exe -fopenmp -O3 conductivitySolver.f90 -o solverRun.exe
! - contact: mirko.siegert@web.de
! - Version 1.0 (28.11.23)
! - Changelog: 
!   Version 1.1 (24.01.24): 
!       - code is organized in seperate files 
!       - added flux field computation for postprocessing 
!       - added automated paraview load file for postprocessing 
!       - adjustments for user-friendliness 
!

!-------------------------------------------------------------------------------------
!-----------MAIN----------------------------------------------------------------------
!-------------------------------------------------------------------------------------
program conductivitySolver
    use omp_lib 
    implicit none
    !USER INPUT
    !integer, parameter :: chunkSize=100000 !activate if schedule(dynamic, chunkSize) is activated during parallel-looping
    !integer, parameter :: nProcessors=8
    include './settings.f90'    
    !additional variables
    character(512) :: pathTemperature, pathFlux, pathPost, pathResult                          
    character(512) :: rawDataPath
    character(2000) :: loadDataFile, dimString0, dimString0Vector, dimString1, nxString, nyString, nzString
    integer, parameter :: dx=voxelRes
    integer, parameter :: dy=voxelRes
    integer, parameter :: dz=voxelRes
    integer(kind=bitLength/8), allocatable :: currSample(:)
    
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
    real(kind=4*realPrecision) :: T_P, T_F, T_B, T_E, T_W, T_N, T_S
    real(kind=4*realPrecision) :: currResidual
    real(kind=4*realPrecision) :: totalFluxFront, totalFluxBack, totalFluxEast, totalFluxWest, totalFluxNorth, totalFluxSouth
    real(kind=4*realPrecision) :: aveFluxX, aveFluxY, aveFluxZ, kEffX, kEffY, kEffZ
    integer(kind=1) BC_F, BC_B, BC_E, BC_W, BC_N, BC_S
    character(8) :: patchTypeFront, patchTypeBack, patchTypeEast, patchTypeWest, patchTypeNorth, patchTypeSouth
    real(kind=4*realPrecision), allocatable :: r(:), p(:), xField(:), Apk(:), zVec(:) !variables for conjugate gradient method
    real(kind=4*realPrecision) :: alpha, beta, rDotProductOld, residualFactor
    real(kind=4*realPrecision), allocatable :: qFlux(:) !variables for flux calculation (postprocessing)
    integer :: arrayPos
    
    !initialize path variables
    write(rawDataPath,"(2A)") TRIM(casePath), TRIM(sampleName)
    write(pathTemperature,"(2A)") TRIM(casePath), TRIM('resTemp.raw')
    write(pathFlux,"(2A)") TRIM(casePath), TRIM('resFlux.raw')
    write(pathPost,"(2A)") TRIM(casePath), TRIM('loadData_paraview.xdmf')
    write(pathResult,"(2A)") TRIM(casePath), TRIM('conductivity.csv')


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
    open(10, file=rawDataPath, form='unformatted', access='stream', status='old')
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
    include './src/buildLinearSystem.f90'
    write(*,*) '... done'
    
    !preconditioned conjugate gradient with Jacobi-Diagonal-Preconditioner
    write(*,*) '3) Solving linear system ...'    
    include './src/conjugateGradient.f90'    
    write(*,*) '... done'
    
    write(*,*) 'Final result:'
    open(1, file=pathResult, access='stream', status='replace', form='formatted')
    write(1,*) 'Iteration ', i,' residual (L1-Norm): ', currResidual
    if (evalDirection=='X') then
        write(*,*) 'totalFluxFront= ', totalFluxFront, ', totalFluxBack= ', totalFluxBack, ', kEffX= ', kEffX        
        write(1,*) 'totalFluxFront= ', totalFluxFront, ', totalFluxBack= ', totalFluxBack, ', kEffX= ', kEffX
    else if (evalDirection=='Y') then
        write(*,*) 'totalFluxEast= ', totalFluxEast, ', totalFluxWest= ', totalFluxWest, ', kEffY= ', kEffY
        write(1,*) 'totalFluxEast= ', totalFluxEast, ', totalFluxWest= ', totalFluxWest, ', kEffY= ', kEffY
    else if (evalDirection=='Z') then
        write(*,*) 'totalFluxNorth= ', totalFluxNorth, ', totalFluxSouth= ', totalFluxSouth, ', kEffZ= ', kEffZ
        write(1,*) 'totalFluxNorth= ', totalFluxNorth, ', totalFluxSouth= ', totalFluxSouth, ', kEffZ= ', kEffZ
    end if
    close(1)


    !postprocessing (if activated)
    if (writeTempertureField) then
        write(*,*) 'Writing temperature field ...'
        open(1, file=pathTemperature, access='stream', status='replace', form='unformatted')
        write(1) real(xField,4)
        close(1)
        write(*,*) '... done'
    end if
    
    if (writeFluxField) then
        write(*,*) 'Writing flux field ...'         
        include './src/fluxCalculation.f90'
        open(1, file=pathFlux, access='stream', status='replace', form='unformatted')
        write(1) real(qFlux,4)
        close(1)
        write(*,*) '... done'
    end if    
 
    include './src/writeParaviewFile.f90'    
    open(1, file=pathPost, access='stream', status='replace', form='formatted')
    write(1,'(A)') loadDataFile
    close(1)

    !write(*,*) "Press Enter to exit..."
    !read(*,*)
    
    contains
    include './src/mainFunctions.f90'
end program conductivitySolver