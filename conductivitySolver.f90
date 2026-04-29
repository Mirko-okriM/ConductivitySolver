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
!   Version 1.2 (24.08.25): 
!       - periodic boundary condtions are added
!       - solver computes conductivity tensor
!   Version 1.3 (29.08.25): 
!       - CSR matrix format added
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
    
    integer(kind=4) :: x, y, z, i, j, currCell, nCounter
    integer(kind=4), parameter :: nCells=nx*ny*nz !this is also equal to the number of elements on the main diagonal
    integer(kind=4), allocatable :: colInd(:), rowPtr(:) !arrays for CSR-Format
    integer(kind=4) :: rowPtrEntry
    real(kind=4*realPrecision), allocatable :: A_mat(:), MD(:) !array for CSR-Format     
    !logical :: frontPeriodicBC, backPeriodicBC, eastPeriodicBC, westPeriodicBC, northPeriodicBC, southPeriodicBC !****    
    real(kind=4*realPrecision) :: kF, kB, kE, kW, kN, kS, kP
    real(kind=4*realPrecision) :: aP, aF, aB, aE, aW, aN, aS
    real(kind=4*realPrecision) :: phiP, phiF, phiB, phiE, phiW, phiN, phiS
    real(kind=4*realPrecision) :: currResidual
    real(kind=4*realPrecision) :: xAve, qAveX, qAveY, qAveZ !****
    integer(kind=1) :: Gx, Gy, Gz !****
    real(kind=4*realPrecision) :: totalFluxFront, totalFluxBack, totalFluxEast, totalFluxWest, totalFluxNorth, totalFluxSouth
    real(kind=4*realPrecision) :: aveFluxX, aveFluxY, aveFluxZ, kEffX, kEffY, kEffZ
    real(kind=4*realPrecision), allocatable :: b(:), r(:), p(:), xField(:), Apk(:), zVec(:) !variables for conjugate gradient method
    real(kind=4*realPrecision) :: alpha, beta, rDotProductOld, residualFactor
    real(kind=4*realPrecision), allocatable :: qFlux(:) !variables for flux calculation (postprocessing)
    integer(kind=4), allocatable :: rockVoxelCellNumber(:) !assign a chronologically sorted "cell number" to the rock cells
    integer,parameter :: lengthSortArray=7 !every cell connects at maximum 1(self)+6(faces)=7 cells, self-looping cells have less connections
    integer(kind=4) :: arrayCellNeighbourNumber(lengthSortArray) !auxiloary array to sort values for matrix composition
    real(kind=4*realPrecision) :: arrayCellNeighbourCoeff(lengthSortArray) !auxiloary array to sort values for matrix composition
    
    integer :: arrayPos
    integer :: currRockCell
    integer :: nCellsRock
    integer :: outputInfoCounter
    integer,parameter :: voidCellValue=0
    integer,parameter :: emptyArrayValue=0
    
    integer :: xSearchFront, xSearchBack, ySearchEast, ySearchWest, zSearchNorth, zSearchSouth
    integer :: xSearchFrontOld, xSearchBackOld, ySearchEastOld, ySearchWestOld, zSearchNorthOld, zSearchSouthOld

    
    !initialize path variables
    write(rawDataPath,"(2A)") TRIM(casePath), TRIM(sampleName)
    write(pathTemperature,"(2A)") TRIM(casePath), TRIM('resTemp.raw')
    write(pathFlux,"(2A)") TRIM(casePath), TRIM('resFlux.raw')
    write(pathPost,"(2A)") TRIM(casePath), TRIM('loadData_paraview.xdmf')
    write(pathResult,"(5A)") TRIM(casePath), TRIM(sampleName), TRIM('_conductivity_'), TRIM(evalDirection), TRIM('.csv')

    !openmp might not work if giant arrays are not defined as allocatable variables
    allocate(currSample(nCells))
    allocate(rockVoxelCellNumber(nCells))
    
    !call omp_set_num_threads(nProcessors)        
    
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
    if (searchDirection=='forwardSearch') then
        include './src/buildLinearSystem_forwardSearch.f90'
    else if (searchDirection=='backwardSearch') then   
        include './src/buildLinearSystem_backwardSearch.f90' 
    end if        
    write(*,*) '... done'
    
    !preconditioned conjugate gradient with Jacobi-Diagonal-Preconditioner
    write(*,*) '3) Solving linear system ...'    
    include './src/conjugateGradient.f90'    
    write(*,*) '... done'
    write(*,*) ''
    
    write(*,*) 'Final result:'
    open(1, file=pathResult, access='stream', status='replace', form='formatted')
    write(1,*) 'Iteration ', i,' residual (L1-Norm): ', currResidual
    if (evalDirection=='X') then
        write(*,*) 'k_xx= ', qAveX, new_line('a'), ' k_yx= ', qAveY, new_line('a'), ' k_zx= ', qAveZ        
        write(1,*) 'k_xx= ', qAveX, new_line('a'), ' k_yx= ', qAveY, new_line('a'), ' k_zx= ', qAveZ 
    else if (evalDirection=='Y') then
        write(*,*) 'k_xy= ', qAveX, new_line('a'), ' k_yy= ', qAveY, new_line('a'), ' k_zy= ', qAveZ
        write(1,*) 'k_xy= ', qAveX, new_line('a'), ' k_yy= ', qAveY, new_line('a'), ' k_zy= ', qAveZ
    else if (evalDirection=='Z') then
        write(*,*) 'k_xz= ', qAveX, new_line('a'), ' k_yz= ', qAveY, new_line('a'), ' k_zz= ', qAveZ
        write(1,*) 'k_xz= ', qAveX, new_line('a'), ' k_yz= ', qAveY, new_line('a'), ' k_zz= ', qAveZ
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
        if (searchDirection=='forwardSearch') then
            include './src/fluxCalculation_forwardSearch.f90' 
        else if (searchDirection=='backwardSearch') then   
            include './src/fluxCalculation_backwardSearch.f90' 
        end if  
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
