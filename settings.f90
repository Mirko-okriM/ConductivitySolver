!-----------------------------------------------------------------------------------------------------------------------------------------!
!--------- USER INPUT FILE ---------------------------------------------------------------------------------------------------------------! 
!-----------------------------------------------------------------------------------------------------------------------------------------!
! - compile solver with: gfortran.exe -fopenmp -O3 conductivitySolver.f90 -o solverRun.exe
! - number of cpus/threads has to be adjusted with environment variable OMP_NUM_THREADS, e.g. on windows 11 --> $env:OMP_NUM_THREADS=8     

integer, parameter :: nx = 400                                                                                                 !Image width
integer, parameter :: ny = 400                                                                                                  !Image height
integer, parameter :: nz = 400                                                                                                  !Image depth
integer, parameter :: voxelRes=1 !do not change! edge length of a cell dx=dy=dz                                                 !VoxelResolution of input-raw data --> SHOULD STAY 1, VALUE DOES NOT REALLY HAVE INFLUENCE ON RESULT, MIGHT BE CHANGED TO TYPE REAL IN FUTURE VERSION
integer, parameter :: realPrecision=2                                                                                           !floating-point precision: 1 (single) or 2 (double)
integer, parameter :: nIteration=100000                                                                                         !max. number of iterations        
integer, parameter :: calcFluxStep=30                                                                                           !compute flux each "calcFluxStep-iteration" (flux calcuation at every step is expansive) 
real(kind=4*realPrecision), parameter :: targetResidual=1e-10                                                                   !target residual    
integer(kind=1), parameter :: bitLength = 8                                                                                     !Bites per pixel, typically grey-values in raw images are stored as either 8-bit or 16-bit integers
integer, dimension(2) :: rawValues=[25,255]                                                                                      !bit-grayvalues of each phase (unsigned format)
real(kind=4*realPrecision), dimension(2) :: correspConductivity=[1.0,0.1]                                                  !conductivity of each phase (only double value, e.g. 1.0)
integer, parameter :: emptyPhaseRawValue = 0 !v2 USE POSITIVE INTEGER
character(512) :: casePath = 'C:\Users\Mirko\Desktop\solver\solver_conductivity\paperNew\'      !END STRING WITH SLASH             !path to case folder (place raw-data in this folder, results will as be saved there as well)
character(512) :: sampleName = 'Berea_400cube_oneSphere.raw'                                                                    !input file (veldsteen_poreus_250cube, Berea_2d25um_binary_400cube, spheres_150cube)
character(1) :: evalDirection = 'Z'                                                                                             !direction of calculation (X or Y or Z)
character(512) :: searchDirection = 'forwardSearch' !forwardSearch or backwardSearch
logical :: writeTempertureField = .False.                                                                                         !write calculated temperature-field (.TRUE. or .FALSE.)
logical :: writeFluxField = .False.                                                                                               !write calculated flux-field (.TRUE. or .FALSE.)
        
!----------------------------------------------------------------------------------------------------------------------------------------!
!--------- If errors occur, please feel free to contact me. I will be happy to help: mirko.siegert@web.de -------------------------------!
!----------------------------------------------------------------------------------------------------------------------------------------!
