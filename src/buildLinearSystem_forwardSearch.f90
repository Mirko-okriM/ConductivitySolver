!determination of matrix elements for each cell
!periodic boundary conditions are used
include './src/fillCellNumberRockArray.f90' !assign a chronologically sorted "cell number" to the rock cells   
nCounter=0
A_mat=0.0d0
colInd=0
rowPtrEntry=1 !start with 1, so first value that is calculated in the following logic is stored at array position 2
rowPtr=1 !first entry has to be 1, so 1 is used for initilization
do z=1,nz    
    do y=1,ny
        do x=1,nx
            currCell=cellPos(x,y,z)
            currRockCell=rockVoxelCellNumber(currCell)
            if (currRockCell .NE. voidCellValue) then
                !- go through each cell ...
                !... and check if front/back/east/west/north/south-face of current cell belongs to neighbour-cell or boundary
                !- if face belongs to boundary, then implement periodic boundary condition (facenormal-forward-search)
                aP=0; aF=0; aB=0; aE=0; aW=0; aN=0; aS=0; !coefficient initilization
            !evalFront
                xSearchFront=x+1 !generally, xSearchFront is the next x-value
                if (x == nx) then !if x was already last position, use loop around logic for xSearchFront
                    xSearchFront=1 
                end if
                !find neighbour with facenormal-forward-search                      
                do i=1,nx
                    !check if currSearchCell is rockCell: if yes, connect both faces
                    if (rockVoxelCellNumber(cellPos(xSearchFront,y,z)) .NE. voidCellValue) then
                        kF=kEval(x,y,z,xSearchFront,y,z) 
                        ! if (xSearchFront==x+1) then !set boundary condition flag
                            ! frontPeriodicBC=.False. 
                        ! else
                            ! frontPeriodicBC=.True.
                        ! end if
                        EXIT
                    else 
                        xSearchFront=xSearchFront+1
                        if (xSearchFront > nx) then !loop around condition
                            xSearchFront=1
                        end if
                    end if
                end do 
            !evalBack
                xSearchBack=x-1 !generally, xSearchBack is previous x-value
                if (x == 1) then !if x was already first position, use loop around logic for xSearchBack
                    xSearchBack=nx 
                end if
                !find neighbour with facenormal-forward-search                      
                do i=1,nx
                    !check if currSearchCell is rockCell: if yes, connect both faces
                    if (rockVoxelCellNumber(cellPos(xSearchBack,y,z)) .NE. voidCellValue) then
                        kB=kEval(x,y,z,xSearchBack,y,z) 
                        EXIT
                    else 
                        xSearchBack=xSearchBack-1
                        if (xSearchBack < 1) then !loop around condition
                            xSearchBack=nx
                        end if
                    end if
                end do
            !evalEast
                ySearchEast=y+1 !generally, ySearchEast is the next y-value
                if (y == ny) then !if y was already last position, use loop around logic for ySearchEast
                    ySearchEast=1 
                end if
                !find neighbour with facenormal-forward-search                      
                do i=1,ny
                    !check if currSearchCell is rockCell: if yes, connect both faces
                    if (rockVoxelCellNumber(cellPos(x,ySearchEast,z)) .NE. voidCellValue) then
                        kE=kEval(x,y,z,x,ySearchEast,z) 
                        EXIT
                    else 
                        ySearchEast=ySearchEast+1
                        if (ySearchEast > ny) then !loop around condition
                            ySearchEast=1
                        end if
                    end if
                end do
            !evalWest
                ySearchWest=y-1 !generally, ySearchWest is previous y-value
                if (y == 1) then !if y was already first position, use loop around logic for ySearchWest
                    ySearchWest=ny
                end if
                !find neighbour with facenormal-forward-search                      
                do i=1,ny
                    !check if currSearchCell is rockCell: if yes, connect both faces
                    if (rockVoxelCellNumber(cellPos(x,ySearchWest,z)) .NE. voidCellValue) then
                        kW=kEval(x,y,z,x,ySearchWest,z)
                        EXIT
                    else 
                        ySearchWest=ySearchWest-1
                        if (ySearchWest < 1) then !loop around condition
                            ySearchWest=ny
                        end if
                    end if
                end do
            !evalNorth
                zSearchNorth=z+1 !generally, zSearchNorth is the next z-value
                if (z == nz) then !if z was already last position, use loop around logic for zSearchNorth
                    zSearchNorth=1 
                end if
                !find neighbour with facenormal-forward-search                      
                do i=1,nz
                    !check if currSearchCell is rockCell: if yes, connect both faces
                    if (rockVoxelCellNumber(cellPos(x,y,zSearchNorth)) .NE. voidCellValue) then
                        kN=kEval(x,y,z,x,y,zSearchNorth) 
                        EXIT
                    else 
                        zSearchNorth=zSearchNorth+1
                        if (zSearchNorth > nz) then !loop around condition
                            zSearchNorth=1
                        end if
                    end if
                end do                
            !evalSouth
                zSearchSouth=z-1 !generally, zSearchSouth is previous z-value
                if (z == 1) then !if z was already first position, use loop around logic for zSearchSouth
                    zSearchSouth=nz 
                end if
                !find neighbour with facenormal-forward-search                      
                do i=1,nz
                    !check if currSearchCell is rockCell: if yes, connect both faces
                    if (rockVoxelCellNumber(cellPos(x,y,zSearchSouth)) .NE. voidCellValue) then
                        kS=kEval(x,y,z,x,y,zSearchSouth)  
                        EXIT
                    else 
                        zSearchSouth=zSearchSouth-1
                        if (zSearchSouth < 1) then !loop around condition
                            zSearchSouth=nz
                        end if
                    end if
                end do
            !assign k-coefficients to matrix-coefficients (basically, if i stick to the k-variable-names, this step is not necessary ...    
            !... logic for aP is necessary)            
                aF=kF 
                aB=kB
                aE=kE
                aW=kW
                aN=kN
                aS=kS      
                !self-looping cells have to be considered, i.e. cells that share a periodic boundary with itself
                !remark "self-looping cells": the corresponding coefficient has no impact to the main- and the off-diagonal entries
                !only non-self-looping cells are influenced by all coefficients (see logic below and handwritten notes)
                !for a given direction, self-looping cells exist, if the following condition is met: xSearchFront==xSearchBack==x, ...
                !futhermore, if a self-looping cell exists, the number of A_mat entries is reduced by one for each self-looping direction                
                if ((xSearchFront .NE. x) .AND. (xSearchBack .NE. x)) aP=aP-(kF+kB)  !to be honest, the .AND.-logic is not necessary 
                if ((ySearchEast .NE. y) .AND. (ySearchWest .NE. y)) aP=aP-(kE+kW)
                if ((zSearchNorth .NE. z) .AND. (zSearchSouth .NE. z)) aP=aP-(kN+kS)   
                if (aP==0) then
                    write(*,*) 'ERROR: Cell (', x, ', ', y, ', ', z, ') has no face connections to other cells.'
                    STOP
                end if
            !store coefficients to csr-format (logic explained in handwritten notes)                  
                !prepare arrays for sorting step and perform sorting                 
                arrayCellNeighbourNumber=emptyArrayValue !initialize arrays and assign placeholder values
                arrayCellNeighbourCoeff=0.0d0                            
                if ((xSearchFront .NE. x) .AND. (xSearchBack .NE. x)) then !include only non-self-looping directions (see handwritten notes)   
                    if(xSearchFront == xSearchBack) then !logic necessary for the case "only two voxels in one direction"
                        arrayCellNeighbourNumber(1)=rockVoxelCellNumber(cellPos(xSearchFront,y,z))
                        arrayCellNeighbourCoeff(1)=aF+aB
                    else                    
                        arrayCellNeighbourNumber(1:2)=[rockVoxelCellNumber(cellPos(xSearchFront,y,z)), &
                                                        rockVoxelCellNumber(cellPos(xSearchBack,y,z))]
                        arrayCellNeighbourCoeff(1:2)=[aF, aB]
                    end if
                endif
                if ((ySearchEast .NE. y) .AND. (ySearchWest .NE. y)) then !include only non-self-looping directions (see handwritten notes)   
                    if(ySearchEast == ySearchWest) then !logic necessary for the case "only two voxels in one direction"
                        arrayCellNeighbourNumber(3)=rockVoxelCellNumber(cellPos(x,ySearchEast,z))
                        arrayCellNeighbourCoeff(3)=aE+aW
                    else 
                        arrayCellNeighbourNumber(3:4)=[rockVoxelCellNumber(cellPos(x,ySearchEast,z)), &
                                                        rockVoxelCellNumber(cellPos(x,ySearchWest,z))]
                        arrayCellNeighbourCoeff(3:4)=[aE, aW]
                    end if
                endif
                if ((zSearchNorth .NE. z) .AND. (zSearchSouth .NE. z)) then !include only non-self-looping directions (see handwritten notes)   
                    if(zSearchNorth == zSearchSouth) then !logic necessary for the case "only two voxels in one direction"
                        arrayCellNeighbourNumber(5)=rockVoxelCellNumber(cellPos(x,y,zSearchNorth))
                        arrayCellNeighbourCoeff(5)=aN+aS
                    else
                        arrayCellNeighbourNumber(5:6)=[rockVoxelCellNumber(cellPos(x,y,zSearchNorth)), &
                                                        rockVoxelCellNumber(cellPos(x,y,zSearchSouth))]
                        arrayCellNeighbourCoeff(5:6)=[aN, aS]
                    end if
                endif                
                arrayCellNeighbourNumber(7)=currRockCell !include information of current cell
                arrayCellNeighbourCoeff(7)=aP          
                !sort 
                call bubbleSort(arrayCellNeighbourNumber, arrayCellNeighbourCoeff)  !sorts both arrays according two first array       
                !assign values to csr-format
                do i=1,lengthSortArray
                    if (arrayCellNeighbourNumber(i) .NE. emptyArrayValue) then !emptyArrayValues exist for self-looping cells
                        nCounter=nCounter+1
                        !write(*,*) 'nCounter=', nCounter
                        A_mat(nCounter)=arrayCellNeighbourCoeff(i)
                        colInd(nCounter)=arrayCellNeighbourNumber(i)                        
                    end if
                    if (i==lengthSortArray) then
                        rowPtrEntry=rowPtrEntry+1
                        rowPtr(rowPtrEntry)=nCounter+1
                    end if
                end do              
                !save main diagonal for jacobi diagonal preconditioner                
                MD(currRockCell)=aP 
            end if
        !RHS treatment
            if (evalDirection=='X') then !logic explained in handwritten notes
                b(currRockCell)=dx*(kB-kF)
            else if (evalDirection=='Y') then
                b(currRockCell)=dy*(kW-kE)
            else if (evalDirection=='Z') then
                b(currRockCell)=dz*(kS-kN)
            end if                 
        end do
    end do
end do

! do i=1,size(A_mat)
    ! write(*,*) 'A_mat(', i, ')=', A_mat(i)
! end do
! do i=1,size(colInd)
    ! write(*,*) 'colInd(', i, ')=', colInd(i)
! end do
! do i=1,size(rowPtr)
    ! write(*,*) 'rowPtr(', i, ')=', rowPtr(i)
! end do
! do i=1,size(b)
    ! write(*,*) 'b(', i, ')=', b(i)
! end do
! write(*,*) 'nCounter=', nCounter

