!fluxCalculation
qAveX=0.0d0
qAveY=0.0d0
qAveZ=0.0d0
!deallocate(MD,ODp1,ODp2,ODp3,r,p,Apk,zVec)    
!$omp parallel do collapse(3) reduction(+:qAveX, qAveY, qAveZ) &
!$omp& private(kF, kB, kE, kW, kN, kS, kP, phiP, phiF, phiB, phiE, phiW, phiN, phiS, currCell, currRockCell, arrayPos, &
!$omp&                xSearchFront, xSearchBack, ySearchEast, ySearchWest, zSearchNorth, zSearchSouth, x, y, z) 
do z=1,nz
    do y=1,ny
        do x=1,nx
            currCell=cellPos(x,y,z)
            currRockCell=rockVoxelCellNumber(currCell)
            if (currRockCell .NE. voidCellValue) then
                arrayPos=1+(currRockCell-1)*3
                phiP=xField(currRockCell)
                kP=kLookUp(x,y,z)
                !- go through each cell ...
                !... and check if front/back/east/west/north/south-face of current cell belongs to neighbour-cell or boundary
                !- if face belongs to boundary, then implement periodic boundary condition (facenormal-backward-search)
                aP=0; aF=0; aB=0; aE=0; aW=0; aN=0; aS=0; !coefficient initilization
            !evalFront
                !Step 1: Check if direct neighbour is rockCell
                xSearchFront=x+1 
                if (x == nx) then !if x was already last position, use loop around logic for xSearchFront
                    xSearchFront=1
                end if
                if (rockVoxelCellNumber(cellPos(xSearchFront,y,z)) .NE. voidCellValue) then
                    kF=kEval(x,y,z,xSearchFront,y,z)
                    phiF=xField(rockVoxelCellNumber(cellPos(xSearchFront,y,z)))                     
                !Step 2: Apply backwardSearch logic, if direct neighbour is voidCell
                else
                    xSearchFrontOld=x
                    xSearchFront=x-1 !generally, xSearchFront is the previous x-value
                    if (x == 1) then !if x was already last position, use loop around logic for xSearchFront
                        xSearchFront=nx 
                    end if
                    !find neighbour with facenormal-backward-search                      
                    do i=1,nx
                        !check if currSearchCell is rockCell: if yes, connect both faces
                        if (rockVoxelCellNumber(cellPos(xSearchFront,y,z)) == voidCellValue) then
                            xSearchFront=xSearchFrontOld
                            kF=kEval(x,y,z,xSearchFront,y,z) 
                            phiF=xField(rockVoxelCellNumber(cellPos(xSearchFront,y,z))) 
                            ! if (xSearchFront==x+1) then !set boundary condition flag
                                ! frontPeriodicBC=.False. 
                            ! else
                                ! frontPeriodicBC=.True.
                            ! end if
                            EXIT
                        else 
                            xSearchFrontOld=xSearchFront
                            xSearchFront=xSearchFront-1
                            if (xSearchFront < 1) then !loop around condition
                                xSearchFront=nx
                            end if
                        end if
                    end do 
                end if
            !evalBack
                !Step 1: Check if direct neighbour is rockCell
                xSearchBack=x-1 
                if (x == 1) then !if x was already last position, use loop around logic for xSearchBack
                    xSearchBack=nx 
                end if
                if (rockVoxelCellNumber(cellPos(xSearchBack,y,z)) .NE. voidCellValue) then
                    kB=kEval(x,y,z,xSearchBack,y,z)
                    phiB=xField(rockVoxelCellNumber(cellPos(xSearchBack,y,z)))
                !Step 2: Apply backwardSearch logic, if direct neighbour is voidCell
                else
                    xSearchBackOld=x
                    xSearchBack=x+1 !generally, xSearchBack is the next x-value
                    if (x == nx) then !if x was already first position, use loop around logic for xSearchBack
                        xSearchBack=1 
                    end if
                    !find neighbour with facenormal-backward-search                      
                    do i=1,nx
                        !check if currSearchCell is rockCell: if yes, connect both faces
                        if (rockVoxelCellNumber(cellPos(xSearchBack,y,z)) == voidCellValue) then
                            xSearchBack=xSearchBackOld
                            kB=kEval(x,y,z,xSearchBack,y,z)
                            phiB=xField(rockVoxelCellNumber(cellPos(xSearchBack,y,z)))                            
                            EXIT
                        else 
                            xSearchBackOld=xSearchBack
                            xSearchBack=xSearchBack+1
                            if (xSearchBack > nx) then !loop around condition
                                xSearchBack=1
                            end if
                        end if
                    end do
                end if
            !evalEast
                !Step 1: Check if direct neighbour is rockCell
                ySearchEast=y+1 
                if (y == ny) then 
                    ySearchEast=1 
                end if
                if (rockVoxelCellNumber(cellPos(x,ySearchEast,z)) .NE. voidCellValue) then
                    kE=kEval(x,y,z,x,ySearchEast,z)
                    phiE=xField(rockVoxelCellNumber(cellPos(x,ySearchEast,z)))
                !Step 2: Apply backwardSearch logic, if direct neighbour is voidCell
                else
                    ySearchEastOld=y
                    ySearchEast=y-1 !generally, ySearchEast is the previous y-value
                    if (y == 1) then !if y was already last position, use loop around logic for ySearchEast
                        ySearchEast=ny 
                    end if
                    !find neighbour with facenormal-backward-search                      
                    do i=1,ny
                        !check if currSearchCell is rockCell: if yes, connect both faces
                        if (rockVoxelCellNumber(cellPos(x,ySearchEast,z)) == voidCellValue) then
                            ySearchEast=ySearchEastOld
                            kE=kEval(x,y,z,x,ySearchEast,z)
                            phiE=xField(rockVoxelCellNumber(cellPos(x,ySearchEast,z)))                            
                            EXIT
                        else 
                            ySearchEastOld=ySearchEast
                            ySearchEast=ySearchEast-1
                            if (ySearchEast < 1) then !loop around condition
                                ySearchEast=ny
                            end if
                        end if
                    end do
                end if
            !evalWest
                !Step 1: Check if direct neighbour is rockCell
                ySearchWest=y-1 
                if (y == 1) then 
                    ySearchWest=ny 
                end if
                if (rockVoxelCellNumber(cellPos(x,ySearchWest,z)) .NE. voidCellValue) then
                    kW=kEval(x,y,z,x,ySearchWest,z)
                    phiW=xField(rockVoxelCellNumber(cellPos(x,ySearchWest,z)))
                !Step 2: Apply backwardSearch logic, if direct neighbour is voidCell
                else
                    ySearchWestOld=y
                    ySearchWest=y+1 !generally, ySearchWest is the next y-value
                    if (y == ny) then !if y was already first position, use loop around logic for ySearchWest
                        ySearchWest=1
                    end if
                    !find neighbour with facenormal-backward-search                      
                    do i=1,ny
                        !check if currSearchCell is rockCell: if yes, connect both faces
                        if (rockVoxelCellNumber(cellPos(x,ySearchWest,z)) == voidCellValue) then
                            ySearchWest=ySearchWestOld
                            kW=kEval(x,y,z,x,ySearchWest,z)
                            phiW=xField(rockVoxelCellNumber(cellPos(x,ySearchWest,z)))
                            EXIT
                        else 
                            ySearchWestOld=ySearchWest
                            ySearchWest=ySearchWest+1
                            if (ySearchWest > ny) then !loop around condition
                                ySearchWest=1
                            end if
                        end if
                    end do
                end if
            !evalNorth
                !Step 1: Check if direct neighbour is rockCell
                zSearchNorth=z+1 
                if (z == nz) then 
                    zSearchNorth=1 
                end if
                if (rockVoxelCellNumber(cellPos(x,y,zSearchNorth)) .NE. voidCellValue) then
                    kN=kEval(x,y,z,x,y,zSearchNorth) 
                    phiN=xField(rockVoxelCellNumber(cellPos(x,y,zSearchNorth)))
                !Step 2: Apply backwardSearch logic, if direct neighbour is voidCell
                else
                    zSearchNorthOld=z
                    zSearchNorth=z-1 !generally, zSearchNorth is the previous z-value
                    if (z == 1) then !if z was already last position, use loop around logic for zSearchNorth
                        zSearchNorth=nz 
                    end if
                    !find neighbour with facenormal-backward-search                      
                    do i=1,nz
                        !check if currSearchCell is rockCell: if yes, connect both faces
                        if (rockVoxelCellNumber(cellPos(x,y,zSearchNorth)) == voidCellValue) then
                            zSearchNorth=zSearchNorthOld
                            kN=kEval(x,y,z,x,y,zSearchNorth) 
                            phiN=xField(rockVoxelCellNumber(cellPos(x,y,zSearchNorth)))
                            EXIT
                        else 
                            zSearchNorthOld=zSearchNorth
                            zSearchNorth=zSearchNorth-1
                            if (zSearchNorth < 1) then !loop around condition
                                zSearchNorth=nz
                            end if
                        end if
                    end do    
                end if                    
            !evalSouth
                !Step 1: Check if direct neighbour is rockCell
                zSearchSouth=z-1 
                if (z == 1) then 
                    zSearchSouth=nz 
                end if
                if (rockVoxelCellNumber(cellPos(x,y,zSearchSouth)) .NE. voidCellValue) then
                    kS=kEval(x,y,z,x,y,zSearchSouth) 
                    phiS=xField(rockVoxelCellNumber(cellPos(x,y,zSearchSouth))) 
                !Step 2: Apply backwardSearch logic, if direct neighbour is voidCell
                else
                    zSearchSouthOld=z
                    zSearchSouth=z+1 !generally, zSearchSouth is the next z-value
                    if (z == nz) then !if z was already first position, use loop around logic for zSearchSouth
                        zSearchSouth=1 
                    end if
                    !find neighbour with facenormal-backward-search                      
                    do i=1,nz
                        !check if currSearchCell is rockCell: if yes, connect both faces
                        if (rockVoxelCellNumber(cellPos(x,y,zSearchSouth)) == voidCellValue) then
                            zSearchSouth=zSearchSouthOld
                            kS=kEval(x,y,z,x,y,zSearchSouth)  
                            phiS=xField(rockVoxelCellNumber(cellPos(x,y,zSearchSouth))) 
                            EXIT
                        else 
                            zSearchSouthOld=zSearchSouth
                            zSearchSouth=zSearchSouth+1
                            if (zSearchSouth > nz) then !loop around condition
                                zSearchSouth=1
                            end if
                        end if
                    end do
                end if 
            !calc qFlux components
                qFlux(arrayPos)=(-kF*((phiF-phiP)/dx+Gx)+kB*((phiB-phiP)/dx-Gx))/2         !x-component of fluxVector for currCell
                qAveX=qAveX-qFlux(arrayPos)
                qFlux(arrayPos+1)=(-kE*((phiE-phiP)/dy+Gy)+kW*((phiW-phiP)/dy-Gy))/2      !y-component of fluxVector for currCell
                qAveY=qAveY-qFlux(arrayPos+1)
                qFlux(arrayPos+2)=(-kN*((phiN-phiP)/dz+Gz)+kS*((phiS-phiP)/dz-Gz))/2      !z-component of fluxVector for currCell 
                qAveZ=qAveZ-qFlux(arrayPos+2)    
            end if                
        end do
    end do
end do
qAveX=qAveX/nCellsRock
qAveY=qAveY/nCellsRock
qAveZ=qAveZ/nCellsRock
