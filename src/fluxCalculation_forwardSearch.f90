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
                        phiF=xField(rockVoxelCellNumber(cellPos(xSearchFront,y,z)))                      
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
                        phiB=xField(rockVoxelCellNumber(cellPos(xSearchBack,y,z)))                        
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
                        phiE=xField(rockVoxelCellNumber(cellPos(x,ySearchEast,z)))
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
                        phiW=xField(rockVoxelCellNumber(cellPos(x,ySearchWest,z)))
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
                        phiN=xField(rockVoxelCellNumber(cellPos(x,y,zSearchNorth)))
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
                        phiS=xField(rockVoxelCellNumber(cellPos(x,y,zSearchSouth)))                        
                        EXIT
                    else 
                        zSearchSouth=zSearchSouth-1
                        if (zSearchSouth < 1) then !loop around condition
                            zSearchSouth=nz
                        end if
                    end if
                end do
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
