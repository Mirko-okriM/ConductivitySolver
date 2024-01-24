!fluxCalculation
deallocate(MD,ODp1,ODp2,ODp3,r,p,Apk,zVec)    
allocate(qFlux(3*nCells))
!$omp parallel do collapse(3) private(kF, kB, kE, kW, kN, kS, T_P, T_F, T_B, T_E, T_W, T_N, T_S, currCell, arrayPos)
do z=1,nz
    do y=1,ny
        do x=1,nx
            currCell=cellPos(x,y,z)
            arrayPos=1+(currCell-1)*3
            T_P=xField(currCell)
            if (isInnerCell(x,y,z)) then !hardcoded equation result
                    kF=kEval(x,y,z,x+1,y,z)
                    kB=kEval(x,y,z,x-1,y,z)
                    kE=kEval(x,y,z,x,y+1,z)
                    kW=kEval(x,y,z,x,y-1,z)
                    kN=kEval(x,y,z,x,y,z+1)
                    kS=kEval(x,y,z,x,y,z-1)
                    T_F=xField(cellPos(x+1,y,z))
                    T_B=xField(cellPos(x-1,y,z))
                    T_E=xField(cellPos(x,y+1,z))
                    T_W=xField(cellPos(x,y-1,z))
                    T_N=xField(cellPos(x,y,z+1))
                    T_S=xField(cellPos(x,y,z-1))
                    qFlux(arrayPos)=-1/dx*(kF*(T_F-T_P)+kB*(T_P-T_B))/2         !x-component of fluxVector for currCell
                    qFlux(arrayPos+1)=-1/dy*(kE*(T_E-T_P)+kW*(T_P-T_W))/2       !y-component of fluxVector for currCell
                    qFlux(arrayPos+2)=-1/dz*(kN*(T_N-T_P)+kS*(T_P-T_S))/2       !z-component of fluxVector for currCell
            else !(placeholder)   
                !evalFront
                qFlux(arrayPos)=0.0
                qFlux(arrayPos+1)=0.0
                qFlux(arrayPos+2)=0.0
                if (x<nx) then !NB is cell
                    kF=kEval(x,y,z,x+1,y,z)
                    T_F=xField(cellPos(x+1,y,z))
                    qFlux(arrayPos)=qFlux(arrayPos)-kF*(T_F-T_P)/dx/2
                else !NB is boundary                            
                    !qFlux(arrayPos)=qFlux(arrayPos)+0.0        !isWall
                    if (patchTypeFront=='patch') then
                        kF=kLookUp(x,y,z)
                        qFlux(arrayPos)=qFlux(arrayPos)-kF*(BC_F-T_P)/dx
                    end if
                end if
                !evalBack
                if (x>1) then !NB is cell
                    kB=kEval(x,y,z,x-1,y,z)
                    T_B=xField(cellPos(x-1,y,z))
                    qFlux(arrayPos)=qFlux(arrayPos)-kB*(T_P-T_B)/dx/2
                else !NB is boundary
                    !qFlux(arrayPos)=qFlux(arrayPos)+0.0        !isWall
                    if (patchTypeBack=="patch") then
                        kB=kLookUp(x,y,z)
                        qFlux(arrayPos)=qFlux(arrayPos)-kB*(T_P-BC_B)/dx
                    end if
                end if
                !evalEast
                if (y<ny) then !NB is cell
                    kE=kEval(x,y,z,x,y+1,z)
                    T_E=xField(cellPos(x,y+1,z))
                    qFlux(arrayPos+1)=qFlux(arrayPos+1)-kE*(T_E-T_P)/dy/2
                else !NB is boundary
                    !qFlux(arrayPos)=qFlux(arrayPos)+0.0        !isWall
                    if (patchTypeEast=="patch") then
                        kE=kLookUp(x,y,z)
                        qFlux(arrayPos+1)=qFlux(arrayPos+1)-kE*(BC_E-T_P)/dy
                    end if
                end if
                !evalWest
                if (y>1) then !NB is cell
                    kW=kEval(x,y,z,x,y-1,z)
                    T_W=xField(cellPos(x,y-1,z))
                    qFlux(arrayPos+1)=qFlux(arrayPos+1)-kW*(T_P-T_W)/dy/2
                else !NB is boundary
                    !qFlux(arrayPos)=qFlux(arrayPos)+0.0        !isWall
                    if (patchTypeWest=="patch") then
                        kW=kLookUp(x,y,z)
                        qFlux(arrayPos+1)=qFlux(arrayPos+1)-kW*(T_P-BC_W)/dy
                    end if
                end if
                !evalNorth
                if (z<nz) then !NB is cell
                    kN=kEval(x,y,z,x,y,z+1)
                    T_N=xField(cellPos(x,y,z+1))
                    qFlux(arrayPos+2)=qFlux(arrayPos+2)-kN*(T_N-T_P)/dz/2
                else !NB is boundary
                    !qFlux(arrayPos)=qFlux(arrayPos)+0.0        !isWall
                    if (patchTypeNorth=="patch") then
                        kN=kLookUp(x,y,z)
                        qFlux(arrayPos+2)=qFlux(arrayPos+2)-kN*(BC_N-T_P)/dz
                    end if
                end if
                !evalSouth
                if (z>1) then !NB is cell
                    kS=kEval(x,y,z,x,y,z-1)
                    T_S=xField(cellPos(x,y,z-1))
                    qFlux(arrayPos+2)=qFlux(arrayPos+2)-kS*(T_P-T_S)/dz/2
                else !NB is boundary
                    !qFlux(arrayPos)=qFlux(arrayPos)+0.0        !isWall
                    if (patchTypeSouth=="patch") then
                        kS=kLookUp(x,y,z)
                        qFlux(arrayPos+2)=qFlux(arrayPos+2)-kS*(T_P-BC_S)/dz
                    end if
                end if                            
            end if                
        end do
    end do
end do