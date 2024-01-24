!determination of matrix elements for each cell
!$omp parallel do collapse(3) private(kF, kB, kE, kW, kN, kS, aP, aF, aB, aE, aW, aN, aS, bCell, currCell)
    do z=1,nz
        do y=1,ny
            do x=1,nx
                currCell=cellPos(x,y,z)
                if (isInnerCell(x,y,z)) then !hardcoded equation result
                    kF=kEval(x,y,z,x+1,y,z)
                    kB=kEval(x,y,z,x-1,y,z)
                    kE=kEval(x,y,z,x,y+1,z)
                    kW=kEval(x,y,z,x,y-1,z)
                    kN=kEval(x,y,z,x,y,z+1)
                    kS=kEval(x,y,z,x,y,z-1)
                    aP=-(kB+kW+kS+kF+kE+kN)
                    MD(currCell)=aP
                    b(currCell)=0
                    ODp1(currCell)=kF
                    !ODm1(currCell+posODm1)=kB; !not necessary since matrix is symmetric
                    ODp2(currCell)=kE
                    !ODm2(currCell+posODm2)=kW; !not necessary since matrix is symmetric
                    ODp3(currCell)=kN
                    !ODm3(currCell+posODm3)=kS; !not necessary since matrix is symmetric
                else !boundaryCell (code of this part might no be very elegant)
                    aP=0; aF=0; aB=0; aE=0; aW=0; aN=0; aS=0;
                    bCell=0
                    !evalFront
                    if (x<nx) then !NB is cell
                        kF=kEval(x,y,z,x+1,y,z)
                        aP=aP-kF
                        aF=kF
                    else !NB is boundary
                        !aP=aP+0 !default case (not necessary to evaluate)
                        !aF=0    !default case (not necessary to evaluate)
                        if (patchTypeFront=='patch') then
                            kF=kLookUp(x,y,z) !kF is equal to kP since cell is at boundary
                            aP=aP-2*kF
                            !aF=0                 !default case (not necessary to evaluate)
                            bCell=-2*kF*BC_F
                        end if
                    end if
                    !evalBack
                    if (x>1) then !NB is cell
                        kB=kEval(x,y,z,x-1,y,z)
                        aP=aP-kB
                        aB=kB
                    else !NB is boundary
                        !aP=aP+0
                        !aB=0
                        if (patchTypeBack=="patch") then
                            kB=kLookUp(x,y,z) !kB is equal to kP since cell is at boundary
                            aP=aP-2*kB
                            !aB=0                !default case (not necessary to evaluate)
                            bCell=-2*kB*BC_B
                        end if
                    end if
                    !evalEast
                    if (y<ny) then !NB is cell
                        kE=kEval(x,y,z,x,y+1,z)
                        aP=aP-kE
                        aE=kE
                    else !NB is boundary
                        !aP=aP+0
                        !aE=0
                        if (patchTypeEast=="patch") then
                            kE=kLookUp(x,y,z) !kE is equal to kP since cell is at boundary
                            aP=aP-2*kE
                            !aE=0                 !default case (not necessary to evaluate)
                            bCell=-2*kE*BC_E
                        end if
                    end if
                    !evalWest
                    if (y>1) then !NB is cell
                        kW=kEval(x,y,z,x,y-1,z)
                        aP=aP-kW
                        aW=kW
                    else !NB is boundary
                        !aP=aP+0
                        !aW=0
                        if (patchTypeWest=="patch") then
                            kW=kLookUp(x,y,z) !kN is equal to kP since cell is at boundary
                            aP=aP-2*kW
                            !aW=0                 !default case (not necessary to evaluate)
                            bCell=-2*kW*BC_W
                        end if
                    end if
                    !evalNorth
                    if (z<nz) then !NB is cell
                        kN=kEval(x,y,z,x,y,z+1)
                        aP=aP-kN
                        aN=kN
                    else !NB is boundary
                        !aP=aP+0
                        !aN=0
                        if (patchTypeNorth=="patch") then
                            kN=kLookUp(x,y,z) !kN is equal to kP since cell is at boundary
                            aP=aP-2*kN
                            !aN=0                !default case (not necessary to evaluate)
                            bCell=-2*kN*BC_N
                        end if
                    end if
                    !evalSouth
                    if (z>1) then !NB is cell
                        kS=kEval(x,y,z,x,y,z-1)
                        aP=aP-kS
                        aS=kS
                    else !NB is boundary
                        !aP=aP+0
                        !aS=0
                        if (patchTypeSouth=="patch") then
                            kS=kLookUp(x,y,z) !kS is equal to kP since cell is at boundary
                            aP=aP-2*kS
                            !aS=0                 !default case (not necessary to evaluate)
                            bCell=-2*kS*BC_S
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