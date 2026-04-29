!assigns a chronologically sorted "cell number value" to the rock cells
!iterates through all cells: 
!if cell doesn't belong to rock (i.e. phaseValue=emptyPhaseValue) rockVoxelCellNumber(currCell) gets voidCellValue (for example zero)
!else (i.e. phaseValue!=emptyPhaseValue) rockVoxelCellNumber(currCell) gets nCellsRock
!first rockCell gets value one
!logic see figure "fig. rockCellNumber"
nCellsRock=0
do z=1,nz
    do y=1,ny
        do x=1,nx
            currCell=cellPos(x,y,z)
            if (currSample(currCell) == emptyPhaseRawValue) then
                rockVoxelCellNumber(currCell)=voidCellValue
            else
                nCellsRock=nCellsRock+1
                rockVoxelCellNumber(currCell)=nCellsRock                
            end if
        end do
    end do
end do

!info
write(*,*) 'Found ', nCellsRock, ' rock cells.'

!allocate variables 
allocate(b(nCellsRock)) 
allocate(r(nCellsRock))
allocate(p(nCellsRock))
allocate(xField(nCellsRock))
allocate(Apk(nCellsRock))
allocate(zVec(nCellsRock))
allocate(qFlux(3*nCellsRock))
allocate(colInd(nCellsRock*7)) !array might be to long if self-lopper-cells are present... solver deals with this problem automatically
allocate(rowPtr(nCellsRock+1))
allocate(A_mat(nCellsRock*7)) !array might be to long if self-lopper-cells are present... solver deals with this problem automatically
allocate(MD(nCellsRock))