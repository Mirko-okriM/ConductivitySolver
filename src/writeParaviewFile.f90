!strings that are problem specific
write(nxString,*) '   <Attribute Name="nx" Dimensions="', nx, '" />'
write(nyString,*) '   <Attribute Name="ny" Dimensions="', ny, '" />'
write(nzString,*) '   <Attribute Name="nz" Dimensions="', nz, '" />'
write(dimString0,*) '			Dimensions="', nx, ny, nz, '"'
write(dimString0Vector,*) '			Dimensions="', nx, ny, nz, 3, '"'
write(dimString1,*) '			Dimensions="', nx+1, ny+1, nz+1, '"'

!build loadFile
loadDataFile='<?xml version="1.0" ?>' // NEW_LINE('A') // &
'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>' // NEW_LINE('A') // &
'<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">' // NEW_LINE('A') // &
'  <Domain>' // NEW_LINE('A') // &
'    <!-- Define Dimension as an attribute -->' // NEW_LINE('A') // &
trim(nxString) // NEW_LINE('A') // &
trim(nyString) // NEW_LINE('A') // &
trim(nzString) // NEW_LINE('A') // &
'    <!-- Define Dimension as an attribute -->' // NEW_LINE('A') // &
'    <Grid Name="3DMesh" GridType="Uniform">' // NEW_LINE('A') // &
'      <Topology ' // NEW_LINE('A') // &
'		TopologyType="3DCORECTMESH" ' // NEW_LINE('A') // &
trim(dimString1) // NEW_LINE('A') // &
'	  />' // NEW_LINE('A') // &
'      <Geometry GeometryType="ORIGIN_DXDYDZ">' // NEW_LINE('A') // &
'         <DataItem Name="Origin" NumberType="Float" Dimensions="3" Format="XML">0. 0. 0.</DataItem>' // NEW_LINE('A') // &
'         <DataItem Name="Spacing" NumberType="Float" Dimensions="3" Format="XML">1. 1. 1.</DataItem>' // NEW_LINE('A') // &
'      </Geometry>' // NEW_LINE('A') // &
'      <Attribute '

if (writeTempertureField) then
loadDataFile= trim(loadDataFile) // NEW_LINE('A') // &
'		Name="Temperature" ' // NEW_LINE('A') // &
'		Active="1" ' // NEW_LINE('A') // &
'		AttributeType="Scalar" ' // NEW_LINE('A') // &
'		Center="Cell">' // NEW_LINE('A') // &
'          <DataItem ' // NEW_LINE('A') // &
trim(dimString0) // NEW_LINE('A') // &
'			Endian="Little" ' // NEW_LINE('A') // &
'			NumberType="Float" ' // NEW_LINE('A') // &
'			Precision="4" ' // NEW_LINE('A') // &
'			Format="Binary"' // NEW_LINE('A') // &
'			>resTemp.raw' // NEW_LINE('A') // &
'		  </DataItem>' // NEW_LINE('A') // &
'      </Attribute>' // NEW_LINE('A') // &
'      <Attribute '
end if

if (writeFluxField) then
loadDataFile= trim(loadDataFile) // NEW_LINE('A') // &
'		Name="fluxVector" ' // NEW_LINE('A') // &
'		Active="1" ' // NEW_LINE('A') // &
'		AttributeType="Vector" ' // NEW_LINE('A') // &
'		Center="Cell">' // NEW_LINE('A') // &
'        <!-- Dimensions = nx, ny, nz, components_per_vector-->' // NEW_LINE('A') // &
'          <DataItem            ' // NEW_LINE('A') // &
trim(dimString0Vector) // NEW_LINE('A') // &
'			Endian="Little" ' // NEW_LINE('A') // &
'			NumberType="Float" ' // NEW_LINE('A') // &
'			Precision="4" ' // NEW_LINE('A') // &
'			Format="Binary"' // NEW_LINE('A') // &
'			>resFlux.raw' // NEW_LINE('A') // &
'		  </DataItem>' // NEW_LINE('A') // &
'      </Attribute>' // NEW_LINE('A') // &
'      <Attribute '
end if

loadDataFile= trim(loadDataFile) // NEW_LINE('A') // &
'		Name="greyvalue_phases" ' // NEW_LINE('A') // &
'		Active="1" ' // NEW_LINE('A') // &
'		AttributeType="Scalar" ' // NEW_LINE('A') // &
'		Center="Cell">' // NEW_LINE('A') // &
'          <DataItem ' // NEW_LINE('A') // &
trim(dimString0) // NEW_LINE('A') // &
'			Endian="Little" ' // NEW_LINE('A') // &
'			NumberType="uchar" ' // NEW_LINE('A') // &
'			Format="Binary"' // NEW_LINE('A') // &
'			>' // NEW_LINE('A') // &
trim(sampleName) // NEW_LINE('A') // &
'		  </DataItem>' // NEW_LINE('A') // &
'      </Attribute>' // NEW_LINE('A') // &
'    </Grid>' // NEW_LINE('A') // &
'  </Domain>' // NEW_LINE('A') // &
'</Xdmf>'

