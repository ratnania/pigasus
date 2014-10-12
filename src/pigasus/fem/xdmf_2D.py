# -*- coding: UTF-8 -*-
#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="ratnani"
__date__ ="$Feb 4, 2012 10:12:44 AM$"


class xdmf_2D:
   def __init__(self, prefix, index=0):
      self.index = index
      self.prefix = prefix
      self._f = open("%s%05d.xmf"%(self.prefix,self.index), "w")
      (self._f).write("<?xml version=\"1.0\" ?>\n")
      (self._f).write("<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n")
      (self._f).write("<Xdmf Version=\"2.0\">\n")
      (self._f).write(" <Domain>\n")
      (self._f).write("   <Grid Name=\"Mesh\" GridType=\"Uniform\">\n");
   def grid(self,nx, ny, meshprefix="Phi_3D_d"):
      self._nx = nx
      self._ny = ny
      self._mp = meshprefix
      (self._f).write("     <Topology TopologyType=\"2DRectMesh\" NumberOfElements=\"%d %d\"/>\n"%(self._ny,self._nx));
      (self._f).write("     <Geometry GeometryType=\"VXVY\">\n");
      (self._f).write("       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n"% self._nx)
      (self._f).write("        %s%05d.h5:/rg\n" % (self._mp,self.index))
      (self._f).write("       </DataItem>\n")
      (self._f).write("       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n"% self._ny)
      (self._f).write("        %s%05d.h5:/thetag\n" % (self._mp,self.index))
      (self._f).write("       </DataItem>\n")
      (self._f).write("     </Geometry>\n")
   def field(self,name):
      assert self._nx>0
      assert self._ny>0
      (self._f).write("     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Node\">\n" % name)
      (self._f).write("       <DataItem Dimensions=\"%d %d\" NumberType=\"Float\" "%(self._ny,self._nx))
      (self._f).write("          Precision=\"8\" Format=\"HDF\">\n" )
      (self._f).write("        %s%05d.h5:/%s\n" % (self.prefix,self.index,name))
      (self._f).write("       </DataItem>\n")
      (self._f).write("     </Attribute>\n")
   def close(self):
      (self._f).write("   </Grid>\n")
      (self._f).write(" </Domain>\n")
      (self._f).write("</Xdmf>\n")
      (self._f).close()
