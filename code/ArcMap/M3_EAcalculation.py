# -*- coding: utf-8 -*-

'''
Standalone channel shifting toolbox (SCS Toolbox)
Created on 17 MAY 2024
Last update on 17 MAY 2024
@author: Milos Rusnak

@devoloped at: CNRS - UMR5600 Environnement Ville Societe
               15 Parvis Rene Descartes, BP 7000, 69342 Lyon Cedex 07, France 
               
@contact: geogmilo@savba.sk
          Institute of geography SAS
          Stefanikova 49, 814 73 Bratislava, Slovakia 

                   
@note: Standalone channel shifting toolbox (SCS Toolbox) was developed as extension of the FluvialCorridor toolbox with implemented the centerline 
       extraction approach and segmentation of DGO from FluvialCorridor toolbox.
       For each use of the Channel toolbox leading to a publication, report, presentation or any other
       document, please refer also to the following article :
       Roux, C., Alber, A., Bertrand, M., Vaudor, L., Piegay, H., 2015. "FluvialCorridor": A new ArcGIS 
       package for multiscale riverscape exploration. Geomorphology, 29-37, 242.
       doi: 10.1016/j.geomorph.2014.04.018

@summary: Modul3_EAcalculation is an open-source python and arcPy code.
          Generate channel migration (E - Erosion, A - Accumulation) from channel polygons
          Calculation erosion and deposition of river channel
          
'''

# required libraries and packages 
import arcpy
import os
import sys
import math

#-----------------------------------------------------
# Local variables and input
# input
output_folder = arcpy.GetParameterAsText(0)
inputLayer = arcpy.GetParameterAsText(1)
inputLayer2 = arcpy.GetParameterAsText(2)
field_year = arcpy.GetParameterAsText(3)
centerline_year = arcpy.GetParameterAsText(4)
statistics = arcpy.GetParameterAsText(5)
interval = arcpy.GetParameter(6)
deleteTF = arcpy.GetParameter(7)

#local
ws = output_folder.replace(os.sep, '/')
channel_layer = inputLayer.split(";")
centerline_layer = inputLayer2.split(";")

arcpy.env.overwriteOutput = True
arcpy.env.workspace = ws

for fc in channel_layer:
    desc = arcpy.Describe(fc)
SR = desc.spatialReference

year_list = []
year_check = []
EA_layer = []
CEN_layer = []
UNI_polygon = []
EA_island = []
row_count = []
EAprocess = []
UNIyy = []
EArateList = []

#===============================================================================
# CODING
#===============================================================================

# MAIN PROGRAM DEFINITION
# Orientation detection DEF
def OrientationMask (polygon_old, polygon_young, centerline_old, centerline_young, year_old, year_young):
    """
    This function calculates orientation mask for the channel layer by combination of the channel polygon and centreline. \n
    Vars:\n
    \t polygon_old = channel polygon envelop (without hollows) for older year \n
    \t polygon_young = channel polygon envelop (without hollows) for younger year \n
    \t centerline_old = centerline for older year \n
    \t centerline_young = centerline for younger year \n
    \t year_old = info about year for old polygon \n
    \t year_young = info about year for young polygon \n
    RETURNS: SIDEMASk = channel mask polygon with information about orientation to LEFT and RIGHT side of channel
    """
    # input
    chanOlder = polygon_old
    chanYounger = polygon_young
    cenOlder = centerline_old
    cenYounger = centerline_young
    y1 = year_old
    y2 = year_young

    #A combine old and young layer to mask polygon
    uni = arcpy.analysis.Union([chanOlder, chanYounger], "%ScratchWorkspace%\\uni")
    arcpy.management.AddField(uni, "DISS", "SHORT")
    with arcpy.da.UpdateCursor(uni, "DISS") as cursor:
       for row in cursor:
         row[0] = 1
         cursor.updateRow(row)
    uni2 = arcpy.management.Dissolve(uni, "%ScratchWorkspace%\\uni2", "DISS")
    
    #A remove hollows in the layer
    with arcpy.da.UpdateCursor(uni2, ["SHAPE@"]) as updateCursor:
       for updateRow in updateCursor:
         shape = updateRow[0]
         new_shape = arcpy.Array()
         for part in shape:
            new_part = arcpy.Array()
            #get the first None point index
            first_null_point_index = 0
            for i in range(len(part)):
               if part[i] == None:
                  first_null_point_index = i
                  break
            if first_null_point_index == 0:
               new_shape.add(part)
            else:
               for j in range(first_null_point_index):
                  new_part.add(part[j])
               new_shape.add(new_part)
         if len(new_shape) > 0:
            new_poly = arcpy.Polygon(new_shape)
            name_pol= "UNI_{}_{}.shp".format(y1,y2)
            arcpy.management.CopyFeatures (new_poly,name_pol)
         else:
            arcpy.management.CopyFeatures (uni2,name_pol)
         UNIyy.append(name_pol)
         arcpy.management.DefineProjection(name_pol, SR)
         updateCursor.updateRow(updateRow)

    #B check centerline to touch union boundary
    polygon = arcpy.da.SearchCursor("UNI_{}_{}.shp".format(y1,y2), ('SHAPE@')).next()[0]
    boundary = polygon.boundary()
    
    #B check centerline to borders OLD polygon
    lst_feats = []
    cnt = 0
    with arcpy.da.SearchCursor(cenOlder, ('SHAPE@')) as cursor:
      for row in cursor:
         cnt += 1
         polyline = row[0]
         pnt1 = polyline.firstPoint
         pnt2 = polyline.lastPoint

         pntg1_snap = boundary.snapToLine(pnt1)
         pntg2_snap = boundary.snapToLine(pnt2)
         
         pnt1 = pntg1_snap.firstPoint
         pnt2 = pntg2_snap.firstPoint
         lst_pnts = []
         for part in polyline:
            for pnt in part:
               lst_pnts.append(pnt)

         lst_pnts.insert(0, pnt1)
         lst_pnts.append(pnt2)

         polyline_out = arcpy.Polyline(arcpy.Array(lst_pnts))
         lst_feats.append(polyline_out)
    cenOlder2 = arcpy.CopyFeatures_management(lst_feats, "%ScratchWorkspace%\\cenOlder2")
    
    #B check centerline to borders YOUNG polygon
    lst_feats = []
    cnt = 0
    with arcpy.da.SearchCursor(cenYounger, ('SHAPE@')) as cursor:
      for row in cursor:
         cnt += 1
         polyline = row[0]
         pnt1 = polyline.firstPoint
         pnt2 = polyline.lastPoint

         pntg1_snap = boundary.snapToLine(pnt1)
         pntg2_snap = boundary.snapToLine(pnt2)
         
         pnt1 = pntg1_snap.firstPoint
         pnt2 = pntg2_snap.firstPoint
         lst_pnts = []
         for part in polyline:
            for pnt in part:
               lst_pnts.append(pnt)

         lst_pnts.insert(0, pnt1)
         lst_pnts.append(pnt2)

         polyline_out = arcpy.Polyline(arcpy.Array(lst_pnts))
         lst_feats.append(polyline_out)
    cenYounger2 = arcpy.CopyFeatures_management(lst_feats, "%ScratchWorkspace%\\cenYounger2")

    #C Identification of side mask orientation OLDER
    bufLo = arcpy.analysis.Buffer(cenOlder2, "%ScratchWorkspace%\\bufLo", 1, "LEFT", "ROUND")
    arcpy.management.AddField(bufLo, "SIDEL", "TEXT")
    with arcpy.da.UpdateCursor(bufLo, "SIDEL") as cursor:
      for row in cursor:
         row[0] = "LEFT"
         cursor.updateRow(row)
    
    bufRo = arcpy.analysis.Buffer(cenOlder2, "%ScratchWorkspace%\\bufRo", 1, "RIGHT", "ROUND")
    arcpy.management.AddField(bufRo, "SIDER", "TEXT")
    with arcpy.da.UpdateCursor(bufRo, "SIDER") as cursor:
      for row in cursor:
         row[0] = "RIGHT"
         cursor.updateRow(row)
   
    bufLRo = arcpy.analysis.Union([bufLo, bufRo], "%ScratchWorkspace%\\bufLRo")
    bufLRclipo = arcpy.analysis.Clip(bufLRo, "UNI_{}_{}.shp".format(y1,y2))
    arcpy.management.AddField(bufLRclipo, "SIDE_{}".format(y1), "TEXT")
    with arcpy.da.UpdateCursor(bufLRclipo, ("SIDE_{}".format(y1), "SIDEL", "SIDER")) as cursor:
      for row in cursor:
         row[0] = row[1]+row[2]
         cursor.updateRow(row)
    
    polcnto = arcpy.FeatureToPolygon_management(["UNI_{}_{}.shp".format(y1,y2),cenOlder2], "%ScratchWorkspace%\\polcnto")
    fm = arcpy.FieldMappings()
    fm.addTable(bufLRclipo)
    fm.addTable(polcnto)
    keepers = ["SIDE_{}".format(y1)]
    for field in fm.fields:
      if field.name not in keepers:
         fm.removeFieldMap(fm.findFieldMapIndex(field.name))
    sideMaskOlder = arcpy.analysis.SpatialJoin(polcnto, bufLRclipo, "%ScratchWorkspace%\\sideMaskOlder", "JOIN_ONE_TO_MANY", "",fm, "INTERSECT")
   
    #C Identification of side mask orientation YOUNG
    bufLy = arcpy.analysis.Buffer(cenYounger2, "%ScratchWorkspace%\\bufLy", 1, "LEFT", "ROUND")
    arcpy.management.AddField(bufLy, "SIDEL", "TEXT")
    with arcpy.da.UpdateCursor(bufLy, "SIDEL") as cursor:
      for row in cursor:
         row[0] = "LEFT"
         cursor.updateRow(row)
   
    bufRy = arcpy.analysis.Buffer(cenYounger2, "%ScratchWorkspace%\\bufRy", 1, "RIGHT", "ROUND")
    arcpy.management.AddField(bufRy, "SIDER", "TEXT")
    with arcpy.da.UpdateCursor(bufRy, "SIDER") as cursor:
      for row in cursor:
         row[0] = "RIGHT"
         cursor.updateRow(row)
   
    bufLRy = arcpy.analysis.Union([bufLy, bufRy], "%ScratchWorkspace%\\bufLRy")
    bufLRclipy = arcpy.analysis.Clip(bufLRy, "UNI_{}_{}.shp".format(y1,y2))
    arcpy.management.AddField(bufLRclipy, "SIDE_{}".format(y2), "TEXT")
    with arcpy.da.UpdateCursor(bufLRclipy, ("SIDE_{}".format(y2), "SIDEL", "SIDER")) as cursor:
      for row in cursor:
         row[0] = row[1]+row[2]
         cursor.updateRow(row)
   
    polcnty = arcpy.FeatureToPolygon_management(["UNI_{}_{}.shp".format(y1,y2),cenYounger2], "%ScratchWorkspace%\\polcnt")
    fm = arcpy.FieldMappings()
    fm.addTable(bufLRclipy)
    fm.addTable(polcnty)
    keepers = ["SIDE_{}".format(y2)]
    for field in fm.fields:
      if field.name not in keepers:
         fm.removeFieldMap(fm.findFieldMapIndex(field.name))
    sideMaskYounger = arcpy.analysis.SpatialJoin(polcnty, bufLRclipy, "%ScratchWorkspace%\\sideMaskYounger", "JOIN_ONE_TO_MANY", "",fm, "INTERSECT")

    #D Create SIDE MASK 
    SIDEMASk = arcpy.analysis.Union([sideMaskOlder, sideMaskYounger])
    SIDEMASk2 = arcpy.management.CopyFeatures (SIDEMASk,"SIDEMASk2.shp")
    
    #===============================================================================
    # DELETING TEMPORARY FILES
    #===============================================================================
    arcpy.Delete_management(uni)
    arcpy.Delete_management(uni2)
    arcpy.Delete_management(cenOlder2)
    arcpy.Delete_management(cenYounger2)
    arcpy.Delete_management(bufLo)
    arcpy.Delete_management(bufRo)
    arcpy.Delete_management(bufLRo)
    arcpy.Delete_management(bufLRclipo)
    arcpy.Delete_management(polcnto)
    arcpy.Delete_management(bufLy)
    arcpy.Delete_management(bufRy)
    arcpy.Delete_management(bufLRy)
    arcpy.Delete_management(bufLRclipy)
    arcpy.Delete_management(polcnty)
    arcpy.Delete_management(sideMaskOlder)
    arcpy.Delete_management(sideMaskYounger)
    arcpy.Delete_management(SIDEMASk)

    return SIDEMASk2

#----------------------------------------------------
#----------------------------------------------------
#----------------------------------------------------
# MAIN PROGRAM

#STEP 1 copy layers with the name of the year extracted from atribute table
arcpy.AddMessage("STEP 1 Preprocessing channel polygons")
for fclist in channel_layer:
   fields_search = [f.name for f in arcpy.ListFields(fclist)]
   for field in fields_search:
        if field == field_year:
            field_check = field
   with arcpy.da.SearchCursor(fclist, field_check) as cursor:
      for row in cursor:
        year = row [0]
   newName = "CH_"+ str(year) + ".shp"
   year_list.append(year)
   newNamepath = os.path.join(ws, newName)
   if arcpy.Exists(newNamepath):
        arcpy.AddMessage("{} exists, not copying".format(newNamepath))
   else:
        arcpy.management.CopyFeatures (fclist,newName)
   EA_layer.append(newName)

#STEP 2 read centerline layer year and create centerline list
arcpy.AddMessage("STEP 2 Preprocessing centerlines")
for fclist in centerline_layer:
   fields_search = [f.name for f in arcpy.ListFields(fclist)]
   for field in fields_search:
        if (field.find(centerline_year) !=-1 ):
            field_check = field
   with arcpy.da.SearchCursor(fclist, field_check) as cursor:
      for row in cursor:
        year = row [0]
   newName = "centro_"+ str(year) + ".shp"
   year_check.append(year)
   newNamepath = os.path.join(ws, newName)
   if arcpy.Exists(newNamepath):
        arcpy.AddMessage("{} exists, not copying".format(newNamepath))
   else:
        arcpy.management.CopyFeatures (fclist,newName)
   CEN_layer.append(newName)

#STEP 3 sort layers from younger to older and check
EA_layer_sort = sorted(EA_layer)
year_sort = sorted(year_list) 
CEN_layer_sort = sorted(CEN_layer)
year_check_sort = sorted(year_check) 

if year_sort == year_check_sort:
   arcpy.AddMessage("Checked out polygons and centerlines")
else:
   arcpy.AddMessage("!!!!! Chanel polygons years do not match centerline years !!!!")

#STEP 4 simplify channel atribute table 
for fc in EA_layer_sort:
    fields_to_delete = [field.name for field in arcpy.ListFields(fc) if not field.required]
    fields_to_delete.pop() 
    for field in fields_to_delete:
      arcpy.DeleteField_management(fc, field)

#STEP 5 create new field with year fieldname and year value
for i in range(len(EA_layer_sort)):
   arcpy.management.AddField(EA_layer_sort[i], "y_{}".format(year_sort[i]), "LONG")
   with arcpy.da.UpdateCursor(EA_layer_sort[i], "y_{}".format(year_sort[i])) as cursor:
      for row in cursor:
         row[0] = year_sort[i]
         cursor.updateRow(row)

#STEP 6 fill holow (create channel without holow polygon)
arcpy.AddMessage("STEP 3 Converting polygons to polygon without hollows")
for n in range(len(EA_layer_sort)):
   with arcpy.da.UpdateCursor(EA_layer_sort[n], ["SHAPE@"]) as updateCursor:
      for updateRow in updateCursor:
         shape = updateRow[0]
         new_shape = arcpy.Array()
         for part in shape:
            new_part = arcpy.Array()
            #get the first None point index
            first_null_point_index = 0
            for i in range(len(part)):
               if part[i] == None:
                  first_null_point_index = i
                  break
            if first_null_point_index == 0:
               new_shape.add(part)
            else:
               for j in range(first_null_point_index):
                  new_part.add(part[j])
               new_shape.add(new_part)
         if len(new_shape) > 0:
            new_poly = arcpy.Polygon(new_shape)
            name_pol= "POL_{}.shp".format(year_sort[n])
            arcpy.management.CopyFeatures (new_poly,name_pol)
            arcpy.management.DefineProjection(name_pol, SR)
         else:
            arcpy.management.CopyFeatures (EA_layer_sort[n],name_pol)
            arcpy.management.DefineProjection(name_pol, SR)
         UNI_polygon.append(name_pol)
         updateCursor.updateRow(updateRow)

#STEP 7 import hollows as islands
arcpy.AddMessage("STEP 4 Create polygons with atribute channel and island")
for i in range(len(UNI_polygon)):
   inter_out = "EA_island_{}.shp".format(year_sort[i])
   arcpy.analysis.Union([UNI_polygon[i], EA_layer_sort[i]], inter_out)
   arcpy.management.AddField(inter_out, "TYP_{}".format(year_sort[i]), "TEXT")
   with arcpy.da.UpdateCursor(inter_out, ("y_{}".format(year_sort[i]), "TYP_{}".format(year_sort[i]))) as cursor:
      for row in cursor:
         if row[0] > 0:
            row[1] = "channel"
         else:
            row[1] = "island" 
         row[0] = year_sort[i]
         cursor.updateRow(row)
   EA_island.append(inter_out) 
   i = i+1

#STEP 8 check centerline processing
for fclist in CEN_layer_sort:
    row_num = 0
    with arcpy.da.UpdateCursor(fclist, ["SHAPE@"]) as updateCursor:
        for updateRow in updateCursor:
            row_count = row_num +1
            if row_count == 1:
               arcpy.AddMessage("layer {} is OK".format(fclist))
            else: 
               arcpy.AddMessage("layer {} has issue with centerline. Check centerline topology.".format(fclist))
            updateCursor.updateRow(updateRow)    

#STEP 9 calculate in-channel process and side orientation labeling 
for i in range(len(EA_island)-1):
   input1 = UNI_polygon[i]
   input2 = UNI_polygon[i+1]
   input3 = CEN_layer_sort[i]
   input4 = CEN_layer_sort[i+1]
   input5 = year_sort[i]
   input6 = year_sort[i+1]

   arcpy.AddMessage("STEP 5 Calculate in-channel proces of erosion and deposition for years {} and {}".format(input5, input6))

   sideMask = OrientationMask (input1, input2, input3, input4, input5, input6)
   
   inputEA1 = EA_island[i]
   inputEA2  = EA_island[i+1]
   y1 = input5
   y2 = input6
   name = "EA_processes{}_{}.shp".format(y1,y2)

   fld = ["y_{}".format(y1), "TYP_{}".format(y1), "y_{}".format(y2), "TYP_{}".format(y2)]
   
   unionEA = arcpy.analysis.Union([inputEA1,inputEA2], "%ScratchWorkspace%\\unionEA")
   arcpy.management.AddField(unionEA, "EA", "TEXT")
   with arcpy.da.UpdateCursor(unionEA, ["EA"] + fld) as cursor:
      for row in cursor:
         if row[1] != y1 and row[4] == "channel":
            row[0] = "erosion"
         elif row[3] != y2 and row[2] == "channel":
            row[0] = "deposition"
         elif row[1] != y1 and row[3] != y2:
            row[0] = "hollow"
         elif row[2] == "island" and row[4] == "channel":
            row[0] = "island_erosion"
         elif row[2] == "channel" and row[4] == "island":
            row[0] = "island_deposition"
         elif (row[2] == row[4]) or (row[1] != y1 and row[4] == "island") or (row[3] != y2 and row[2] == "island"):
            row[0] = "stable"
         cursor.updateRow(row)  

   unionEAmask = arcpy.analysis.Union([unionEA,sideMask], "%ScratchWorkspace%\\unionEAmask")
   
   fld1 = ["EA", "SIDE_{}".format(y1), "SIDE_{}".format(y2),"y_{}".format(y1), "y_{}".format(y2)]
   
   arcpy.management.AddField(unionEAmask, "direction", "TEXT")
   with arcpy.da.UpdateCursor(unionEAmask, ["direction"] + fld1) as cursor:
      for row in cursor:
         if row[1] == "deposition":
            row[0] = row[3]
         elif row[1] == "erosion":
            row[0] = row[2]
         elif row[4] != y1 and row[5] != y2:
            row[1] = "hollow"
            row[0] = row[2]
         else:
            row[0] = "in-channel process"
         cursor.updateRow(row)  
   
   fld2 = ["EA", "direction"]
   
   arcpy.management.AddField(unionEAmask, "migration", "TEXT")
   with arcpy.da.UpdateCursor(unionEAmask, ["migration"] + fld2) as cursor:
      for row in cursor:
         if row[2] == "in-channel process":
            row[0] = "in-channel process"
         elif row[1] == "erosion" or row[1] == "hollow":
            row[0] = "erosion_{}".format(row[2])
         elif row[1] == "deposition":
            row[0] = "deposition_{}".format(row[2])
         cursor.updateRow(row)  

   arcpy.management.AddField(unionEAmask, "period", "TEXT")
   with arcpy.da.UpdateCursor(unionEAmask, "period") as cursor:  
      for row in cursor:
         row[0] = "{}_{}".format(y1,y2) 
         cursor.updateRow(row)
   
   arcpy.management.AddField(unionEAmask, "span_year", "SHORT")
   y3 = str(y1)
   y4 = str(y2) 
   y5 = int(y3[:4])
   y6 = int(y4[:4])
   #year_older = int(str(y2[:3]))
   with arcpy.da.UpdateCursor(unionEAmask, "span_year") as cursor:  
      for row in cursor:
         if y6 - y5 > 0:
            row[0] = y6 - y5
         else:
            row[0] = 1
         cursor.updateRow(row)

   fld3 = ["y_{}".format(y1), "TYP_{}".format(y1), "y_{}".format(y2), "TYP_{}".format(y2), "EA", "direction", "period", "span_year", "migration"]
   fields_to_delete = [field.name for field in arcpy.ListFields(unionEAmask) if not field.required and field.name not in fld3]
   fields_to_delete.pop() 
   for field in fields_to_delete:
      arcpy.DeleteField_management(unionEAmask, field)

   #STEP 10 final data export
   unionEAdiss = arcpy.management.Dissolve(unionEAmask, "%ScratchWorkspace%\\unionEAdiss", ["EA", "direction", "period", "span_year", "migration"])
   arcpy.management.CopyFeatures (unionEAdiss,name)
   arcpy.management.DefineProjection(name, SR)
   EAprocess.append(name)
   i = i+1

if len(statistics) != 0:
    #Calculate EA statistics and erosion intensity for segment
    arcpy.AddMessage("STEP 6 Calculate EA_process layer with channel segments information: erosion intensity and migration rate for every segment")
    for i in range(len(EAprocess)):
      unionEAseg = arcpy.Intersect_analysis ([EAprocess[i], statistics], "%ScratchWorkspace%\\unionEAseg", "ALL")
      unionEAsegSingle = arcpy.management.MultipartToSinglepart(unionEAseg, "%ScratchWorkspace%\\unionEAsegSingle")
      arcpy.DeleteField_management(unionEAsegSingle, "ORIG_FID")
      name2 = "EAsegments_{}_{}.shp".format(year_sort[i],year_sort[i+1])
      arcpy.management.CopyFeatures (unionEAsegSingle,name2)
      arcpy.management.DefineProjection(name2, SR)

            
    #Calculate erosion intensity for segment and add info to EA statistics
    arcpy.AddMessage("STEP 7 Calculate erosion intensity and migration rate for every segment")
    for i in range(len(EAprocess)):
      EArate = "EA_rate_{}_{}.shp".format(year_sort[i], year_sort[i+1]) 
      diss = arcpy.management.Dissolve(EAprocess[i], "%ScratchWorkspace%\\diss", ["span_year", "migration", "period"])
      arcpy.Intersect_analysis ([diss, statistics], EArate, "ALL")
      arcpy.management.AddField(EArate, "EA_rate_A", "DOUBLE")
      arcpy.management.AddField(EArate, "EA_rate_m", "DOUBLE")
      fld4 = ["migration", "span_year", "EA_rate_A", "EA_rate_m"]
      with arcpy.da.UpdateCursor(EArate, ["SHAPE@AREA"] + fld4) as cursor:  
         for row in cursor:
            if row[1] == "in-channel process":
               row[3] = 0
               row[4] = 0
            elif row[1] == "erosion_LEFT" or row[1] == "erosion_RIGHT":
               row[3] = row[0]/row[2]
               row[4] = (row[0]/row[2])/interval
            elif row[1] == "deposition_LEFT" or row[1] == "deposition_RIGHT":
               row[3] = (row[0]/row[2]) * (-1) 
               row[4] = ((row[0]/row[2]) * (-1))/interval
            cursor.updateRow(row)
      
      EArateList.append(EArate) 
      i = i+1 

    #===============================================================================
    # DELETING TEMPORARY FILES
    #===============================================================================
    arcpy.Delete_management(unionEAseg)
    arcpy.Delete_management(unionEAsegSingle)
    arcpy.Delete_management(diss)
      

#===============================================================================
# DELETING processing FILES
#===============================================================================
if deleteTF == True:
   arcpy.AddMessage("Deleting processing files")
   for i in range(len(EA_layer)):
      arcpy.Delete_management(EA_layer[i])
   for i in range(len(UNI_polygon)):
      arcpy.Delete_management(UNI_polygon[i])
   for i in range(len(EA_island)):
      arcpy.Delete_management(EA_island[i])
   for i in range(len(UNIyy)):
      arcpy.Delete_management(UNIyy[i])
else:
   arcpy.AddMessage("Processing files preserved in output folder")

#===============================================================================
# DELETING TEMPORARY FILES
#===============================================================================
arcpy.Delete_management("SIDEMASk2.shp")
arcpy.Delete_management(unionEA)
arcpy.Delete_management(unionEAmask)
arcpy.Delete_management(unionEAdiss)