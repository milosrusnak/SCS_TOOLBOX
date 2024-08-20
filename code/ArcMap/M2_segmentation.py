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
       Roux, C., Alber, A., Bertrand, M., Vaudor, L., Piegay, H., submitted. "FluvialCorridor": A new ArcGIS 
       package for multiscale riverscape exploration. Geomorphology
       
@summary: Modul2_segmentation is an open-source python and arcPy code.
          Generate channel segments from polygon. Modul using segmentation centerline from Modul1. 
          
'''

# required libraries and packages 
import os
import sys
import math
import arcpy

#-----------------------------------------------------
# Local variables and input
# input
output_folder = arcpy.GetParameterAsText(0)
inputLayer = arcpy.GetParameterAsText(1)
inputCenterline = arcpy.GetParameterAsText(2)
field_year =arcpy.GetParameterAsText(3)
interval = arcpy.GetParameter(4) 
simplification = arcpy.GetParameter(5) 
deleteTF = arcpy.GetParameter(6)

#local
ws = output_folder.replace(os.sep, '/')
channel_layer = inputLayer.split(";")

arcpy.env.overwriteOutput = True
arcpy.env.workspace = ws

for fc in channel_layer:
    desc = arcpy.Describe(fc)
SR = desc.spatialReference

EA_layer = []
year_list = []
UNI_polygon = []



#===============================================================================
# CODING
#===============================================================================

#----------------------------------------------------
# MAIN PROGRAM
#----------------------------------------------------

#STEP 1 copy layers with the name of the year extracted from atribute table
arcpy.AddMessage("STEP 1 Preprocessing polygons")
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
    arcpy.management.CopyFeatures (fclist,newName)
    EA_layer.append(newName)

#STEP 2 union all channel layer
arcpy.AddMessage("STEP 2 Create union of all polygons")
union_pol = arcpy.analysis.Union(EA_layer, "%ScratchWorkspace%\\union_pol", "ALL")
arcpy.management.AddField(union_pol, "DISS", "SHORT")
with arcpy.da.UpdateCursor(union_pol, "DISS") as cursor:
    for row in cursor:
        row[0] = 1
        cursor.updateRow(row)
union_pol2 = arcpy.management.Dissolve(union_pol, "%ScratchWorkspace%\\union_pol2", "DISS")

#STEP 3 fill holow (create union channel without holow polygon)
arcpy.AddMessage("STEP 3 Converting union of polygons to union without hollows")
with arcpy.da.UpdateCursor(union_pol2, ["SHAPE@"]) as updateCursor:
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
            name_pol= "union_channel.shp"
            arcpy.management.CopyFeatures (new_poly,name_pol)
            arcpy.management.DefineProjection(name_pol, SR)
        else:
            arcpy.management.CopyFeatures (union_pol2,name_pol)
            arcpy.management.DefineProjection(name_pol, SR)
        updateCursor.updateRow(updateRow)

channel = name_pol

#STEP 4 Simplification centerline if is defined in output
if simplification == 0:
    centro = arcpy.management.CopyFeatures (inputCenterline, "%ScratchWorkspace%\\centro")
else:
    arcpy.AddMessage("Simplification of centerline")
    centro = arcpy.management.CopyFeatures (inputCenterline, "centro_simple_{}m.shp".format(simplification))
    arcpy.management.Integrate(centro, simplification)
    

#STEP 5 Split centerline in defined interval 
arcpy.AddMessage("STEP 4 Create longitudinal segments....")
clipPoints = arcpy.management.GeneratePointsAlongLines(centro, "%ScratchWorkspace%\\clipPoints", "DISTANCE", interval, "","END_POINTS")
arcpy.management.AddField(clipPoints, "Distance", "LONG")
arcpy.management.AddField(clipPoints, "seg_rev", "SHORT")

all_rows = [i[0] for i in arcpy.da.SearchCursor(clipPoints,"OID@")]
max_val = max(all_rows)
with arcpy.da.UpdateCursor(clipPoints, ["OID@","Distance", "seg_rev"]) as cursor:
    for row in cursor:
        row[1] = interval
        row[2] = max_val - row[0]
        cursor.updateRow(row)

if interval > 2:
    rad = "{} meters".format(1)
else:
    rad = "{} meters".format(interval/5)
centerlinePointsCLIP = arcpy.management.SplitLineAtPoint(centro, clipPoints, "%ScratchWorkspace%\\centerlinePointsCLIP", rad)

#STEP 6 Combine split line centerline with sequenced points
fm = arcpy.FieldMappings()
fm.addTable(clipPoints)
IDIndex = fm.findFieldMapIndex("seg_rev")
fieldmap = fm.getFieldMap(IDIndex)
field = fieldmap.outputField
field.name = "ID_SEQ"
field.aliasName = "ID_SEQ"
fieldmap.outputField = field
fieldmap.mergeRule = "Max"
fm.replaceFieldMap(IDIndex, fieldmap)

centerSeg = arcpy.analysis.SpatialJoin(centerlinePointsCLIP, clipPoints, "%ScratchWorkspace%\\centerSeg", "JOIN_ONE_TO_ONE", "",fm, "CONTAINS")

#STEP 7 Create midpoints of centerline and Thiessen polygon from this midpoints 
centreMidpoint = arcpy.FeatureVerticesToPoints_management(centerSeg, "%ScratchWorkspace%\\centreMidpoint", "MID")

desc = arcpy.Describe(channel)
extent = desc.extent
    
ext = str(extent).split(" ")
Xmin = float(ext[0].replace(",","."))
Ymin = float(ext[1].replace(",","."))
Xmax = float(ext[2].replace(",","."))
Ymax = float(ext[3].replace(",","."))
    
chaine = str(Xmin) + " " + str(Ymin) + " " + str(Xmax) + " " + str(Ymax)

arcpy.env.extent = chaine

midPointThiessen = arcpy.CreateThiessenPolygons_analysis(centreMidpoint, "midPointThiessen", "ALL")

#STEP 8 create segment 
SegmentClip = arcpy.Clip_analysis(midPointThiessen, channel, "%ScratchWorkspace%\\SegmentClip")

fld = ["Distance", "ID_SEQ"]
fields_to_delete = [field.name for field in arcpy.ListFields(SegmentClip) if not field.required and field.name not in fld]
fields_to_delete.pop() 
for field in fields_to_delete:
    arcpy.DeleteField_management(SegmentClip, field)

name = "Segments_{}m.shp".format(interval)
arcpy.management.CopyFeatures (SegmentClip,name)
arcpy.management.DefineProjection(name, SR)

#===============================================================================
# DELETING processing FILES
#===============================================================================
if deleteTF == True:
    arcpy.AddMessage("Deleting processing files")
    arcpy.Delete_management(channel)
    for i in range(len(EA_layer)):
        arcpy.Delete_management(EA_layer[i])
else:
     arcpy.AddMessage("Processing files preserved in output folder")

#===============================================================================
# DELETING TEMPORARY FILES
#===============================================================================
arcpy.Delete_management(union_pol)
arcpy.Delete_management(union_pol2)
arcpy.Delete_management(clipPoints)
arcpy.Delete_management(centerlinePointsCLIP)
arcpy.Delete_management(centerSeg)
arcpy.Delete_management(centreMidpoint)
arcpy.Delete_management(SegmentClip)
arcpy.Delete_management(midPointThiessen)
if simplification == 0:
    arcpy.Delete_management(centro)