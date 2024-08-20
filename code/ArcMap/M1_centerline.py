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
       
@summary: Modul1_Centerline is an open-source python and arcPy code.
          Generate channel centerline from polygon
          
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
field_year = arcpy.GetParameterAsText(2)
selection = arcpy.GetParameter(3) 
deleteTF = arcpy.GetParameter(4)

#local
ws = output_folder.replace(os.sep, '/')
channel_layer = inputLayer.split(";")

arcpy.env.overwriteOutput = True
arcpy.env.workspace = ws
arcpy.env.extent = "MAXOF"

for fc in channel_layer:
    desc = arcpy.Describe(fc)
SR = desc.spatialReference

year_list = []
EA_layer = []
UNI_polygon = []
EA_hol = []
centro_list = []

#===============================================================================
# CODING
#===============================================================================

# MAIN PROGRAM DEFINITION
# Centerline detection DEF
def Centro (channel):
     """
     This function calculates centerline for the polygon evelope. \n
     Vars:\n
     \t channel = polygon \n
     RETURNS: centerline = line feature
     """
     #import and pre-process channel data (densify polygons with regular distribution of vertices)
     poly = arcpy.management.CopyFeatures (channel, "%ScratchWorkspace%\\poly")
     polyPOINT = arcpy.management.FeatureVerticesToPoints(poly, "%ScratchWorkspace%\\polyPOINT")
     vert_count = float(arcpy.GetCount_management(polyPOINT).getOutput(0))
     poly_lenght = float(sum(row[0] for row in arcpy.da.SearchCursor(poly, 'SHAPE@LENGTH')))
     density = (poly_lenght / vert_count)
     arcpy.edit.Densify(poly, "DISTANCE", density)
     polyToLine = arcpy.management.PolygonToLine(poly, "%ScratchWorkspace%\\polyToLine")

     #create Thessen polygons  
     polyPOINT2 = arcpy.management.FeatureVerticesToPoints(poly, "%ScratchWorkspace%\\polyPOINT2")
     ThiessenPOLY = arcpy.analysis.CreateThiessenPolygons(polyPOINT2, "%ScratchWorkspace%\\ThiessenPOLY")
     Thiessen = arcpy.analysis.Clip(ThiessenPOLY, poly, "%ScratchWorkspace%\\Thiessen")

     #selection of the Thiessen line near the centerline of the polygon (with errors and small line)
     ThiessenToline = arcpy.management.PolygonToLine(Thiessen, "%ScratchWorkspace%\\ThiessenToline")
     thiesTL = arcpy.MakeFeatureLayer_management(ThiessenToline)
     polyTL = arcpy.MakeFeatureLayer_management(polyToLine)
     thiessenSelect = arcpy.management.SelectLayerByLocation(thiesTL, "INTERSECT", polyTL, "", "NEW_SELECTION", "INVERT")
     rawCenter = arcpy.MakeFeatureLayer_management(thiessenSelect) 

     #Thiessen cenerline raw cleaning by removing all lines perpedicular to channel banks
     arcpy.management.AddField(rawCenter, "ANGLE", "DOUBLE")
     rows = arcpy.UpdateCursor(rawCenter)
     shapeN = arcpy.Describe(rawCenter).shapeFieldName
     for row in rows:
       fc = row.getValue(shapeN)
       dx = fc.lastPoint.X - fc.firstPoint.X
       dy = fc.lastPoint.Y - fc.firstPoint.Y
       radian = math.atan2(dy,dx)
       degrees = radian * 180 / math.pi
       if degrees >= 0:
         Ang = degrees
       else:
         Ang = 180 - abs(degrees)
       row.ANGLE = Ang
       rows.updateRow(row)     
    
     arcpy.analysis.Near(rawCenter, polyToLine, "", "", "ANGLE",  "PLANAR")   
     arcpy.management.AddField(rawCenter, "NEAR_edit", "DOUBLE")   
     arcpy.management.AddField(rawCenter, "DIFF", "DOUBLE")
     with arcpy.da.UpdateCursor(rawCenter, ("ANGLE","NEAR_ANGLE", "NEAR_edit", "DIFF")) as cursor:
       for row in cursor:
         if row[1] >= 0:
             row[2] = row[1]
         else:
             row[2] = 180 - abs(row[1])
         row[3] = abs(row[0] - row [2])
         cursor.updateRow(row)
     selectCenter = arcpy.management.SelectLayerByAttribute(rawCenter, "NEW_SELECTION",  '"DIFF" > 50 and "DIFF" < 130') 
     cleanCenter = arcpy.MakeFeatureLayer_management(selectCenter)   
    
     #Clean centerline from small unconnected line 
     arcpy.management.AddField(cleanCenter, "DISS", "SHORT")
     with arcpy.da.UpdateCursor(cleanCenter, "DISS") as cursor:
         for row in cursor:
             row[0] = 1
             cursor.updateRow(row)
     centroDiss = arcpy.management.Dissolve(cleanCenter, "%ScratchWorkspace%\\centroDiss", "DISS", "", "SINGLE_PART", "UNSPLIT_LINES")
     with arcpy.da.UpdateCursor(centroDiss, ["SHAPE@"]) as updateCursor:
         for row in updateCursor:
             shape = row[0]
             line = shape.getPart(0)
             ptscount = line.count
             if ptscount < 4:
                 updateCursor.deleteRow()
     tolerance = 1*density
    
     arcpy.management.Integrate(centroDiss, tolerance)
     centroDiss2 = arcpy.management.Dissolve(centroDiss, "%ScratchWorkspace%\\centroDiss2", "DISS", "", "SINGLE_PART", "DISSOLVE_LINES")
    
     #extent line to the borders
     arcpy.management.AddField(centroDiss2, "Centerln", "TEXT")
     with arcpy.da.UpdateCursor(centroDiss2, "Centerln") as cursor:
         for row in cursor:
             row[0] = "centerline"
             cursor.updateRow(row)
     polyCentro = arcpy.Merge_management ([centroDiss2, polyTL], "%ScratchWorkspace%\\polycentro")
     polyCentro2 = arcpy.management.CopyFeatures(polyCentro, "polyCentro2.shp") 
     arcpy.edit.ExtendLine(polyCentro2, "", "FEATURE")
     polyCentroATR = arcpy.MakeFeatureLayer_management(polyCentro2) 
     selectpolyCenter = arcpy.management.SelectLayerByAttribute( polyCentroATR, "NEW_SELECTION", "Centerln = 'centerline'") 
     centerline = arcpy.MakeFeatureLayer_management(selectpolyCenter)
     centerline2 = arcpy.management.CopyFeatures (centerline,"centerline2.shp")

     #===============================================================================
     # DELETING TEMPORARY FILES
     #===============================================================================
     arcpy.Delete_management(poly)
     arcpy.Delete_management(polyPOINT)
     arcpy.Delete_management(polyToLine)
     arcpy.Delete_management(polyPOINT2)
     arcpy.Delete_management(ThiessenPOLY)
     arcpy.Delete_management(Thiessen)
     arcpy.Delete_management(ThiessenToline)
     arcpy.Delete_management(centroDiss)
     arcpy.Delete_management(centroDiss2)
     arcpy.Delete_management(polyCentro)
     arcpy.Delete_management(polyCentro2)
     arcpy.Delete_management(centerline)
        
     return centerline2
           
#----------------------------------------------------
#----------------------------------------------------
#----------------------------------------------------
# MAIN PROGRAM
#----------------------------------------------------

###################################################
#### INDIVIDUAL CENTERLINE (selection = FALSE) ####
###################################################

if selection == False:
    #STEP 1 copy layers with the name of the year extracted from atribute table
    arcpy.AddMessage("Calculation individual centerlines")
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

    #STEP 2 sort layers from younger to older 
    EA_layer_sort = sorted(EA_layer)
    year_sort = sorted(year_list) 

    #STEP 3 fill holow (create channel without holow polygon)
    arcpy.AddMessage("STEP 2 Converting input polygons to polygons without hollows")
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

    #STEP 4 create centerline
    arcpy.AddMessage("STEP 3 Create centerline for every single channel....")
    for i in range(len(UNI_polygon)):
        inter_out = Centro(UNI_polygon[i])
        name_out = "centro_{}.shp".format(year_sort[i])
        arcpy.management.CopyFeatures (inter_out,name_out)
        arcpy.management.DefineProjection(name_out, SR)
        centro_list.append(name_out) 
        
        fields_to_delete = [field.name for field in arcpy.ListFields(name_out) if not field.required]
        fields_to_delete.pop() 
        for field in fields_to_delete:
            arcpy.DeleteField_management(name_out, field)

        arcpy.management.AddField(name_out, "cnt", "LONG")
        with arcpy.da.UpdateCursor(name_out, "cnt") as cursor:
            for row in cursor:
                row[0] = year_sort[i]
                cursor.updateRow(row)
        i = i+1

        arcpy.management.AddGeometryAttributes(name_out, "LENGTH")
        
    #===============================================================================
    # DELETING processing FILES
    #===============================================================================
    if deleteTF == True:
        arcpy.AddMessage("Deleting processing files")
        for i in range(len(EA_layer_sort)):
            arcpy.Delete_management(EA_layer_sort[i])
        for i in range(len(UNI_polygon)):
            arcpy.Delete_management(UNI_polygon[i])

    else:
        arcpy.AddMessage("Processing files preserved in output folder")

    #===============================================================================
    # DELETING TEMPORARY FILES
    #===============================================================================
    arcpy.Delete_management("centerline2.shp")


  
####################################################
#### SEGMENTATION CENTERLINE (selection = TRUE) ####
####################################################
if selection == True:
    #STEP 1 copy layers with the name of the year extracted from atribute table
    arcpy.AddMessage("Calculation segmentation centerline")
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
        newNamepath = output_folder + "\CH_"+ str(year) + ".shp" 
        checkName = fclist.replace(os.sep, '\\')
        if newNamepath != checkName:
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
    
    #STEP 3 fill holow in union polygon (create union channel without holow polygon)
    arcpy.AddMessage("STEP 3 Converting input polygons")
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
    
    #STEP 4 create layer centerline for union polygon
    arcpy.AddMessage("STEP 4 Create centerline for union of all polygons....")
    channel = "union_channel.shp"
    center_out = Centro(channel)
    name_out ="SegCenterline.shp"
    arcpy.management.CopyFeatures (center_out,name_out)
    arcpy.management.DefineProjection(name_out, SR)

    fields_to_delete = [field.name for field in arcpy.ListFields(name_out) if not field.required]
    fields_to_delete.pop() 
    for field in fields_to_delete:
        arcpy.DeleteField_management(name_out, field)

    arcpy.management.AddGeometryAttributes(name_out, "LENGTH")
    
    #===============================================================================
    # DELETING processing FILES
    #===============================================================================
    if deleteTF == True:
        arcpy.AddMessage("Deleting processing files")
        arcpy.Delete_management("union_channel.shp")
        for i in range(len(EA_layer)):
            arcpy.Delete_management(EA_layer[i])

    else:
        arcpy.AddMessage("Processing files saved in output folder")

    #===============================================================================
    # DELETING TEMPORARY FILES
    #===============================================================================
    arcpy.Delete_management(union_pol)
    arcpy.Delete_management(union_pol2)
    arcpy.Delete_management("centerline2.shp")