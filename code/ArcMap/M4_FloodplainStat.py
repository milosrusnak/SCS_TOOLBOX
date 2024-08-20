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
       
@summary: Modul4_FloodplainStat is an open-source python and arcPy code.
          Generate floodplain age map (FAM), height above channel (HACH) and vegetation canopy height model (CHM)
          Adding segments from Modul2 create statistic file about longitudinal position, floodplain age dating, height above channel and 
          vegetation height statistics of floodplain
          
'''

# required libraries and packages 
import os
import sys
import math
import arcpy
from arcpy import env
from arcpy.sa import *

#-----------------------------------------------------
# Local variables and input
# input
output_folder = arcpy.GetParameterAsText(0)
inputLayer = arcpy.GetParameterAsText(1)
field_year = arcpy.GetParameterAsText(2)
dem = arcpy.GetParameterAsText(3)
dsm = arcpy.GetParameterAsText(4)
flow = arcpy.GetParameterAsText(5)
segments = arcpy.GetParameterAsText(6)
deleteTF = arcpy.GetParameter(7)

#local
ws = output_folder.replace(os.sep, '/') 
channel_layer = inputLayer.split(";")

arcpy.env.overwriteOutput = True
arcpy.env.workspace = ws

for fc in channel_layer:
    desc = arcpy.Describe(fc)
SR = desc.spatialReference

fam = []
year_list = []
U_layer = []
fields_fam  = []



#===============================================================================
# CODING
#===============================================================================

#----------------------------------------------------
# MAIN PROGRAM
#----------------------------------------------------
# PART 0 selecting optional layer input

if len(dem) == 0 or len(flow) == 0:
    arcpy.AddMessage("Not calculate height above channel (HACH) and canopy height model (CHM)")
    k = 3
elif len(dem) != 0 and len(flow) != 0 and len(dsm) == 0:
    arcpy.AddMessage("Not calculate canopy height model (CHM)")
    k = 2
    cell = arcpy.GetRasterProperties_management(dem, "CELLSIZEX")
    cellsz = cell.getOutput(0)
    cellsize = int(cellsz)
    outTrend = "%ScratchWorkspace%\\outTrend.tif"
    detrended = "%ScratchWorkspace%\\detrended.tif"
    hachTab = "%ScratchWorkspace%\\hachTab"
elif len(dem) != 0 and len(flow) != 0 and len(dsm) != 0:
    arcpy.AddMessage("Calculate all statistics included HACH and CHM")
    k = 1
    cell = arcpy.GetRasterProperties_management(dem, "CELLSIZEX")
    cellsz = cell.getOutput(0)
    cellsize = int(cellsz)
    outTrend = "%ScratchWorkspace%\\outTrend.tif"
    detrended = "%ScratchWorkspace%\\detrended.tif"
    hachTab = "%ScratchWorkspace%\\hachTab"
    chm = "%ScratchWorkspace%\\chm.tif"
    chm_clear = "%ScratchWorkspace%\\chm_clear.tif"
    vegTab = "%ScratchWorkspace%\\vegTab"

    #####################################
    #### ALL data statistics (k = 1) ####
    #####################################
if k == 1 :
    #PART A calculation FAM
    #union all channel layer
    arcpy.AddMessage("STEP 1 Union all channel polygons")
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
        U_layer.append(newName)
    
    for i in range(len(U_layer)):
        fldnames = [field.name for field in arcpy.ListFields(U_layer[i])]
        fi = "y{}".format(year_list[i])
        if fi not in fldnames:
            arcpy.management.AddField(U_layer[i], "y{}".format(year_list[i]), "LONG")
            with arcpy.da.UpdateCursor(U_layer[i], "y{}".format(year_list[i])) as cursor:
                for row in cursor:
                    row[0] = year_list[i]
                    cursor.updateRow(row)
        fields_fam.append("y{}".format(year_list[i]))

    union = arcpy.Union_analysis (U_layer, "%ScratchWorkspace%\\union", "ALL")

    #cleaning fields
    fields_to_delete = [field.name for field in arcpy.ListFields(union) if not field.required and field.name not in fields_fam]
    fields_to_delete.pop() 
    for field in fields_to_delete:
        arcpy.DeleteField_management(union, field)

    #calculate floodplain age (FAM)
    arcpy.AddMessage("STEP 2 Create Floodplain Age Map layer (FAM)")
    arcpy.management.AddField(union, "FAM", "LONG")
    fields_fam.append("FAM")

    with arcpy.da.UpdateCursor(union, fields_fam) as cursor:
        for row in cursor:
            maxF=row[0]
            for i in range(len(row)-1):
                if row[i] > maxF:
                    maxF = row[i]
            row[len(row)-1]=maxF  
            cursor.updateRow(row)    
    name = "fam_layer.shp"
    union2 = arcpy.management.MultipartToSinglepart(union, "%ScratchWorkspace%\\union2")
    arcpy.DeleteField_management(union2, "ORIG_FID")
    arcpy.management.CopyFeatures (union2,name)

    fld = ["FAM"]
    fields_to_delete = [field.name for field in arcpy.ListFields(union2) if not field.required and field.name not in fld]
    fields_to_delete.pop() 
    for field in fields_to_delete:
        arcpy.DeleteField_management(union2, field)

    #PART B calculation HACH
    arcpy.AddMessage("STEP 3 Create Height Above Channel layer (HACH)")
    arcpy.env.extent = dem
    arcpy.env.snapRaster = dem
    arcpy.env.mask = dem

    Fpath_dis = arcpy.management.CopyFeatures(flow, "%ScratchWorkspace%\\Fpath_dis")
    distance= cellsize*5
    arcpy.Densify_edit(Fpath_dis, "DISTANCE", distance)
    Fpath_point = arcpy.FeatureVerticesToPoints_management(Fpath_dis, "%ScratchWorkspace%\\Fpath_point", "All")
    Fpath_point_Z = arcpy.sa.ExtractValuesToPoints(Fpath_point, dem, "%ScratchWorkspace%\\Fpath_point_Z", "NONE", "VALUE_ONLY")
    Fpath_point_Z2 = arcpy.MakeFeatureLayer_management(Fpath_point_Z,"Fpath_point_Z.shp")
    inpt =  "{} RASTERVALU PointElevation".format(Fpath_point_Z2)
    outTrend = arcpy.ddd.TopoToRaster(inpt, "outTrend.tif", cellsize)   
    detrended = Minus (dem, outTrend)
    detrended.save(output_folder + "/" + "DED.tif")

    #PART C calcualte CHM
    arcpy.AddMessage("STEP 4 Create Canopy Height Model layer (CHM)")
    chm = Minus (dsm, dem)
    whereClause = "VALUE <= 0"
    chm_clear = arcpy.sa.SetNull(chm, chm, whereClause)
    chm_clear.save(output_folder + "/" + "veget_CHM.tif")

    #PART D calculate floodplain zone data properties
    arcpy.AddMessage("STEP 5 Create floodplain zone statistic with channel segments")
    unionFAMseg = arcpy.Intersect_analysis ([union2, segments], "%ScratchWorkspace%\\unionFAMseg", "ALL")

    unionFAMsegSingle = arcpy.management.MultipartToSinglepart(unionFAMseg, "%ScratchWorkspace%\\unionFAMsegSingle")
    arcpy.DeleteField_management(unionFAMsegSingle, "ORIG_FID")

    # HACH creation
    hachTab = arcpy.sa.ZonalStatisticsAsTable(unionFAMsegSingle, "FID", detrended, "hachTab","DATA", "ALL")

    fldlst = ["MIN", "MAX", "RANGE", "MEAN", "STD", "SUM"]
    fieldList = [field.name for field in arcpy.ListFields(hachTab) if field.name in fldlst]
   
    fldlstE = []
    for field in fieldList:
        arcpy.AddField_management(hachTab, "e_" + field, "DOUBLE")
        with arcpy.da.UpdateCursor(hachTab, (field,"e_{}".format(field))) as cursor:
            for row in cursor:
                row[1] = row[0]
                cursor.updateRow(row) 
        arcpy.DeleteField_management(hachTab, field)
        fldlstE.append("e_" + field)

    # VEG creation   
    vegTab = arcpy.sa.ZonalStatisticsAsTable(unionFAMsegSingle, "FID", chm_clear, "vegTab","DATA", "ALL")

    fldlst = ["MIN", "MAX", "RANGE", "MEAN", "STD", "SUM"]
    fieldList = [field.name for field in arcpy.ListFields(vegTab) if field.name in fldlst]
   
    fldlstV = []
    for field in fieldList:
        arcpy.AddField_management(vegTab, "v_" + field, "DOUBLE")
        with arcpy.da.UpdateCursor(vegTab, (field,"v_{}".format(field))) as cursor:
            for row in cursor:
                row[1] = row[0]
                cursor.updateRow(row) 
        arcpy.DeleteField_management(vegTab, field)
        fldlstV.append("v_" + field)

    #DATA UNION
    arcpy.management.JoinField(unionFAMsegSingle, "OBJECTID", hachTab, "OBJECTID",fldlstE)
    arcpy.management.JoinField(unionFAMsegSingle, "OBJECTID", vegTab, "OBJECTID",fldlstV)
    name2 = "M4stattistics_all.shp"
    arcpy.management.CopyFeatures (unionFAMsegSingle,name2)
    arcpy.management.DefineProjection(name2, SR)

    ##############################################
    #### FAM and HACH data statistics (k = 2) ####
    ##############################################
if k == 2 :
    #PART A calculation FAM
    #union all channel layer
    arcpy.AddMessage("STEP 1 Union all channel polygons")
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
        U_layer.append(newName)
    
    for i in range(len(U_layer)):
        fldnames = [field.name for field in arcpy.ListFields(U_layer[i])]
        fi = "y{}".format(year_list[i])
        if fi not in fldnames:
            arcpy.management.AddField(U_layer[i], "y{}".format(year_list[i]), "LONG")
            with arcpy.da.UpdateCursor(U_layer[i], "y{}".format(year_list[i])) as cursor:
                for row in cursor:
                    row[0] = year_list[i]
                    cursor.updateRow(row)
        fields_fam.append("y{}".format(year_list[i]))

    union = arcpy.Union_analysis (U_layer, "%ScratchWorkspace%\\union", "ALL")

    #cleaning fields
    
    fields_to_delete = [field.name for field in arcpy.ListFields(union) if not field.required and field.name not in fields_fam]
    fields_to_delete.pop() 
    for field in fields_to_delete:
        arcpy.DeleteField_management(union, field)

    #calculate floodplain age (FAM)
    arcpy.AddMessage("STEP 2 Create Floodplain Age Map layer (FAM)")
    arcpy.management.AddField(union, "FAM", "SHORT")
    fields_fam.append("FAM")

    with arcpy.da.UpdateCursor(union, fields_fam) as cursor:
        for row in cursor:
            maxF=row[0]
            for i in range(len(row)-1):
                if row[i] > maxF:
                    maxF = row[i]
            row[len(row)-1]=maxF  
            cursor.updateRow(row)    
    name = "fam_layer.shp"
    union2 = arcpy.management.MultipartToSinglepart(union, "%ScratchWorkspace%\\union2")
    arcpy.DeleteField_management(union2, "ORIG_FID")
    arcpy.management.CopyFeatures (union2,name)

    fld = ["FAM"]
    fields_to_delete = [field.name for field in arcpy.ListFields(union2) if not field.required and field.name not in fld]
    fields_to_delete.pop() 
    for field in fields_to_delete:
        arcpy.DeleteField_management(union2, field)

    #PART B calculation HACH
    arcpy.AddMessage("STEP 3 Create Height Above Channel layer (HACH)")
    arcpy.env.extent = dem
    arcpy.env.snapRaster = dem
    arcpy.env.mask = dem
    Fpath_dis = arcpy.management.CopyFeatures(flow, "%ScratchWorkspace%\\Fpath_dis")
    distance= cellsize*5
    arcpy.Densify_edit(Fpath_dis, "DISTANCE", distance)
    Fpath_point = arcpy.FeatureVerticesToPoints_management(Fpath_dis, "%ScratchWorkspace%\\Fpath_point", "All")
    Fpath_point_Z = arcpy.sa.ExtractValuesToPoints(Fpath_point, dem, "%ScratchWorkspace%\\Fpath_point_Z", "NONE", "VALUE_ONLY")
    Fpath_point_Z2 = arcpy.MakeFeatureLayer_management(Fpath_point_Z,"Fpath_point_Z.shp")
    inpt =  "{} RASTERVALU PointElevation".format(Fpath_point_Z2)
    outTrend = arcpy.ddd.TopoToRaster(inpt, "outTrend", cellsize)   
    detrended = Minus (dem, outTrend)
    detrended.save(output_folder + "/" + "DED.tif")

    #PART D calculate floodplain zone data properties
    arcpy.AddMessage("STEP 4 Create floodplain zone statistic with channel segments")
    unionFAMseg = arcpy.Intersect_analysis ([union2, segments], "%ScratchWorkspace%\\unionFAMseg", "ALL")

    unionFAMsegSingle = arcpy.management.MultipartToSinglepart(unionFAMseg, "%ScratchWorkspace%\\unionFAMsegSingle")
    arcpy.DeleteField_management(unionFAMsegSingle, "ORIG_FID")

    # HACH creation
    hachTab = arcpy.sa.ZonalStatisticsAsTable(unionFAMsegSingle, "FID", detrended, "hachTab","DATA", "ALL")

    fldlst = ["MIN", "MAX", "RANGE", "MEAN", "STD", "SUM"]
    fieldList = [field.name for field in arcpy.ListFields(hachTab) if field.name in fldlst]
    
    fldlstE = []
    for field in fieldList:
        arcpy.AddField_management(hachTab, "e_" + field, "DOUBLE")
        with arcpy.da.UpdateCursor(hachTab, (field,"e_{}".format(field))) as cursor:
            for row in cursor:
                row[1] = row[0]
                cursor.updateRow(row) 
        arcpy.DeleteField_management(hachTab, field)
        fldlstE.append("e_" + field)

    #DATA UNION
    arcpy.management.JoinField(unionFAMsegSingle, "OBJECTID", hachTab, "OBJECTID",fldlstE)
    name2 = "M4stattistics_hach.shp"
    arcpy.management.CopyFeatures (unionFAMsegSingle,name2)
    arcpy.management.DefineProjection(name2, SR)

    ##########################################
    #### Only FAM data statistics (k = 3) ####
    ##########################################
if k == 3:
    #PART A calculation FAM
    #union all channel layer
    arcpy.AddMessage("STEP 1 Union all channel polygons")
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
        U_layer.append(newName)
    
    for i in range(len(U_layer)):
        fldnames = [field.name for field in arcpy.ListFields(U_layer[i])]
        fi = "y{}".format(year_list[i])
        if fi not in fldnames:
            arcpy.management.AddField(U_layer[i], "y{}".format(year_list[i]), "LONG")
            with arcpy.da.UpdateCursor(U_layer[i], "y{}".format(year_list[i])) as cursor:
                for row in cursor:
                    row[0] = year_list[i]
                    cursor.updateRow(row)
        fields_fam.append("y{}".format(year_list[i]))

    union = arcpy.Union_analysis (U_layer, "%ScratchWorkspace%\\union", "ALL")

    #cleaning fields
    
    fields_to_delete = [field.name for field in arcpy.ListFields(union) if not field.required and field.name not in fields_fam]
    fields_to_delete.pop() 
    for field in fields_to_delete:
        arcpy.DeleteField_management(union, field)

    #calculate floodplain age (FAM)
    arcpy.AddMessage("STEP 2 Create Floodplain Age Map layer (FAM)")
    arcpy.management.AddField(union, "FAM", "LONG")
    fields_fam.append("FAM")

    with arcpy.da.UpdateCursor(union, fields_fam) as cursor:
        for row in cursor:
            maxF=row[0]
            for i in range(len(row)-1):
                if row[i] > maxF:
                    maxF = row[i]
            row[len(row)-1]=maxF  
            cursor.updateRow(row)    
    name = "fam_layer.shp"
    union2 = arcpy.management.MultipartToSinglepart(union, "%ScratchWorkspace%\\union2")
    arcpy.DeleteField_management(union2, "ORIG_FID")
    arcpy.management.CopyFeatures (union2,name)

    fld = ["FAM"]
    fields_to_delete = [field.name for field in arcpy.ListFields(union2) if not field.required and field.name not in fld]
    fields_to_delete.pop() 
    for field in fields_to_delete:
        arcpy.DeleteField_management(union2, field)

    #PART D calculate floodplain zone data properties
    arcpy.AddMessage("STEP 3 Create floodplain zone statistic with channel segments")
    unionFAMseg = arcpy.Intersect_analysis ([union2, segments], "%ScratchWorkspace%\\unionFAMseg", "ALL")

    unionFAMsegSingle = arcpy.management.MultipartToSinglepart(unionFAMseg, "%ScratchWorkspace%\\unionFAMsegSingle")
    arcpy.DeleteField_management(unionFAMsegSingle, "ORIG_FID")

    #DATA UNION
    name2 = "M4stattistics_FAM.shp"
    arcpy.management.CopyFeatures (unionFAMsegSingle,name2)
    arcpy.management.DefineProjection(name2, SR)



#===============================================================================
# DELETING processing FILES
#===============================================================================
if deleteTF == True:
   arcpy.AddMessage("Deleting processing files")
   for i in range(len(U_layer)):
      arcpy.Delete_management(U_layer[i])
else:
   arcpy.AddMessage("Processing files preserved in output folder")


#===============================================================================
# DELETING TEMPORARY FILES
#===============================================================================
if k == 1:
    arcpy.AddMessage("Deleting temporary files")
    arcpy.Delete_management(outTrend)
    arcpy.Delete_management(chm)
    arcpy.Delete_management(hachTab)
    arcpy.Delete_management(vegTab)
    arcpy.Delete_management(union)
    arcpy.Delete_management(union2)
    arcpy.Delete_management(Fpath_dis)
    arcpy.Delete_management(Fpath_point)
    arcpy.Delete_management(Fpath_point_Z)
    arcpy.Delete_management(Fpath_point_Z2)
    arcpy.Delete_management(unionFAMseg)
    arcpy.Delete_management(unionFAMsegSingle)

if k == 2:
    arcpy.AddMessage("Deleting temporary files")
    arcpy.Delete_management(outTrend)
    arcpy.Delete_management(hachTab)
    arcpy.Delete_management(union)
    arcpy.Delete_management(union2)
    arcpy.Delete_management(Fpath_dis)
    arcpy.Delete_management(Fpath_point)
    arcpy.Delete_management(Fpath_point_Z)
    arcpy.Delete_management(Fpath_point_Z2)
    arcpy.Delete_management(unionFAMseg)
    arcpy.Delete_management(unionFAMsegSingle)
    
if k == 3:
    arcpy.AddMessage("Deleting temporary files")
    arcpy.Delete_management(union)
    arcpy.Delete_management(union2)
    arcpy.Delete_management(unionFAMseg)
    arcpy.Delete_management(unionFAMsegSingle)