# tnz
Code for a project examining whether metabolism constrains bird and mammal ranges

Files descriptions and work flow:

Assemble data:
  MRelevation_AssembleDataTraits

Use range shapefiles to extract Tmin and Tmax across range:
  MRelevation_ExtractTminTmax_birds
  MRelevation_ExtractTminTmax_mammals
  (Also MRelevation_SplitMammalShapefiles)
  
Add Tmin and Tmax to dataset:
  MRelevation_AddTminScope
  
Project Ranges:
  MRelevation_ProjectRangesWorldClim_Birds
  MRelevation_ProjectRangesWorldClim_Mammals
