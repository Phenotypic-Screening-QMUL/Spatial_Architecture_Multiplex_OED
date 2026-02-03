# Spatial_Architecture_Multiplex_OED
Analysis for Sheffield/QMUL Project

Manuscript Working Title: : Spatial architecture of the immune microenvironment predicts malignant transformation in oral epithelial dysplasia

Overview:
Scripts are numbered and should be run consecutivley - raw data is available on request

1)	Import Data 
2)	Make maps for all markers and produce scaled intensity data
3)	Stats Using HALO Classifications (Binary / Intensity) 
  a)	Using HALO Classifications (marker positivity) per Annotation Layer
  b)	Using HALO Classifications (marker positivity) per transformation status
  c)	Using HALO Intensities per transformation status
  d)	Using HALO Classifications (marker positivity) per transformation status + Analysis Region 
    i)	Stats and Box Plot 
    ii)	Dot Plot
4)	Cluster Based Cell typing
  a)	Clustering and UMAPs 
  b)	Heatmaps 
5)	Stats using Cluster Based Cell Typing
6)	Find Nearest Neighbour – input for CN analysis 
7)	Neighbourhood Analysis – individual plot maps
8)	Interpreting Neighbourhoods 
9)	CN Stats and Composition 
  a)	  Stats and Heatmap 
  b)	  Dot plot per region 
10)	Distance Changes 
  a)  Stats 
  b)	Interactions 
11)	Models with Cross validation 
12)	Random Forrest PCA
