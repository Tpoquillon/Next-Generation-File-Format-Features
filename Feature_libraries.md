Table of content
- [Basic features](#basic-features)
  * [[Pyradiomics](https://pyradiomics.readthedocs.io/en/latest/features.html)]
    + [First Order Statistics (19 features)](#first-order-statistics--19-features-)
    + [Shape-based (2D, 3D) (26 features)](#shape-based--2d--3d---26-features-)
    + [Gray Level Co-occurrence Matrix (24 features)](#gray-level-co-occurrence-matrix--24-features-)
    + [Gray Level Run Length Matrix (16 features)](#gray-level-run-length-matrix--16-features-)
    + [Gray Level Size Zone Matrix (16 features)](#gray-level-size-zone-matrix--16-features-)
    + [Neighbouring Gray Tone Difference Matrix (5 features)](#neighbouring-gray-tone-difference-matrix--5-features-)
    + [Gray Level Dependence Matrix (16 features)](#gray-level-dependence-matrix--16-features-)
  * [[Nyxus](https://github.com/PolusAI/nyxus)]
    + [Pixel intensity features](#pixel-intensity-features)
    + [Morphology features](#morphology-features)
    + [Texture features](#texture-features)
    + [Radial intensity distribution features](#radial-intensity-distribution-features)
    + [Frequency and orientational features](#frequency-and-orientational-features)
    + [2D image moments](#2d-image-moments)
    + [Neighbor features](#neighbor-features)
  * [[Pyfeats](https://github.com/giakou4/pyfeats)]
    + [Textural Features](#textural-features)
    + [Morphological Features](#morphological-features)
    + [Histogram Based Features](#histogram-based-features)
    + [Multi-scale Features](#multi-scale-features)
    + [Other Features](#other-features)
  * [[FeatureExtraction](https://github.com/ProfBressan/FeatureExtraction)]
    + [BIC (Border/Interior Pixel Classification) 128 features](#--bic--border-interior-pixel-classification----128-features)
    + [TAS (Threshold Adjacency Statistics) 162 features](#--tas--threshold-adjacency-statistics----162-features)
    + [LBP (Local Binary Part) 352 features](#--lbp--local-binary-part----352-features)
    + [FOM (First Order Measures) 8 features (gray) | 24 features (color)](#--fom--first-order-measures----8-features--gray----24-features--color-)
    + [Zernike 72 features](#--zernike---72-features)
    + [Haralick 13 features](#--haralick---13-features)
    + [GCH (Global Color Histogram) 30 features](#--gch--global-color-histogram----30-features)
  * [[Mahotas](https://mahotas.readthedocs.io/en/latest/features.html)]
    + [SURF](#surf)
    + [Zernike](#zernike)
    + [Haralick](#haralick)
    + [LBP](#lbp)
  * [[scMultiplex](https://github.com/fmi-basel/gliberal-scMultipleX/tree/main)]
    + [Standard Pixel Features](#standard-pixel-features)
    + [Standard Object Features](#standard-object-features)
    + [Convex Hull Features](#convex-hull-features)
    + [Skeleton Features](#skeleton-features)
  * [[Scikit-image](https://scikit-image.org/docs/dev/api/skimage.feature.html)]
    + [Features](#features)
    + [Filters](#filters)
  * [[Xrayimage](https://github.com/vatsalsaglani/xrayimage_extractfeatures/tree/master)]
    + [Entropy features](#entropy-features)
    + [GLCM Features](#glcm-features)
    + [Moments](#moments)
    + [Region Properties](#region-properties)
- [Spatial Features](#spatial-features)
  * [[SquidPy](https://squidpy.readthedocs.io/en/stable/)]
    + [Graph features](#graph-features)
    + [Image features](#image-features)
    + [Spatial plots](#spatial-plots)
  * [[DypFISH](https://github.com/cbib/dypfish)]
    + [Spatial features](#spatial-features)
  * [[Giotto](https://giottosuite.readthedocs.io/en/master/documentation.html#cell-neighborhood)]
    + [Cell neighborhood](#cell-neighborhood)
    + [Cell cell](#cell-cell)
    + [Cell cell communication](#cell-cell-communication)
  * [[Big-FISH](https://big-fish.readthedocs.io/en/stable/classification/features.html)]
  * [[IMCDataAnalysis (R)](https://bodenmillergroup.github.io/IMCDataAnalysis/performing-spatial-analysis.html)]
    + [Spatial interaction graphs](#spatial-interaction-graphs)
    + [Spatial community analysis](#spatial-community-analysis)
    + [Cellular neighborhood analysis](#cellular-neighborhood-analysis)
    + [Spatial context analysis](#spatial-context-analysis)
    + [Patch detection](#patch-detection)
    + [Interaction analysis](#interaction-analysis)
  * [[SpatialEpiApp (R)](https://github.com/Paula-Moraga/SpatialEpiApp)]
    + [Areal](#areal)
    + [Geostatistical features](#geostatistical-features)
    + [Spatial point patterns](#spatial-point-patterns)
- [Graph features](#graph-features-1)
  * [[NetworkX](https://networkx.org/documentation/stable/index.html)]
- [Other resources](#other-resources)

## Basic features
### [Pyradiomics](https://pyradiomics.readthedocs.io/en/latest/features.html)

#### First Order Statistics (19 features)
```
getEnergyFeatureValue
getTotalEnergyFeatureValue
getEntropyFeatureValue
getMinimumFeatureValue
get10PercentileFeatureValue
get90PercentileFeatureValue
getMaximumFeatureValue
getMeanFeatureValue
getMedianFeatureValue
getInterquartileRangeFeatureValue
getRangeFeatureValue
getMeanAbsoluteDeviationFeatureValue
getRobustMeanAbsoluteDeviationFeatureValue
getRootMeanSquaredFeatureValue
getStandardDeviationFeatureValue
getSkewnessFeatureValue
getKurtosisFeatureValue
getVarianceFeatureValue
getUniformityFeatureValue
```
#### Shape-based (2D, 3D) (26 features)
```
getMeshSurfaceFeatureValue
getPixelSurfaceFeatureValue
getPerimeterFeatureValue
getPerimeterSurfaceRatioFeatureValue
getSphericityFeatureValue
getSphericalDisproportionFeatureValue
getMaximumDiameterFeatureValue
getMajorAxisLengthFeatureValue
getMinorAxisLengthFeatureValue
getElongationFeatureValue
getMeshVolumeFeatureValue
getVoxelVolumeFeatureValue
getSurfaceAreaFeatureValue
getSurfaceVolumeRatioFeatureValue
getSphericityFeatureValue
getCompactness1FeatureValue
getCompactness2FeatureValue
getSphericalDisproportionFeatureValue
getMaximum3DDiameterFeatureValue
getMaximum2DDiameterSliceFeatureValue
getMaximum2DDiameterColumnFeatureValue
getMaximum2DDiameterRowFeatureValue
getMajorAxisLengthFeatureValue
getMinorAxisLengthFeatureValue
getLeastAxisLengthFeatureValue
getElongationFeatureValue
getFlatnessFeatureValue
```

#### Gray Level Co-occurrence Matrix (24 features)
```
getAutocorrelationFeatureValue
getJointAverageFeatureValue
getClusterProminenceFeatureValue
getClusterShadeFeatureValue
getClusterTendencyFeatureValue
getContrastFeatureValue
getCorrelationFeatureValue
getDifferenceAverageFeatureValue
getDifferenceEntropyFeatureValue
getDifferenceVarianceFeatureValue
getDissimilarityFeatureValue
getJointEnergyFeatureValue
getJointEntropyFeatureValue
getHomogeneity1FeatureValue
getHomogeneity2FeatureValue
getImc1FeatureValue
getImc2FeatureValue
getIdmFeatureValue
getMCCFeatureValue
getIdmnFeatureValue
getIdFeatureValue
getIdnFeatureValue
getInverseVarianceFeatureValue
getMaximumProbabilityFeatureValue
getSumAverageFeatureValue
getSumVarianceFeatureValue
getSumEntropyFeatureValue
getSumSquaresFeatureValue
```
#### Gray Level Run Length Matrix (16 features)
```
getShortRunEmphasisFeatureValue
getLongRunEmphasisFeatureValue
getGrayLevelNonUniformityFeatureValue
getGrayLevelNonUniformityNormalizedFeatureValue
getRunLengthNonUniformityFeatureValue
getRunLengthNonUniformityNormalizedFeatureValue
getRunPercentageFeatureValue
getGrayLevelVarianceFeatureValue
getRunVarianceFeatureValue
getRunEntropyFeatureValue
getLowGrayLevelRunEmphasisFeatureValue
getHighGrayLevelRunEmphasisFeatureValue
getShortRunLowGrayLevelEmphasisFeatureValue
getShortRunHighGrayLevelEmphasisFeatureValue
getLongRunLowGrayLevelEmphasisFeatureValue
getLongRunHighGrayLevelEmphasisFeatureValue
```
#### Gray Level Size Zone Matrix (16 features)
```
getSmallAreaEmphasisFeatureValue
getLargeAreaEmphasisFeatureValue
getGrayLevelNonUniformityFeatureValue
getGrayLevelNonUniformityNormalizedFeatureValue
getSizeZoneNonUniformityFeatureValue
getSizeZoneNonUniformityNormalizedFeatureValue
getZonePercentageFeatureValue
getGrayLevelVarianceFeatureValue
getZoneVarianceFeatureValue
getZoneEntropyFeatureValue
getLowGrayLevelZoneEmphasisFeatureValue
getHighGrayLevelZoneEmphasisFeatureValue
getSmallAreaLowGrayLevelEmphasisFeatureValue
getSmallAreaHighGrayLevelEmphasisFeatureValue
getLargeAreaLowGrayLevelEmphasisFeatureValue
getLargeAreaHighGrayLevelEmphasisFeatureValue
```
#### Neighbouring Gray Tone Difference Matrix (5 features)
```
getCoarsenessFeatureValue
getContrastFeatureValue
getBusynessFeatureValue
getComplexityFeatureValue
getStrengthFeatureValue
```
#### Gray Level Dependence Matrix (16 features)
```
getSmallDependenceEmphasisFeatureValue
getLargeDependenceEmphasisFeatureValue
getGrayLevelNonUniformityFeatureValue
getGrayLevelNonUniformityNormalizedFeatureValue
getDependenceNonUniformityFeatureValue
getDependenceNonUniformityNormalizedFeatureValue
getGrayLevelVarianceFeatureValue
getDependenceVarianceFeatureValue
getDependenceEntropyFeatureValue
getDependencePercentageFeatureValue
getLowGrayLevelEmphasisFeatureValue
getHighGrayLevelEmphasisFeatureValue
getSmallDependenceLowGrayLevelEmphasisFeatureValue
getSmallDependenceHighGrayLevelEmphasisFeatureValue
getLargeDependenceLowGrayLevelEmphasisFeatureValue
getLargeDependenceHighGrayLevelEmphasisFeatureValue
```
### [Nyxus](https://github.com/PolusAI/nyxus)
Complete feature list : [https://github.com/PolusAI/nyxus/blob/main/docs/source/featurelist.rst](https://github.com/PolusAI/nyxus/blob/main/docs/source/featurelist.rst)

#### Pixel intensity features
```
INTEGRATED\_INTENSITY : Integrated intensity of the region of interest (ROI)
MEAN : Mean intensity value of the ROI
MEDIAN : The median value of pixels in the ROI
MIN : Minimum intensity value in the ROI
MAX : Maximum intensity value in the ROI
RANGE : Range between the maximmu and minimum
COVERED\_IMAGE\_INTENSITY\_RANGE : intensity range of the ROI to intensity range of all the ROIs
STANDARD\_DEVIATION : Standard deviation (unbiased)
STANDARD\_DEVIATION\_BIASED : Biased standard deviation
COV : Coefficient of variation
STANDARD\_ERROR : Standard error
SKEWNESS : Skewness - the 3rd standardized moment
KURTOSIS : Kurtosis - the 4th standardized moment (Pearson formula)
EXCESS\_KURTOSIS : Excess kurtosis - the 4th standardized moment (Fisher-corrected formula, IBSI feature IPH6)
HYPERSKEWNESS : Hyperskewness - the 5th standardized moment
HYPERFLATNESS : Hyperflatness - the 6th standardized moment
MEAN\_ABSOLUTE\_DEVIATION : Mean absolute deviation
MEDIAN\_ABSOLUTE\_DEVIATION : Median absolute deviation
ENERGY : ROI energy
ROOT\_MEAN\_SQUARED : Root of mean squared deviation
ENTROPY : ROI entropy - a measure of the amount of information (that is, randomness) in the ROI
MODE : The mode value of pixels in the ROI - the value that appears most often in a set of ROI intensity values
UNIFORMITY : Uniformity - measures how uniform the distribution of ROI intensities is
UNIFORMITY\_PIU : Percent image uniformity, another measure of intensity distribution uniformity
P01, P10, P25, P75, P90, P99 : 1%, 10%, 25%, 75%, 90%, and 99% percentiles of intensity distribution
QCOD : quantile coefficient of dispersion
INTERQUARTILE\_RANGE : Distribution's interquartile range
ROBUST\_MEAN\_ABSOLUTE\_DEVIATION : Robust mean absolute deviation
MASS\_DISPLACEMENT : ROI mass displacement
INTEGRATED\_INTENSITY : Integrated intensity of the region of interest (ROI)
MEAN : Mean intensity value of the ROI
MEDIAN : The median value of pixels in the ROI
MIN : Minimum intensity value in the ROI
MAX : Maximum intensity value in the ROI
RANGE : Range between the maximmu and minimum
COVERED\_IMAGE\_INTENSITY\_RANGE : intensity range of the ROI to intensity range of all the ROIs
STANDARD\_DEVIATION : Standard deviation (unbiased)
STANDARD\_DEVIATION\_BIASED : Biased standard deviation
COV : Coefficient of variation
STANDARD\_ERROR : Standard error
SKEWNESS : Skewness - the 3rd standardized moment
KURTOSIS : Kurtosis - the 4th standardized moment (Pearson formula)
EXCESS\_KURTOSIS : Excess kurtosis - the 4th standardized moment (Fisher-corrected formula, IBSI feature IPH6)
HYPERSKEWNESS : Hyperskewness - the 5th standardized moment
HYPERFLATNESS : Hyperflatness - the 6th standardized moment
MEAN\_ABSOLUTE\_DEVIATION : Mean absolute deviation
MEDIAN\_ABSOLUTE\_DEVIATION : Median absolute deviation
ENERGY : ROI energy
ROOT\_MEAN\_SQUARED : Root of mean squared deviation
ENTROPY : ROI entropy - a measure of the amount of information (that is, randomness) in the ROI
MODE : The mode value of pixels in the ROI - the value that appears most often in a set of ROI intensity values
UNIFORMITY : Uniformity - measures how uniform the distribution of ROI intensities is
UNIFORMITY\_PIU : Percent image uniformity, another measure of intensity distribution uniformity
P01, P10, P25, P75, P90, P99 : 1%, 10%, 25%, 75%, 90%, and 99% percentiles of intensity distribution
QCOD : quantile coefficient of dispersion
INTERQUARTILE\_RANGE : Distribution's interquartile range
ROBUST\_MEAN\_ABSOLUTE\_DEVIATION : Robust mean absolute deviation
MASS\_DISPLACEMENT : ROI mass displacement
```
#### Morphology features
```
AREA\_PIXELS\_COUNT : ROI area in the number of pixels
AREA\_UM2 : ROI area in metric units
CENTROID\_X : X-coordinate of the enter point of the ROI
CENTROID\_Y : Y-coordinate of the center point of the ROI
COMPACTNESS : Mean squared distance of the object's pixels from the centroid divided by the area. Compactness of a filled circle is 1, compactness of irregular objects or objects with holes is greater than 1
BBOX\_YMIN : Y-position and size of the smallest axis-aligned box containing the ROI
BBOX\_XMIN : X-position and size of the smallest axis-aligned box containing the ROI
BBOX\_HEIGHT : Height of the smallest axis-aligned box containing the ROI
BBOX\_WIDTH : Width of the smallest axis-aligned box containing the ROI
MAJOR\_AXIS\_LENGTH : Length (in pixels) of the major axis of the ellipse that has the same normalized second central moments as the region
MINOR\_AXIS\_LENGTH : Length (in pixels) of the minor axis of the ellipse that has the same normalized second central moments as the region
ECCENTRICITY : Ratio of ROI's inertia ellipse focal distance over the major axis length
ORIENTATION : Angle between the 0th axis and the major axis of the ellipse that has same second moments as the region
ROUNDNESS : Represents how similar a ROI's inertia ellipse is to circle. Calculated based on the major and minor exis lengths
EXTENT : Proportion of the pixels (2D) or voxels (3D) in the bounding box that are also in the region. Computed as the area/volume of the object divided by the area/volume of the bounding box
ASPECT\_RATIO : The ratio of the major axis to the minor axis of ROI's inertia ellipse
CONVEX\_HULL\_AREA : Area of ROI's convex hull
SOLIDITY : Ratio of pixels in the ROI common with its convex hull image
PERIMETER : Number of pixels in ROI's contour
EQUIVALENT\_DIAMETER : Diameter of a circle with the same area as the ROI
EDGE\_MEAN\_INTENSITY : Mean intensity of ROI's contour pixels
EDGE\_STDDEV\_INTENSITY : Standard deviation of ROI's contour pixels
EDGE\_MAX\_INTENSITY : Maximum intensity of ROI's contour pixels
EDGE\_MIN\_INTENSITY : Minimum intensity of ROI's contour pixels
CIRCULARITY : Represents how similar a shape is to circle. Clculated based on ROI's area and its convex perimeter
EROSIONS\_2\_VANISH : Number of erosion operations for a ROI to vanish in its axis aligned bounding box
EROSIONS\_2\_VANISH\_COMPLEMENT : Number of erosion operations for a ROI to vanish in its convex hull
FRACT\_DIM\_BOXCOUNT : Fractal dimension determined by the box counting method according to ISO 9276-6. If C is a fractal set, with fractal dimension DF \< D, then the number N of boxes of size R needed to cover the set scales as R^(-DF). DF is known as the Hausdorf dimension, or Kolmogorov capacity, or Kolmogorov dimension, or simply box-counting dimension
FRACT\_DIM\_PERIMETER : Fractal dimension determined by the perimeter method according to ISO 9276-6. If we approximate ROI's contour with rulers of length lambda, the perimeter based fractal dimension is the slope of the best fit line of log ROI perimeter versus log lambda, subtracted from 1
WEIGHTED\_CENTROID\_Y : X-coordinate of centroid
WEIGHTED\_CENTROID\_X : Y-coordinate of centroid
MIN\_FERET\_DIAMETER : Feret diameter (or maximum caliber diameter) is the longest distance between any two ROI points along the same (horizontal) direction. This feature is the minimum Feret diameter for angles ranging 0 to 180 degrees
MAX\_FERET\_DIAMETER : Maximum Feret diameter for angles ranging 0 to 180 degrees
MIN\_FERET\_ANGLE : Angle of the minimum Feret diameter
MAX\_FERET\_ANGLE : Angle of the maximum Feret diameter
STAT\_FERET\_DIAM\_MIN : Minimum of Feret diameters of the ROI rotated at angles 0-180 degrees
STAT\_FERET\_DIAM\_MAX : Maximum of Feret diameters of the ROI rotated at angles 0-180 degrees
STAT\_FERET\_DIAM\_MEAN : Mean Feret diameter of the ROI rotated at angles 0-180 degrees
STAT\_FERET\_DIAM\_MEDIAN : Median value of Feret diameters of the ROI rotated at angles 0-180 degrees
STAT\_FERET\_DIAM\_STDDEV : Standard deviation of Feret diameter of the ROI rotated at angles 0-180 degrees
STAT\_FERET\_DIAM\_MODE : Histogram mode of Feret diameters of the ROI rotated at angles 0-180 degrees
STAT\_MARTIN\_DIAM\_MIN : Minimum of Martin diameters of the ROI rotated at angles 0-180 degrees
STAT\_MARTIN\_DIAM\_MAX : Maximum of Martin diameters of the ROI rotated at angles 0-180 degrees
STAT\_MARTIN\_DIAM\_MEAN : Mean of Martin diameter of the ROI rotated at angles 0-180 degrees
STAT\_MARTIN\_DIAM\_MEDIAN : Median value of Martin diameters of the ROI rotated at angles 0-180 degrees
STAT\_MARTIN\_DIAM\_STDDEV : Standard deviation of Martin diameter of the ROI rotated at angles 0-180 degrees
STAT\_MARTIN\_DIAM\_MODE : Histogram mode of Martin diameters of the ROI rotated at angles 0-180 degrees
STAT\_NASSENSTEIN\_DIAM\_MIN : Minimum of Nassenstein diameters of the ROI rotated at angles 0-180 degrees
STAT\_NASSENSTEIN\_DIAM\_MAX : Maximum of Nassenstein diameters of the ROI rotated at angles 0-180 degrees
STAT\_NASSENSTEIN\_DIAM\_MEAN : Mean of Nassenstein diameter of the ROI rotated at angles 0-180 degrees
STAT\_NASSENSTEIN\_DIAM\_MEDIAN : Median value of Nassenstein diameters of the ROI rotated at angles 0-180 degrees
STAT\_NASSENSTEIN\_DIAM\_STDDEV : Standard deviation of Nassenstein diameter of the ROI rotated at angles 0-180 degrees
STAT\_NASSENSTEIN\_DIAM\_MODE : Histogram mode of Nassenstein diameters of the ROI rotated at angles 0-180 degrees
MAXCHORDS\_MAX : Maximum of ROI's longest chords built at angles 0-180 degrees
MAXCHORDS\_MAX\_ANG : Angle of the chord referenced in MAXCHORDS\_MAX
MAXCHORDS\_MIN : Minimum of ROI's longest chords built at angles 0-180 degrees
MAXCHORDS\_MIN\_ANG : Angle of the chord referenced in MAXCHORDS\_MIN
MAXCHORDS\_MEDIAN : Median value of ROI's longest chords built at angles 0-180 degrees
MAXCHORDS\_MEAN : Mean value of ROI's longest chords built at angles 0-180 degrees
MAXCHORDS\_MODE : Histogram mode of ROI's longest chords built at angles 0-180 degrees
MAXCHORDS\_STDDEV : Sndard deviation of ROI's longest chords built at angles 0-180 degrees
ALLCHORDS\_MAX : Maximum of all the ROI's chords built at angles 0-180 degrees
ALLCHORDS\_MAX\_ANG : Angle of the chord referenced in ALLCHORDS\_MAX
ALLCHORDS\_MIN : Minimum of all the ROI's chords built at angles 0-180 degrees
ALLCHORDS\_MIN\_ANG : Angle of the chord referenced in ALLCHORDS\_MIN
ALLCHORDS\_MEDIAN : Median value of all the ROI's chords built at angles 0-180 degrees
ALLCHORDS\_MEAN : Mean value of all the ROI's chords built at angles 0-180 degrees
ALLCHORDS\_MODE : Histogram mode of all the ROI's chords built at angles 0-180 degrees
ALLCHORDS\_STDDEV : Sndard deviation of all the ROI's chords built at angles 0-180 degrees
EULER\_NUMBER : Euler characteristic of the ROI - the number of objects in the ROI minus the number of holes assuming the 8-neighbor connectivity of ROI's pixels
EXTREMA\_P1\_X : X-ccordinate of ROI's axis aligned bounding box extremum point #1
EXTREMA\_P1\_Y : Y-ccordinate of ROI's axis aligned bounding box extremum point #1
EXTREMA\_P2\_X : X-ccordinate of ROI's axis aligned bounding box extremum point #2
EXTREMA\_P2\_Y :
EXTREMA\_P3\_X : X-ccordinate of ROI's axis aligned bounding box extremum point #3
EXTREMA\_P3\_Y :
EXTREMA\_P4\_X : X-ccordinate of ROI's axis aligned bounding box extremum point #4
EXTREMA\_P4\_Y :
EXTREMA\_P5\_X : X-ccordinate of ROI's axis aligned bounding box extremum point #5
EXTREMA\_P5\_Y :
EXTREMA\_P6\_X : X-ccordinate of ROI's axis aligned bounding box extremum point #6
EXTREMA\_P6\_Y :
EXTREMA\_P7\_X : X-ccordinate of ROI's axis aligned bounding box extremum point #7
EXTREMA\_P7\_Y :
EXTREMA\_P8\_X : X-ccordinate of ROI's axis aligned bounding box extremum point #8
EXTREMA\_P8\_Y :
POLYGONALITY\_AVE : The score ranges from $ -infty $ to 10. Score 10 indicates the object shape is polygon and score $ -infty $ indicates the ROI shape is not polygon
HEXAGONALITY\_AVE : The score ranges from $ -infty $ to 10. Score 10 indicates the object shape is hexagon and score $ -infty $ indicates the ROI shape is not hexagon
HEXAGONALITY\_STDDEV : Standard deviation of hexagonality\_score relative to its mean
DIAMETER\_MIN\_ENCLOSING\_CIRCLE : Diameter of the minimum enclosing circle
DIAMETER\_CIRCUMSCRIBING\_CIRCLE : Diameter of the circumscribing circle
DIAMETER\_INSCRIBING\_CIRCLE : Diameter of inscribing circle
GEODETIC\_LENGTH : Geodetic length approximated by a rectangle with the same area and perimeter: $ area = geodeticlength \* thickness$;
THICKNESS : Thickness approximated by a rectangle with the same area and perimeter: $ area = geodeticlength \* thickness$;
ROI\_RADIUS\_MEAN : Mean centroid to edge distance
ROI\_RADIUS\_MAX : Maximum of centroid to edge distances
ROI\_RADIUS\_MEDIAN : Median value of centroid to edge distances
```
#### Texture features
```
GLCM\_ASM : GLCM, Angular second moment, IBSI # 8ZQL
GLCM\_ACOR : GLCM, Autocorrelation, IBSI # QWB0
GLCM\_CLUPROM : GLCM, Cluster prominence, IBSI # AE86
GLCM\_CLUSHADE : GLCM, Cluster shade, IBSI # 7NFM
GLCM\_CLUTEND : GLCM, Cluster tendency, IBSI # DG8W
GLCM\_CONTRAST : GLCM, Contrast, IBSI # ACUI
GLCM\_CORRELATION : GLCM, Correlation, IBSI # NI2N
GLCM\_DIFAVE : GLCM, Difference average, IBSI # TF7R
GLCM\_DIFENTRO : GLCM, Difference entropy, IBSI # NTRS
GLCM\_DIFVAR : GLCM, Difference variance, IBSI # D3YU
GLCM\_DIS : GLCM, Dissimilarity, IBSI # 8S9J
GLCM\_ENERGY : GLCM, Energy
GLCM\_ENTROPY : GLCM, Entropy
GLCM\_HOM1 : GLCM, Homogeneity-1
GLCM\_HOM2 : GLCM, Homogeneity-2
GLCM\_ID : GLCM, Inverse difference, IBSI # IB1Z
GLCM\_IDN : GLCM, Inverse difference normalized, IBSI # NDRX
GLCM\_IDM : GLCM, Inverse difference moment, IBSI # WF0Z
GLCM\_IDMN : GLCM, Inverse difference moment normalized, IBSI # 1QCO
GLCM\_INFOMEAS1 : GLCM, Information measure of correlation 1, IBSI # R8DG
GLCM\_INFOMEAS2 : GLCM, Information measure of correlation 2, IBSI # JN9H
GLCM\_IV : GLCM, Inverse variance, IBSI # E8JP
GLCM\_JAVE : GLCM, Joint average, IBSI # 60VM
GLCM\_JE : GLCM, Joint entropy, IBSI # TU9B
GLCM\_JMAX : GLCM, Joint maximum (aka max probability), IBSI # GYBY
GLCM\_JVAR : GLCM, Joint variance (aka sum of squares), IBSI # UR99
GLCM\_SUMAVERAGE : GLCM, Sum average, IBSI # ZGXS
GLCM\_SUMENTROPY : GLCM, Sum entropy, IBSI # P6QZ
GLCM\_SUMVARIANCE : GLCM, Sum variance, IBSI # OEEB
GLCM\_VARIANCE : GLCM, Variance
GLRLM\_SRE : Grey level run-length matrix (GLRLM) based feature, Short Run Emphasis
GLRLM\_LRE : GLRLM, Long Run Emphasis
GLRLM\_GLN : GLRLM, Grey Level Non-Uniformity
GLRLM\_GLNN : GLRLM, Grey Level Non-Uniformity Normalized
GLRLM\_RLN : GLRLM, Run Length Non-Uniformity
GLRLM\_RLNN : GLRLM, Run Length Non-Uniformity Normalized
GLRLM\_RP : GLRLM, Run Percentage
GLRLM\_GLV : GLRLM, Grey Level Variance
GLRLM\_RV : GLRLM, Run Variance
GLRLM\_RE : GLRLM, Run Entropy
GLRLM\_LGLRE : GLRLM, Low Grey Level Run Emphasis
GLRLM\_HGLRE : GLRLM, High Grey Level Run Emphasis
GLRLM\_SRLGLE : GLRLM, Short Run Low Grey Level Emphasis
GLRLM\_SRHGLE : GLRLM, Short Run High Grey Level Emphasis
GLRLM\_LRLGLE : GLRLM, Long Run Low Grey Level Emphasis
GLRLM\_LRHGLE : GLRLM, Long Run High Grey Level Emphasis
GLDZM\_SDE : GLDZM, Small Distance Emphasis
GLDZM\_LDE : GLDZM, Large Distance Emphasis
GLDZM\_LGLE : GLDZM, Low Grey Level Emphasis
GLDZM\_HGLE : GLDZM, High GreyLevel Emphasis
GLDZM\_SDLGLE : GLDZM, Small Distance Low Grey Level Emphasis
GLDZM\_SDHGLE : GLDZM, Small Distance High GreyLevel Emphasis
GLDZM\_LDLGLE : GLDZM, Large Distance Low Grey Level Emphasis
GLDZM\_LDHGLE : GLDZM, Large Distance High Grey Level Emphasis
GLDZM\_GLNU : GLDZM, Grey Level Non Uniformity
GLDZM\_GLNUN : GLDZM, Grey Level Non Uniformity Normalized
GLDZM\_ZDNU : GLDZM, Zone Distance Non Uniformity
GLDZM\_ZDNUN : GLDZM, Zone Distance Non Uniformity Normalized
GLDZM\_ZP : GLDZM, Zone Percentage
GLDZM\_GLM : GLDZM, Grey Level Mean
GLDZM\_GLV : GLDZM, Grey Level Variance
GLDZM\_ZDM : GLDZM, Zone Distance Mean
GLDZM\_ZDV : GLDZM, Zone Distance Variance
GLDZM\_ZDE : GLDZM, Zone Distance Entropy
GLSZM\_SAE : GLDZM, Grey level size zone matrix (GLSZM) based feature, Small Area Emphasis
GLSZM\_LAE : Large Area Emphasis
GLSZM\_GLN : Grey Level Non - Uniformity
GLSZM\_GLNN : Grey Level Non - Uniformity Normalized
GLSZM\_SZN : Size - Zone Non - Uniformity
GLSZM\_SZNN : Size - Zone Non - Uniformity Normalized
GLSZM\_ZP : Zone Percentage
GLSZM\_GLV : Grey Level Variance
GLSZM\_ZV : Zone Variance
GLSZM\_ZE : Zone Entropy
GLSZM\_LGLZE : Low Grey Level Zone Emphasis
GLSZM\_HGLZE : High Grey Level Zone Emphasis
GLSZM\_SALGLE : Small Area Low Grey Level Emphasis
GLSZM\_SAHGLE : Small Area High Grey Level Emphasis
GLSZM\_LALGLE : Large Area Low Grey Level Emphasis
GLSZM\_LAHGLE : Large Area High Grey Level Emphasis
GLDM\_SDE : Grey level dependency matrix (GLDM) based feature, Small Dependence Emphasis(SDE)
GLDM\_LDE : Large Dependence Emphasis (LDE)
GLDM\_GLN : Grey Level Non-Uniformity (GLN)
GLDM\_DN : Dependence Non-Uniformity (DN)
GLDM\_DNN : Dependence Non-Uniformity Normalized (DNN)
GLDM\_GLV : Grey Level Variance (GLV)
GLDM\_DV : Dependence Variance (DV)
GLDM\_DE : Dependence Entropy (DE)
GLDM\_LGLE : Low Grey Level Emphasis (LGLE)
GLDM\_HGLE : High Grey Level Emphasis (HGLE)
GLDM\_SDLGLE : Small Dependence Low Grey Level Emphasis (SDLGLE)
GLDM\_SDHGLE : Small Dependence High Grey Level Emphasis (SDHGLE)
GLDM\_LDLGLE : Large Dependence Low Grey Level Emphasis (LDLGLE)
GLDM\_LDHGLE : Large Dependence High Grey Level Emphasis (LDHGLE)
NGLDM\_LDE : Low Dependence Emphasis
NGLDM\_HDE : High Dependence Emphasis
NGLDM\_LGLCE : Low Grey Level Count Emphasis
NGLDM\_HGLCE : High Grey Level Count Emphasis
NGLDM\_LDLGLE : Low Dependence Low Grey Level Emphasis
NGLDM\_LDHGLE : Low Dependence High Grey Level Emphasis
NGLDM\_HDLGLE : High Dependence Low Grey Level Emphasis
NGLDM\_HDHGLE : High Dependence High Grey Level Emphasis
NGLDM\_GLNU : Grey Level Non-Uniformity
NGLDM\_GLNUN : Grey Level Non-Uniformity Normalised
NGLDM\_DCNU : Dependence Count Non-Uniformity
NGLDM\_DCNUN : Dependence Count Non-Uniformity Normalised
NGLDM\_GLM : Grey Level Mean
NGLDM\_GLV : Grey Level Variance
NGLDM\_DCM : Dependence Count Mean
NGLDM\_DCV : Dependence Count Variance
NGLDM\_DCE : Dependence Count Entropy
NGLDM\_DCENE : Dependence Count Energy
NGTDM\_COARSENESS : Neighbouring Grey Tone Difference Matrix (NGTDM) Features, Coarseness
NGTDM\_CONTRAST : NGTDM, Contrast
NGTDM\_BUSYNESS : NGTDM, Busyness
NGTDM\_COMPLEXITY : NGTDM, Complexity
NGTDM\_STRENGTH : NGTDM, Strength
```
#### Radial intensity distribution features
```
ZERNIKE2D : Zernike features
FRAC\_AT\_D : Fraction of total intensity at a given radius
MEAN\_FRAC : Mean fractional intensity at a given radius
RADIAL\_CV : Coefficient of variation of intensity within a ring (band), calculated across
ZERNIKE2D : Zernike features
FRAC\_AT\_D : Fraction of total intensity at a given radius
MEAN\_FRAC : Mean fractional intensity at a given radius
RADIAL\_CV : Coefficient of variation of intensity within a ring (band), calculated across n slices
```
#### Frequency and orientational features
```
GABOR : A set of Gabor filters of varying frequencies and orientations
```
#### 2D image moments
```
SPAT\_MOMENT\_00 : Spatial (raw) moments
SPAT\_MOMENT\_01 : of order 00, 01, 02, etc
SPAT\_MOMENT\_02 :
SPAT\_MOMENT\_03 :
SPAT\_MOMENT\_10 :
SPAT\_MOMENT\_11 :
SPAT\_MOMENT\_12 :
SPAT\_MOMENT\_20 :
SPAT\_MOMENT\_21 :
SPAT\_MOMENT\_30 :
WEIGHTED\_SPAT\_MOMENT\_00 : Spatial moments weighted by pixel distance to ROI edge
WEIGHTED\_SPAT\_MOMENT\_01 :
WEIGHTED\_SPAT\_MOMENT\_02 :
WEIGHTED\_SPAT\_MOMENT\_03 :
WEIGHTED\_SPAT\_MOMENT\_10 :
WEIGHTED\_SPAT\_MOMENT\_11 :
WEIGHTED\_SPAT\_MOMENT\_12 :
WEIGHTED\_SPAT\_MOMENT\_20 :
WEIGHTED\_SPAT\_MOMENT\_21 :
WEIGHTED\_SPAT\_MOMENT\_30 :
CENTRAL\_MOMENT\_02 : Central moments
CENTRAL\_MOMENT\_03 :
CENTRAL\_MOMENT\_11 :
CENTRAL\_MOMENT\_12 :
CENTRAL\_MOMENT\_20 :
CENTRAL\_MOMENT\_21 :
CENTRAL\_MOMENT\_30 :
WEIGHTED\_CENTRAL\_MOMENT\_02 : Central moments weighted by pixel distance to ROI edge
WEIGHTED\_CENTRAL\_MOMENT\_03 :
WEIGHTED\_CENTRAL\_MOMENT\_11 :
WEIGHTED\_CENTRAL\_MOMENT\_12 :
WEIGHTED\_CENTRAL\_MOMENT\_20 :
WEIGHTED\_CENTRAL\_MOMENT\_21 :
WEIGHTED\_CENTRAL\_MOMENT\_30 :
NORM\_CENTRAL\_MOMENT\_02 : Normalized central moments
NORM\_CENTRAL\_MOMENT\_03 :
NORM\_CENTRAL\_MOMENT\_11 :
NORM\_CENTRAL\_MOMENT\_12 :
NORM\_CENTRAL\_MOMENT\_20 :
NORM\_CENTRAL\_MOMENT\_21 :
NORM\_CENTRAL\_MOMENT\_30 :
NORM\_SPAT\_MOMENT\_00 : Normalized (standardized) spatial moments
NORM\_SPAT\_MOMENT\_01 :
NORM\_SPAT\_MOMENT\_02 :
NORM\_SPAT\_MOMENT\_03 :
NORM\_SPAT\_MOMENT\_10 :
NORM\_SPAT\_MOMENT\_20 :
NORM\_SPAT\_MOMENT\_30 :
HU\_M1 : Hu's moment 1
HU\_M2 : Hu's moment 2
HU\_M3 : Hu's moment 3
HU\_M4 : Hu's moment 4
HU\_M5 : Hu's moment 5
HU\_M6 : Hu's moment 6
HU\_M7 : Hu's moment 7
WEIGHTED\_HU\_M1 : Weighted Hu's moment 1
WEIGHTED\_HU\_M2 : Weighted Hu's moment 2
WEIGHTED\_HU\_M3 : Weighted Hu's moment 3
WEIGHTED\_HU\_M4 : Weighted Hu's moment 4
WEIGHTED\_HU\_M5 : Weighted Hu's moment 5
WEIGHTED\_HU\_M6 : Weighted Hu's moment 6
WEIGHTED\_HU\_M7 : Weighted Hu's moment 7
```
#### Neighbor features
```
NUM\_NEIGHBORS : The number of neighbors bordering the ROI's perimeter within proximity radius specified by command line argument --pixelDistance. (Default value of --pixelDistance is 5.) Algorithmically calculating this feature invilves solving the nearest neighbors search problem that in turn involves the proximity measure and the proximity threshold. Particularly, this plugin uses the L\_2 norm measure over Cartesian space of pixel coordinates and parameter --pixelDistance
PERCENT\_TOUCHING : Percent of ROI's contour pixels located at distance 0 from neighboring other ROIs's contour
CLOSEST\_NEIGHBOR1\_DIST : Distance in pixels from ROI's centroid to the nearest neighboring ROI's centroid
CLOSEST\_NEIGHBOR1\_ANG : Angle in degrees between ROI's centroid and its nearest neighboring ROI's centroid
CLOSEST\_NEIGHBOR2\_DIST : Distance in pixels from ROI's centroid to the second nearest neighboring ROI's centroid
CLOSEST\_NEIGHBOR2\_ANG : Angle in degrees between ROI's centroid and its second nearest neighboring ROI's centroid
ANG\_BW\_NEIGHBORS\_MEAN : Mean angle in degrees between ROI's centroid and centroids of its neighboring ROIs
ANG\_BW\_NEIGHBORS\_STDDEV : Standard deviation in degrees of angles between ROI's centroid and centroids of its neighboring ROIs
ANG\_BW\_NEIGHBORS\_MODE : Mode value in degrees of angles between ROI's centroid and centroids of its neighboring ROIs
```

### [Pyfeats](https://github.com/giakou4/pyfeats)

#### Textural Features
```
First Order Statistics/Statistical Features (FOS/SF)
Gray Level Co-occurence Matrix (GLCM/SGLDM)
Gray Level Difference Statistics (GLDS)
Neighborhood Gray Tone Difference Matrix (NGTDM)
Statistical Feature Matrix (SFM)
Law's Texture Energy Measures (LTE/TEM)
Fractal Dimension Texture Analysis (FDTA)
Gray Level Run Length Matrix (GLRLM)
Fourier Power Spectrum (FPS)
Shape Parameters
Gray Level Size Zone Matrix (GLSZM)
Higher Order Spectra (HOS)
Local Binary Pattern (LPB)
```
#### Morphological Features
```
Grayscale Morphological Analysis
Multilevel Binary Morphological Analysis
```
#### Histogram Based Features
```
Histogram
Multi-region histogram
Correlogram
```
#### Multi-scale Features
```
Fractal Dimension Texture Analysis (FDTA)
Amplitude Modulation – Frequency Modulation (AM-FM)
Discrete Wavelet Transform (DWT)
Stationary Wavelet Transform (SWT)
Wavelet Packets (WP)
Gabor Transform (GT)
```
#### Other Features
```
Zernikes' Moments
Hu's Moments
Threshold Adjacency Matrix (TAS)
Histogram of Oriented Gradients (HOG)
```
### [FeatureExtraction](https://github.com/ProfBressan/FeatureExtraction)

#### BIC (Border/Interior Pixel Classification) 128 features

#### TAS (Threshold Adjacency Statistics) 162 features

#### LBP (Local Binary Part) 352 features

#### FOM (First Order Measures) 8 features (gray) | 24 features (color)

#### Zernike 72 features

#### Haralick 13 features

#### GCH (Global Color Histogram) 30 features

### [Mahotas](https://mahotas.readthedocs.io/en/latest/features.html)

#### SURF

#### Zernike

#### Haralick

#### LBP

### [scMultiplex](https://github.com/fmi-basel/gliberal-scMultipleX/tree/main)
```
spacing
set\_spacing
fixed\_percentiles
skewness
kurtos
stdv
bounding\_box\_ratio
convex\_hull\_ratio
convex\_hull\_area\_resid
convex\_hull\_centroid\_dif
circularity
aspect\_ratio
minor\_major\_axis\_ratio
concavity\_count
disconnected\_component
surface\_area\_marchingcube
flag\_touching
is\_touching\_border\_xy
is\_touching\_border\_z
centroid\_weighted\_correct
centroid\_weighted\_correct
```

###[Ilastik](https://github.com/ilastik/ilastik/tree/main)

#### Standard Pixel Features
```
Color/Intensity: these features should be selected if the color or brightness can be used to discern objects
Edge: should be selected if brightness or color gradients can be used to discern objects.
Texture: this might be an important feature if the objects in the image have a special textural appearance.
```
#### Standard Object Features
```
Bounding Box Maximum
The coordinates of the upper right corner of the object's bounding box. The first axis is x, then y, then z (if available).
Bounding Box Minimum
The coordinates of the lower left corner of the object's bounding box. The first axis is x, then y, then z (if available).
Size in pixels
Total size of the object in pixels. No correction for anisotropic resolution or anything else.
Covariance of Channel Intensity
For multi-channel images this feature computes the covariance between the channels inside the object.
Covariance of Channel Intensity in neighborhood
For multi-channel images this feature computes the covariance between the channels in the object neighborhood. The size of the neighborhood is determined from the controls in the lower part of the dialogue.
Histogram of Intensity
Histogram of the intensity distribution inside the object. The histogram has 64 bins and its range is computed from the global minimum and maximum intensity values in the whole image.
Histogram of Intensity in neighborhood
Histogram of the intensity distribution in the object neighborhood. The histogram has 64 bins and its range is computed from the global minimum and maximum intensity values in the whole image. The size of the neighborhood is determined from the controls in the lower part of the dialogue.
Kurtosis of Intensity
Kurtosis of the intensity distribution inside the object, also known as the fourth standardized moment. This feature measures the heaviness of the tails for the distribution of intensity over the object's pixels. For multi-channel data, this feature is computed channel-wise.
Kurtosis of Intensity in neighborhood
Kurtosis of the intensity distribution in the object neighborhood, also known as the fourth standardized moment. This feature measures the heaviness of the tails for the distribution of intensity over the object's pixels. For multi-channel data, this feature is computed channel-wise. The size of the neighborhood is determined from the controls in the lower part of the dialogue.
Maximum intensity
Maximum intensity value inside the object. For multi-channel data, this feature is computed channel-wise.
Maximum intensity in neighborhood
Maximum intensity value in the object neighborhood. For multi-channel data, this feature is computed channel-wise. The size of the neighborhood is determined from the controls in the lower part of the dialogue.
Mean Intensity
Mean intensity inside the object. For multi-channel data, this feature is computed channel-wise.
Mean Intensity in neighborhood
Mean intensity in the object neighborhood. For multi-channel data, this feature is computed channel-wise. The size of the neighborhood is determined from the controls in the lower part of the dialogue.
Minimum intensity
Minimum intensity value inside the object. For multi-channel data, this feature is computed channel-wise.
Minimum intensity in neighborhood
Minimum intensity value in the object neighborhood. For multi-channel data, this feature is computed channel-wise. The size of the neighborhood is determined from the controls in the lower part of the dialogue.
PrincipalAxes
PrincipalAxes, stay tuned for more details
Quantiles of Intensity
Quantiles of the intensity distribution inside the object, in the following order: 0%, 10%, 25%, 50%, 75%, 90%, 100%.
Principal components of the object
Eigenvectors of the PCA on the coordinates of the object's pixels. Very roughly, this corresponds to the axes of an ellipse fit to the object. The axes are ordered starting from the one with the largest eigenvalue.
Center of the object
Average of the coordinates of this object's pixels.
Radii of the object
Eigenvalues of the PCA on the coordinates of the object's pixels. Very roughly, this corresponds to the radii of an ellipse fit to the object. The radii are ordered, with the largest value as first.
Skewness of Intensity
Skewness of the intensity distribution inside the object, also known as the third standardized moment. This feature measures the asymmetry of the intensity distribution inside the object. For multi-channel data, this feature is computed channel-wise.
Skewness of Intensity in neighborhood
Skewness of the intensity distribution in the object neighborhood, also known as the third standardized moment. This feature measures the asymmetry of the intensity distribution in the object neighborhood. For multi-channel data, this feature is computed channel-wise. The size of the neighborhood is determined from the controls in the lower part of the dialogue.
Total Intensity
Sum of intensity values for all the pixels inside the object. For multi-channel images, computed channel-wise.
Total Intensity in neighborhood
Sum of intensity values for all the pixels in the object neighborhood. For multi-channel images, computed channel-wise. The size of the neighborhood is determined from the controls in the lower part of the dialogue.
Variance of Intensity
Variance of the intensity distribution inside the object. For multi-channel data, this feature is computed channel-wise.
Variance of Intensity in neighborhood
Variance of the intensity distribution in the object neighborhood. For multi-channel data, this feature is computed channel-wise. The size of the neighborhood is determined from the controls in the lower part of the dialogue.
```
#### Convex Hull Features
```
Convexity
The ratio between the areas of the object and its convex hull (\<= 1)
Defect Center
Combined centroid of convexity defects, which are defined as areas of the convex hull, not covered by the original object.
Number of Defects
Total number of defects, i.e. number of connected components in the area of the convex hull, not covered by the original object
Mean Defect Displacement
Mean distance between the centroids of the original object and the centroids of the defects, weighted by defect area.
Kurtosis of Defect Area
Kurtosis (4th standardized moment, measure of tails' heaviness) of the distribution of the areas of convexity defects. Defects are defined as connected components in the area of the convex hull, not covered by the original object.
Mean Defect Area
Average of the areas of convexity defects. Defects are defined as connected components in the area of the convex hull, not covered by the original object.
Skewness of Defect Area
Skewness (3rd standardized moment, measure of asymmetry) of the distribution of the areas of convexity defects. Defects are defined as connected components in the area of the convex hull, not covered by the original object.
Variance of Defect Area
Variance of the distribution of areas of convexity defects. Defects are defined as connected components in the area of the convex hull, not covered by the original object.
Convex Hull Center
Centroid of the convex hull of this object. The axes order is x, y, z
Convex Hull Area
Area of the convex hull of this object
Object Center
Centroid of this object. The axes order is x, y, z
Object Area
Area of this object, computed from the interpixel contour (can be slightly larger than simple size of the object in pixels). This feature is used to compute convexity.
```

#### Skeleton Features
```
Average Branch Length
Average length of a branch in the skeleton
Number of Branches
Total number of branches in the skeleton of this object.
Diameter
The longest path between two endpoints on the skeleton.
Euclidean Diameter
The Euclidean distance between the endpoints (terminals) of the longest path on the skeleton
Number of Holes
The number of cycles in the skeleton (i.e. the number of cavities in the region)
Center of the Skeleton
The coordinates of the midpoint on the longest path between the endpoints of the skeleton.
Length of the Skeleton
Total length of the skeleton in pixels
```

### [Scikit-image](https://scikit-image.org/docs/dev/api/skimage.feature.html)

#### Features
```
blob\_dog
blob\_doh
blob\_log
canny
corner\_fast
corner\_foerstner
corner\_harris
corner\_kitchen\_rosenfeld
corner\_moravec
corner\_orientations
corner\_peaks
corner\_shi\_tomasi
corner\_subpix
daisy
draw\_haar\_like\_feature
draw\_multiblock\_lbp
fisher\_vector
graycomatrix
graycoprops
haar\_like\_feature
haar\_like\_feature\_coord
hessian\_matrix
hessian\_matrix\_det
hessian\_matrix\_eigvals
hog
learn\_gmm
local\_binary\_pattern
match\_descriptors
match\_template
multiblock\_lbp
multiscale\_basic\_features
peak\_local\_max
plot\_matches
structure\_tensor
structure\_tensor\_eigenvalues
BRIEF
CENSURE
Cascade
ORB
SIFT
```

[https://scikit-image.org/docs/dev/api/skimage.filters.html](https://scikit-image.org/docs/dev/api/skimage.filters.html)

#### Filters
```
apply\_hysteresis\_threshold
butterworth
correlate\_sparse
difference\_of\_gaussians
farid
farid\_h
farid\_v
filter\_forward
filter\_inverse
frangi
gabor
gabor\_kernel
gaussian
hessian
laplace
median
meijering
prewitt
prewitt\_h
prewitt\_v
rank\_order
roberts
roberts\_neg\_diag
roberts\_pos\_diag
sato
scharr
scharr\_h
scharr\_v
sobel
sobel\_h
sobel\_v
threshold\_isodata
threshold\_li
threshold\_local
threshold\_mean
threshold\_minimum
threshold\_multiotsu
threshold\_niblack
threshold\_otsu
threshold\_sauvola
threshold\_triangle
threshold\_yen
try\_all\_threshold
unsharp\_mask
wiener
window
LPIFilter2D
```

### [Xrayimage](https://github.com/vatsalsaglani/xrayimage_extractfeatures/tree/master)

#### Entropy features
```
Shannon's Entropy
Simple Entropy
```
#### GLCM Features
```
Correlation
Homogeneity
Energy
Contrast
All GLCM Features
```
#### Moments
```
24 Variant Moments
7 Hu Moments (opencv)
```
#### Region Properties
```
max\_area
eccentricity
euler\_number
solidity
perimeter
mean\_area
std\_area
thresh\_img
bb
bb\_area
centroid\_r
convex\_area\_r
coordinates\_r
eq\_diameter
extent\_r
filled\_area\_r
inertia\_tensor\_r
inertia\_tensor\_eigvals\_r
label\_r
local\_centroid\_r
maj\_ax\_len
min\_ax\_len
```

## Spatial Features

### [SquidPy](https://squidpy.readthedocs.io/en/stable/)

#### Graph features
```
spatial\_neighbors: Create a graph from spatial coordinates.
nhood\_enrichment : Compute neighborhood enrichment by permutation test.
co\_occurrence :Compute co-occurrence probability of clusters.
centrality\_scores : Compute centrality scores per cluster or cell type.
interaction\_matrix : Compute interaction matrix for clusters.
ripley:Calculate various Ripley's statistics for point processes.
ligrec :Perform the permutation test as described in [Efremova et al., 2020].
spatial\_autocorr : Calculate Global Autocorrelation Statistic (Moran's I or Geary's C).
sepal: Identify spatially variable genes with Sepal.
```
#### Image features
```
process: Process an image by applying a transformation.
segment: Segment an image.
calculate\_image\_features : Calculate image features for all observations in adata.
```
#### Spatial plots
```
spatial\_scatter : Plot spatial omics data with data overlayed on top.
spatial\_segment : Plot spatial omics data with segmentation masks on top.
nhood\_enrichment : Plot neighborhood enrichment.
centrality\_scores : Plot centrality scores.
interaction\_matrix : Plot cluster interaction matrix.
ligrec : Plot the result of a receptor-ligand permutation test.
ripley: Plot Ripley's statistics for each cluster.
co\_occurrence : Plot co-occurrence probability ratio for each cluster.
extract(adata[, obsm\_key, prefix]) : Create a temporary anndata.AnnData object for plotting.
var\_by\_distance : Plot a variable using a smooth regression line with increasing distance to an anchor point.
```
### [DypFISH](https://github.com/cbib/dypfish)

#### Spatial features
```
Cytoplasmic total counts
Peripheral distance map
Cell quantization
2D quantization: quadrants and per quadrant statistics
Fine-grained quantization
3D quantization: 3D quadrants and per quadrant statistics
Peripheral fraction and enrichment
Volume corrected noise measure
Cytoplasmic spread
MTOC polarity index
mRNA / protein distribution profile
Colocalization score
Degree of clustering (Ripley-K)
mRNA spatial distribution
```

### [Giotto](https://giottosuite.readthedocs.io/en/master/documentation.html#cell-neighborhood)
#### Cell neighborhood
```
cellProximityEnrichment : Calculate Cell-Cell Interaction Enrichment
cellProximityBarplot : Create Barplot from Cell-Cell Proximity Score
cellProximityHeatmap : Create Heatmap from Cell-Cell Proximity Score
cellProximityNetwork : Create Network from Cell-Cell Proximity Score
cellProximitySpatPlot : Visualize Cell-Cell Interactions (2D)
cellProximitySpatPlot3D : Visualize Cell-Cell Interactions (3D)
```
#### Cell cell 
```
findInteractionChangedGenes : Identify Cell-Cell Interaction Changed Genes (ICGs)
findICG : Identify Cell-Cell Interaction Changed Genes (ICGs)
findCellProximityGenes : Identify Cell-Cell Interaction Changed Genes (ICGs)
findCPG : Identify Cell-Cell Interaction Changed Genes (ICGs)
filterCellProximityGenes : Identify Cell-Cell Interaction Changed Genes (ICGs)
filterInteractionChangedGenes : Filter The Identified Cell-Cell Interaction Changed Genes (ICGs)
filterICG : Filter ICGs
filterCPG : Filter ICGs
combineInteractionChangedGenes : Combine ICG Scores (Pairwise)
combineICG : Combine ICG Scores (Pairwise)
combineCellProximityGenes : Combine ICG Scores (Pairwise)
combineCPG : Combine ICG Scores (Pairwise)
plotInteractionChangedGenes : Visualize ICGs via Barplot
plotICG : Visualize ICGs via Barplot
plotCellProximityGenes : Visualize Cell Proximity Gene Scores
plotCPG : Visualize Cell Proximity Gene Scores
plotCombineInteractionChangedGenes : Visualize Combined ICG Scores
plotCombineICG : Visualize Combined ICG Scores
plotCombineCellProximityGenes : Visualize Combined ICG Scores
plotCombineCPG : Visualize Combined ICG Scores
```
#### Cell cell communication
```
exprCellCellcom : Calculate Cell-Cell Communication Scores
spatCellCellcom : Calculate Spatial Cell-Cell Communication Scores
plotCCcomDotplot : Plot Ligand-Receptor Communication Scores
plotRankSpatvsExpr : Plot Comparison of Ligand-Receptor Rankings
plotRecovery : Plot Comparison of Ligand-Receptor Rankings
```
### [Big-FISH](https://big-fish.readthedocs.io/en/stable/classification/features.html)
```
features\_distance
features\_in\_out\_nucleus
features\_protrusion
features\_dispersion
features\_topography
features\_foci
features\_area
features\_centrosome
```
### [IMCDataAnalysis (R)](https://bodenmillergroup.github.io/IMCDataAnalysis/performing-spatial-analysis.html)

#### Spatial interaction graphs

#### Spatial community analysis

#### Cellular neighborhood analysis

#### Spatial context analysis

#### Patch detection

#### Interaction analysis

### [SpatialEpiApp (R)](https://github.com/Paula-Moraga/SpatialEpiApp)

#### Areal 
```
Spatial neighborhood matrices
Spatial autocorrelation
Bayesian spatial models
Disease risk modeling
Areal data issues
```
#### Geostatistical features
```
Geostatistical data
Spatial interpolation methods
Kriging
Model-based geostatistics
Methods assessment
```
#### Spatial point patterns
```
Spatial point patterns
The spatstat package
Spatial point processes and simulation
Complete spatial randomness
Intensity estimation
The K-function
Point process modeling
```

## Graph features
### [NetworkX](https://networkx.org/documentation/stable/index.html)
```
Approximations and Heuristics
Assortativity
Asteroidal
Bipartite
Boundary
Bridges
Centrality
Chains
Chordal
Clique
Clustering
Coloring
Communicability
Communities
Components
Connectivity
Cores
Covering
Cycles
Cuts
D-Separation
Directed Acyclic Graphs
Distance Measures
Distance-Regular Graphs
Dominance
Dominating Sets
Efficiency
Eulerian
Flows
Graph Hashing
Graphical degree sequence
Hierarchy
Hybrid
Isolates
Isomorphism
Link Analysis
Link Prediction
Lowest Common Ancestor
Matching
Minors
Maximal independent set
non-randomness
Moral
Node Classification
Operators
Planarity
Planar Drawing
Graph Polynomials
Reciprocity
Regular
Rich Club
Shortest Paths
Similarity Measures
Simple Paths
Small-world
s metric
Sparsifiers
Structural holes
Summarization
Swap
Threshold Graphs
Time dependent
Tournament
Traversal
Tree
Triads
Vitality
Voronoi cells
Walks
Wiener index
```

## Other resources
### Ontology
[EDAM-BIOIMAGING](https://bioportal.bioontology.org/ontologies/EDAM-BIOIMAGING/?p=classes&conceptid=root)

General computer vision feature book :
[https://books.google.ch/books?hl=fr&lr=&id=KcW-DwAAQBAJ&oi=fnd&pg=PP1&dq=graph+features+extraction+python+from+image&ots=10qF2rQE3S&sig=YWUElA38FMI8bvFp419Uh9-CyrI#v=onepage&q&f=false](https://books.google.ch/books?hl=fr&lr=&id=KcW-DwAAQBAJ&oi=fnd&pg=PP1&dq=graph+features+extraction+python+from+image&ots=10qF2rQE3S&sig=YWUElA38FMI8bvFp419Uh9-CyrI#v=onepage&q&f=false)

Spatial data features:
[https://geographicdata.science/book/notebooks/12\_feature\_engineering.html](https://geographicdata.science/book/notebooks/12_feature_engineering.html)

--- 

TODO  : Attempt to organize features in categories

1. Color Features:
- Color histograms
- Color moments (mean, variance, skewness, kurtosis)
- Color space transformations (RGB, HSV, Lab, etc.)
2. Texture Features:
- Gray-level co-occurrence matrix (GLCM)
- Local binary patterns (LBP)
- Gabor filters
- Haralick texture features
3. Shape Features:
- Contour-based features
- Area, perimeter, and circularity
- Aspect ratio
- Euler number
4. Histogram Features:
- Histogram of pixel intensity values
- Histogram of gradient orientations (HOG)
5. Statistical Features:
- Mean, median, standard deviation
- Skewness and kurtosis
- Moments
- Percentile values
6. Edge Features:
- Canny edge detection features
- Sobel gradient features
- Laplacian of Gaussian (LoG) features
7. Corner Features:
- Harris corners
- Shi-Tomasi corners
8. Scale-Invariant Features:
- Scale-invariant feature transform (SIFT)
- Speeded-Up Robust Features (SURF)
- Oriented FAST and Rotated BRIEF (ORB)
9. Deep Learning Features:
- Features extracted from pre-trained convolutional neural networks (CNNs)
- Activation maps from intermediate layers
- Feature vectors from neural network embeddings
10. Segmentation Features:
- Region properties (area, centroid, etc.)
- Boundary-based features
- Texture or color features within segmented regions
11. Geometric Features:
- Geometric moments
- Fourier descriptors
- Zernike moments
12. Statistical Texture Features:
- Gray-level run-length matrix (GLRLM)
- Gray-level size zone matrix (GLSZM)
13. Fractal Features:
- Fractal dimension
- Lacunarity
14. Wavelet Transform Features:
- Wavelet coefficients
- Wavelet energy
15. Local Features:
- Histogram of oriented gradients (HOG)
- Local binary patterns (LBP)
- Haar-like features
16. Graph-based Features:
- Graph properties of image regions or objects
- Connectivity and graph invariants
17. Harmonic Features:
- Fourier coefficients
- Discrete cosine transform (DCT) coefficients
18. Statistical Moments:
- Central moments
- Normalized moments
- Hu moments
19. Morphological Features:
- Erosion and dilation statistics
- Skeletonization-based features
- Convexity and solidity
20. Spatial Domain Features:
- Spatial moments
- Autocorrelation


TODO : Biology oriented classification


1. Cell Morphology Features:
- Cell area
- Cell perimeter
- Cell circularity
- Cell eccentricity
- Cell aspect ratio
2. Nucleus Features:
- Nucleus area
- Nucleus perimeter
- Nucleus circularity
- Nucleus eccentricity
- Nucleus aspect ratio
3. Texture Features for Tissue Analysis:
- Gray-level co-occurrence matrix (GLCM) for texture
- Haralick texture features
- Run-length texture features
- Gabor texture features
4. Vascular Features:
- Vessel diameter
- Vessel tortuosity
- Vessel branching patterns
- Vessel density
