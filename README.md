# Next-Generation-File-Format-Features

This GitHub repository dedicated to the discussion around new specifications for feature definition,  extraction and storage within the context of the OME NGFF (Next Generation File Format) standard. This collaborative work is the the combined efforts of @romainGuiet, @fdsteffen, @DillanSaunders, @Tpoquillon, @retogerber, @ebouilho, @nrepina  who participated in the Zurich Next Generation Image Analysis Workflow Hackathon 2023.

Key objectives:

* Standardization: Establish standardized guidelines for defining and storing features, ensuring interoperability across various image analysis tools and platforms.

* Flexibility: Create a framework that accommodates diverse feature types and differents bioimaging context, empowering users to extract and store a wide range of biological and image-based data.

* Efficiency: Optimize feature storage to minimize data redundancy and enhance computational efficiency.

* Community Engagement: Encourage collaboration and feedback from the wider scientific and bioimaging community to refine and improve these specifications continually.

## Multi-Object Feature Relationship Managment

*Biological images frequently encompass various object types and sizes, 
(such as organoids, cell mitochondria, and p-bodies). The process of 
quantifying these objects and deriving relevant features from them also 
encompasses the extraction of descriptors that characterize their 
interactions. This gives rise to challenges in managing the features of 
multiple objects, elucidating, quantifying, and storing their 
interrelationships.*

These issues are curently discussed in:
- https://github.com/Tpoquillon/Next-Generation-File-Format-Features/pull/2 a proposal for tracking data uscase
- https://github.com/Tpoquillon/Next-Generation-File-Format-Features/pull/1 a proposal for a multilevel object uscase


## Feature Atlas
*Many different bioimaging pipeline implement feature measurement, we observe that there is  both a strong redundancy between those, with the same features being re-implemented and  a need to implement more, a challenge to categorize those features and a  need for new meaningful feature descriptor*

- [A compilation of different feature library from @ebouilho](Feature_libraries.md)
- https://github.com/Tpoquillon/Next-Generation-File-Format-Features/pull/3 is a proposal for a unified feature table draft


## Related works

- Nyxus Features https://github.com/PolusAI/nyxus/blob/main/docs/source/featurelist.rst

- scMultipleX https://github.com/fmi-basel/gliberal-scMultipleX/blob/main/src/scmultiplex/features/FeatureFunctions.py

- AnnData [anndata - Annotated data &#8212; anndata 0.11.0.dev24+gaf7a5b7 documentation](https://anndata.readthedocs.io/en/latest/)

- Parquet httpTable n-data storage format for OME-NGFF https://github.com/ome/ngff/pull/64


