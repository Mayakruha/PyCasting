# Installation
```
pip install git+https://github.com/mayakruha/PyCasting.git
```
#
Simple quasi 2D model of a continuous casting billet for calculations of shrinkage, solid shell thickness, stress/temperature distribution. Elastic strains are neglected in the model. Stresses, strains etc are calculated on basis of a creep law and the suggestion that cross-sections stay being plane.
The model fits well for calculations in the upper part of a continuous casting machine where solid shell thickness is small relative to billet cross-section size. The model was used to calculate mold taper of Continuous Casting Machine for a number of metallurgy plants in Russia including JSC Severstal, NLMK.
The model was described in the article: A.A.Makrushin, A.V.Kuklev., Yu.M.Ayzin "Computation of a surface profile for a narrow slab side in the area of a mould" - Stal - Apr 2004 (http://www.corad.ru/pdf/Steel_2004_04.pdf in russian)
