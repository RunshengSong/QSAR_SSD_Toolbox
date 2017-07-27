# QSAR-based SSD ToolBox

Runsheng Song </br>
runsheng@umail.ucsb.edu

---
A framework to create Species [Sensitivity Distributions (SSD)](https://www3.epa.gov/caddis/da_advanced_2.html) using pre-trained QSAR models

QSAR models were developed using Neural Networks in Tensorflow + Keras. 
Descriptors were calculated using [Rdkit](https://github.com/rdkit/rdkit) and [Mordred](https://github.com/mordred-descriptor/mordred), and optimized using tree-based feature selection.

All QSAR models have been cross-validated

## Prerequisite:
* Anaconda Python 2.7
* Recommend using Linux or MacOS. 

## Install 

* download the package
```bash
git clone https://github.com/RunshengSong/QSAR_SSD_Toolbox
```

* create a new conda environment

```bash
conda env create -f qsar.yml
```

```bash
python setup.py install
```

## Basic Usage
Input: List of Chemical SMILEs
Output: List of LC50 in -log10(LC50), mol/L

```python
from QSAR_SSD_Toolbox import qsar

this_q = qsar(model_name='Lepomis_Macrochirus')
prediction = this_q.predict(['CCC'])
```
---
## Available QSAR Models:
* Lepomis_Macrochirus:
R^2 on testing chemicals: 0.51 </br>
![image1]('/QSAR_SSD_Toolbox/models/Lepomis_Macrochirus/0714a_results.png')
</br>

* Oncorhynchus_Mykiss:
R^2 on testing chemicals: 0.72 </br>
![image2]('/QSAR_SSD_Toolbox/models/Oncorhynchus_Mykiss/0714a_results.png')
</br>