# QSAR-based SSD ToolBox

Runsheng Song </br>
runsheng@umail.ucsb.edu

---
A framework to create [Species Sensitivity Distributions (SSD)](https://www3.epa.gov/caddis/da_advanced_2.html) using pre-trained QSAR models

QSAR models were developed using Neural Networks in Tensorflow + Keras. 
Descriptors were calculated using [Rdkit](https://github.com/rdkit/rdkit) and [Mordred](https://github.com/mordred-descriptor/mordred), and optimized using tree-based feature selection.

All QSAR models have been cross-validated

Current toxicity endpoint is [LC50](http://www.businessdictionary.com/definition/lethal-concentration-50-LC50.html)

## Prerequisite:
* Anaconda Python 2.7
* Recommend using Linux or MacOS. 

## Install 

* download the package
```bash
git clone https://github.com/RunshengSong/QSAR_SSD_Toolbox
```
* Install rdkit with conda first(save ur life):
```bash
conda install -c rdkit rdkit=2017.03.1
```

* Install other requirements (this will install packages like tensorflow and keras):
```bash
pip install -r requirements.txt
```
* Install QSAR_SSD_Toolbox with setup.py:
```bash
python setup.py install
```

## Basic Usage
### Single QSAR Model on One Species:
```python
from QSAR_SSD_Toolbox.src.qsar import qsar

SMILEs = ['CCCC'] # The input SMILEs must be a list

this_q = qsar("Lepomis Macrochirus") # the name of the species, see below for avaliable species
print this_q.predict(SMILEs) # return a list of predicted LC50 values for the given species
```

### Run model on all species and prepare for SSD development:
```python
from QSAR_SSD_Toolbox.src.qsar import run_all

SMILEs = ['CCCC'] # The input SMILEs must be a list

this_q = run_all.run(SMILEs) # return a pandas dataframe: index -- species name | 'val' -- LC50 values for the input chemicals on corrosponding species. 
```

### Plot SSD Curves
```python
from QSAR_SSD_Toolbox.src.ssd import ssd_generator

this_ssd = ssd_generator()
this_ssd.generate(this_q, dist=lognorm, run_bootstrap=True, bootstrap_time=1000, display=True) # this will return a plot with bootstrap and baseline SSD curves. For more information about bootstrap in SSD refer to this blog: https://edild.github.io/ssd/
```
---
## Available QSAR Models:
* Lepomis_Macrochirus:
R^2 on testing chemicals: 0.51 </br>
![image1](QSAR_SSD_Toolbox/models/Lepomis Macrochirus/0714a_results.png?raw=true)
</br>

* Oncorhynchus_Mykiss:
R^2 on testing chemicals: 0.72 </br>
![image2](QSAR_SSD_Toolbox/models/Oncorhynchus Mykiss/0713a_results.png?raw=true)
</br>
