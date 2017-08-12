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

* Install rdkit with conda first(save ur life):
```bash
conda install -c rdkit rdkit=2017.03.1
```

* Install QSAR_SSD_Toolbox via pip:
```bash
pip install QSAR_SSD_Toolbox
```

* Install the requirments.txt if some packages are missing via 
```bash
pip install -r requirements.txt
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
from scipy.stats import lognorm

this_ssd = ssd_generator()
this_ssd.generate(this_q, dist=lognorm, run_bootstrap=True, bootstrap_time=1000, display_range=[0.8,100]) # this will return a plot with bootstrap and baseline SSD curves. For more information about bootstrap in SSD refer to this blog: https://edild.github.io/ssd/
```

---
## Available QSAR Models:
* Lepomis Macrochirus:
R^2 on testing chemicals: 0.51 </br>
![](QSAR_SSD_Toolbox/models/Lepomis&#32;Macrochirus/0714a_results.png?raw=true)
</br>

* Oncorhynchus Mykiss:
R^2 on testing chemicals: 0.72 </br>
![](QSAR_SSD_Toolbox/models/Oncorhynchus&#32;Mykiss/0713a_results.png?raw=true)
</br>

* Americamysis bahia:
R^2 on testing chemicals: 0.45 </br>
![](QSAR_SSD_Toolbox/models/Americamysis&#32;bahia/results.png?raw=true)
</br>

* Oncorhynchus Mykiss:
R^2 on testing chemicals: 0.75 </br>
![](QSAR_SSD_Toolbox/models/Oncorhynchus&#32;Mykiss/0713a_results.png?raw=true)
</br>

* Oryzias latipes:
R^2 on testing chemicals: 0.56 </br>
![](QSAR_SSD_Toolbox/models/Oryzias&#32;latipes/results.png?raw=true)
</br>

* Pimephales promelas:
R^2 on testing chemicals: 0.72 </br>
![](QSAR_SSD_Toolbox/models/Pimephales&#32;promelas/results.png?raw=true)
</br>

* Daphnia magna:
R^2 on testing chemicals: 0.77 </br>
![](QSAR_SSD_Toolbox/models/Daphnia&#32;magna/results.png?raw=true)
</br>

* Other water fleas model:
This model include the experimental data (LC50) of different kind of water fleas, except Daphnia magna
R^2 on testing chemicals: 0.61 </br>
![](QSAR_SSD_Toolbox/models/Other&#32;Water&#32;fleas/results.png?raw=true)
</br>
