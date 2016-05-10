Genconv
=======

Gene Conversion Project

This is the first approach with the smallest possible multi-gene family (only two paralogs). 

Data: Yeast

Dependent: [jsonctmctree package](http://jsonctmctree.readthedocs.org/en/latest/) (powerful likelihood  calculation 
engine by Alex Griffing)

Coding Language: Python

Instruction of how to run this code:

1, Install jsonctmctree package:

`
pip install --user git+https://github.com/argriffing/jsonctmctree.git
`

2, Install IGCexpansion package:

`
cd IGCexpansion
`

`
pip install .
`  _(preferred)_

or

`
python setup.py install
`  


3, Similarly install any other python packages


4, use Run.py or RunBootstrap.py to perform analyses with following example commands:

`
python Run.py --model MG94 --paralog1 YBL087C --paralog2 YER117W --no-force --no-clock
`

or


`
python RunBootstrap.py  --paralog1 YBL087C --paralog2 YER117W --bootnum 1
`

To uninstall:

`
pip uninstall IGCexpansion
`

Contact me at:
xji3 _at_ ncsu.edu