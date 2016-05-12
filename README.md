IGCexpansion
=======

Gene Conversion Project

This is the first approach with the smallest possible multi-gene family (only two paralogs). 

Data: Yeast

Dependent: 

[jsonctmctree package](http://jsonctmctree.readthedocs.org/en/latest/) (powerful likelihood  calculation 
engine by Alex Griffing)

[Biopython](http://biopython.org/wiki/Biopython) (you could install it by `pip install --user Biopython`)

Coding Language: Python

Instruction of how to run this code:

0, To install python packages, you need to use [pip](https://pip.pypa.io/en/stable/installing/) (package management). 

1, Install jsonctmctree package:

`
pip install --user git+https://github.com/argriffing/jsonctmctree.git
`

2, Install IGCexpansion package:

`
git clone https://github.com/xji3/Genconv.git
`

`
cd Genconv/IGCexpansion
`

`
pip install --user .
`  _(preferred)_

or

`
python setup.py install
`  


3, Similarly install any other python packages

`
pip install --user networkx
`


4, use Run.py or RunBootstrap.py to perform analyses with following example commands:

`
python Run.py --model MG94 --paralog1 YBL087C --paralog2 YER117W --no-force --no-clock
`


To uninstall:

`
pip uninstall IGCexpansion
`

Some known issues:

1. If you encounter "ValueError: LAPACK function dlange could not be found" when running the code. First, update your scipy packages, it's very likely to be out-dated. You can do this by
`
pip install --user --upgrade scipy numpy
`
If the problem persists, you may solve it by either updating your LAPCK package to [atlas](http://math-atlas.sourceforge.net/) or installing an older version of jsonctmctree by the following code:

	`pip uninstall jsonctmctree`

	`pip install --user git+https://github.com/argriffing/jsonctmctree.git@cb1ba60ee2b57d6703cd9a3987000c2fd4dd68a5`

2. IF you encounter "NotImplementedError: Wrong number or type of arguments for overloaded function 'coo_matvec'". Please update your scipy and numpy packages to the newest versions.



Contact me at:
xji3 _at_ ncsu.edu