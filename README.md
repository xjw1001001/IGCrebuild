IGCexpansion
=======

Gene Conversion Project

This is the first approach with the smallest possible multi-gene family (only two paralogs). Code for IGC inference as described in Ji-Griffing-Thorne 2016 MBE paper is available here.

Data: Yeast

Dependent: 

[jsonctmctree package](http://jsonctmctree.readthedocs.org/en/latest/) (powerful likelihood  calculation 
engine by Alex Griffing)

[Biopython](http://biopython.org/wiki/Biopython) (you could install it by `pip install --user Biopython`)

Coding Language: Python 2.7

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

`
pip install --user Biopython
`


4, use Run.py to perform analyses with following example commands:

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
pip install --user --upgrade numpy scipy
`
If the problem persists, you may solve it by either updating your LAPCK package to [atlas](http://math-atlas.sourceforge.net/) or installing an older version of jsonctmctree by the following code:

	`pip uninstall jsonctmctree`

	`pip install --user git+https://github.com/argriffing/jsonctmctree.git@cb1ba60ee2b57d6703cd9a3987000c2fd4dd68a5`

2. If you encounter "NotImplementedError: Wrong number or type of arguments for overloaded function 'coo_matvec'". Please update your scipy and numpy packages to the newest versions.

3. If you are using windows operating system, and encountered "numpy.distutils.system_info.NotFoundError: no lapack/blas resources found" when building/upgrading scipy package. Please try Enthought Canopy or Anaconda as suggested [Here](http://docs.scipy.org/doc/numpy-1.10.1/user/install.html#id4). They come with built in packages like numpy and scipy. 

4. If you encountered "Fatal error in launcher: Unable to create process using ..." when using pip in cmd (windows os). Please try ` python -m pip install xxx` as suggested [here](http://stackoverflow.com/questions/24627525/fatal-error-in-launcher-unable-to-create-process-using-c-program-files-x86).



Contact me at:
xji3 _at_ ncsu.edu