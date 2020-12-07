## Installation


### 1. Overview

ORFmine is a package that consists of two independent programs ORFmap
and ORFold. The setup.py program will install both ORFmap and ORFold. 
They can be used together or independently. 

If you are not interested in the calculation of the disorder 
and/or aggregation propensities with ORFold, you do not need to 
install further dependencies and simply have to follow the 
installation instructions presented here (xxxx). 
 
Otherwise, you need first to download IUPRED[REF] and TANGO[REF] at 
xxx and yyyy respectively.

Then xxxxx 

Then follow the installation instructions presented here (xxxx)

### 2. Download and uncompress the latest release archive

#### Download the latest release
Here: 
[ ![](img/icons/download_16x16.png "Click to download the latest release")](https://github.com/nchenche/orfmap/releases/latest/)

#### Uncompress the archive
If you downloaded:

* the *.zip* file: ```unzip orfmine-x.x.x.zip```
* the *.tar.gz* file: ```tar xzvf orfmine-x.x.x.tar.gz```


### 3. Create an isolated environment
Although not strictly necessary, this step is highly recommended (it will allow you to work on different projects without having
any conflicting library versions).
 
#### Install virtualenv
``` python
python3 -m pip install virtualenv
```

#### Create a virtual environment
```bash
virtualenv -p python3 my_env
```

#### Activate the created environment
```bash
source my_env/bin/activate
```

Once activated, any python library you will install using pip 
will be installed solely in this isolated environment.
You must activate this environment very time you need libraries installed 
in this environment. 

Once you are done working on your project, 
simply type `deactivate` to exit the environment.


<div class="admonition note">
    <p class="first admonition-title">
        Note
    </p>
    <p class="last">
        To delete definitively your virutal environment, you can simply
        remove the directory with the following instruction:
        <code>rm -r my_env/</code>
    </p>
</div>



### 4. Install ORFMine in your isolated environment

Be sure your virtual environment is activated, 
and then follow the procedure described below.

#### Go to the ORFMap directory
 
```bash
cd orfmap-x.x.x/
```

#### Install 

```python
python setup.py install
```

or 
```python
pip install .
```


<div class="admonition note">
    <p class="first admonition-title">
        Note
    </p>
    <p class="last">
        If you need to remove the package (for example in order to install a newest version), try:
        <br> <code>pip uninstall orfmap</code> </br>
    </p>
</div>


