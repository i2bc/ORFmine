## Installation


### 1. Download and uncompress the latest release archive

#### Download the latest release
Here: 
[ ![](img/icons/download_16x16.png "Click to download the latest release")](https://github.com/nchenche/orfmap/releases/latest/)

#### Uncompress the archive
If you downloaded:

* the *.zip* file: ```unzip orfmap-x.x.x.zip```
* the *.tar.gz* file: ```tar xzvf orfmap-x.x.x.tar.gz```


### 2. Create an isolated environment
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

Once activated, any python library you'll install using pip will be installed solely in this isolated environment.
Every time you'll need to work with libraries installed in this environment, you'll have
to activate it. 

Once you're done working on your project, simply type `deactivate` to exit the environment.


<div class="admonition note">
    <p class="first admonition-title">
        Note
    </p>
    <p class="last">
        If you're sure you won't need your isolated environment anymore, you could clear it simply
        by removing the directory:
        <code>rm -r my_env/</code>
    </p>
</div>



### 3. Install ORFMap in your isolated environment

Be sure you're virtual environment is activated, and then follow the procedure described below.

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
