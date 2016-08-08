# Installing hIPPYlib

The current version of hIPPYlib depends on [FEniCS](http://fenicsproject.org/) version 1.6.
Below you can find how to install the correct version of FEniCS on your system.

- **MacOS 10.10 and 10.11 systems:**

Download FEniCS 1.6.0 from [here](https://fenicsproject.org/download/osx_details.html).
Find  your  MacOS  version  (either  10.11  or  10.10)  and  download  the appropriate  binaries  of  FEniCS  1.6.0.
If  you  are  running bash as default shell, you can add the following line to your profile in your home folder:
```sh
source /Applications/FEniCS.app/Contents/Resources/share/fenics/fenics.conf
```

Alternatively you can just double-click on the FEniCS icon in your Applications directory and that will generate a new shell preconfigured with the paths that FEniCS needs. Just run FEniCS from within this shell.

- **MacOS 10.9:**

Download FEniCS 1.5.0 from [here](http://fenicsproject.org/download/older_releases.html#older-releases).
Find  your  MacOS  version and  download  the appropriate  binaries  of  FEniCS  1.5.0.
Note FEniCS 1.5.0 will not be supported by future releases of hIPPYlib

- **Ubuntu LTS 14.04:**

Open a shell and run the following commands:

```sh
sudo add-apt-repository ppa:fenics-packages/fenics-1.6.x
sudo apt-get update
sudo apt-get install -y fenics
sudo apt-get dist-upgrade
sudo apt-get install -y ipython-notebook
sudo apt-get install -y paraview
sudo apt-get install -y git
```
  
If in the future you decide to uninstall FEniCS and remove all its dependencies, you can run the following commands:
```sh
sudo apt-get purge --auto-remove fenics
sudo ppa-purge ppa:fenics-packages/fenics-1.6.x
```
    
