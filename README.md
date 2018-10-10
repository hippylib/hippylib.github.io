# hIPPYlib / web

This repo contains the [MkDocs](http://mkdocs.org) source files of the hIPPYlib website.

## Prerequisites

* Install [MkDocs](http://mkdocs.org)   
```
#!sh
    sudo pip install mkdocs==0.17.2
```

* Install [bootstrap](http://getbootstrap.com/) and [bootswatch](https://bootswatch.com/) themes
```
#!sh
    sudo pip install mkdocs-bootstrap==0.2.0
    sudo pip install mkdocs-bootswatch==0.5.0
```
    
* Install [MathJax](https://www.mathjax.org/) support using [python-markdown-math](https://github.com/mitya57/python-markdown-math)
```
#!sh
    sudo pip install python-markdown-math==0.6
```

Or, if you have conda installed, simply type

```conda-env create --file mkdocs-env.txt```



## To make changes to the website

* Clone this repository

```
#!sh
    git clone git@github.com:hippylib/web.git
```

* Edit or add some `.md` files in the src folder (you may also need to update the `mkdocs.yml` config); 
* Preview locally
```
#!sh
    mkdocs serve
```
* Publish on GitHub
```
#!sh
    git push origin master:master
    mkdocs build
    cd site
    # push this repository to github: https://github.com/hippylib/hippylib.github.io.git
```
