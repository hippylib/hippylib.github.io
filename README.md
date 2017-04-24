# hIPPYlib / web

This repo contains the [MkDocs](http://mkdocs.org) source files of the hIPPYlib website.

## Prerequisites

* Install [MkDocs](http://mkdocs.org)   
```
#!sh
    sudo pip install mkdocs==0.14.0
```

* Install [bootstrap](http://getbootstrap.com/) and [bootswatch](https://bootswatch.com/) themes
```
#!sh
    sudo pip install mkdocs-bootstrap==0.1.1
    sudo pip install mkdocs-bootswatch==0.1.0
```
    
* Install [MathJax](https://www.mathjax.org/) support using [python-markdown-math](https://github.com/mitya57/python-markdown-math)
```
#!sh
    sudo pip install python-markdown-math==0.2
```


## To make changes to the website

* Clone this repository
```
#!sh
    git clone git@bitbucket.org:hippylibdev/hippylib-website.git
```
* Add the github repository as public (i.e. where we will publish the website).
```
#!sh
    git remote add public https://github.com/hippylib/hippylib.github.io.git
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
    # push this repository to github
```